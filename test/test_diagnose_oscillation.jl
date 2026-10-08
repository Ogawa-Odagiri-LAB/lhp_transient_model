# =============================================================================
# test_diagnose_oscillation.jl
#
# 液管の速度の振動の原因を切り分ける診断。
#
# 観察 (transient.jl の 200 s の結果、saveat=1 s のプロット):
#   ・液管の平均速度が振動する (振幅 0.05〜0.2 m/s。定常の流速は 0.0076 m/s (20 W))。
#     グルーブ・蒸気管・凝縮器の速度は滑らか
#   ・液管の P−P_cc は、過渡で ±500〜1500 Pa、後半は ±50〜100 Pa のノイズ
#   ・saveat=1 s では、数 Hz の振動は、折り返しで、ランダムなノイズに見える
#
# 候補:
#   H1 許容誤差レベルの雑音が、液体の硬さ (∂P/∂ρ=1.1e6 Pa/(kg/m³), ∂P/∂T=1.7e6 Pa/K)
#      で増幅される。 → 振幅が rtol に比例する。広帯域
#   H2 液柱の物理的な振動 (液管の液柱 = 質量、蒸気側の圧縮性 = ばね。約 1.4 Hz の見積もり)
#      → 狭い周波数のピーク。振幅は rtol に依存しない。凝縮の前線の位置と相関
#   H3 離散的な事象 (凝縮相関の領域切り替え、テーブルの格子線の通過、Newton 失敗) が
#      液柱を叩く → 孤立したスパイク
#
# このスクリプトは、実運転の幾何で、指定した時間窓の中の IDA の全ステップを、コールバックで
# 記録し (全ステップの保存はメモリが膨大)、周波数スペクトルを求める。設定 (rtol, 次数) を変えて
# 繰り返し、振幅が、設定でどう変わるかを見る。
#
# 実行時間の目安: FluidTables の型安定化 (残差評価が約 7 倍速くなった) 後は、1 設定あたり、
# T_RUN=70 s で 1〜2 分、既定の 4 設定で 5〜10 分 (型安定化の前は、1 設定 10 分前後だった)。
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays,
      SparseArrays, LinearAlgebra, Random, Sundials, Statistics
import DiffEqBase      # (`using DiffEqBase: 名前...` は、列挙した名前だけを取り込み、モジュール名 DiffEqBase は束縛しない。
                       #  DiffEqBase.u_modified! のように、モジュール名を直接使うので、import も必要)
using DiffEqBase: DAEProblem, DAEFunction, solve, DiscreteCallback
# u_modified! は、SciMLBase の関数 (DiffEqBase が再エクスポートする)。コールバックで使うので、先に確認する
isdefined(DiffEqBase, :u_modified!) ||
    error("DiffEqBase.u_modified! が見つかりません。SciMLBase を、プロジェクトに追加して、" *
          "`using SciMLBase` し、コールバックの DiffEqBase.u_modified! を SciMLBase.u_modified! に置き換えてください")

# ---------------------------- 設定 ----------------------------
const T_RUN = 70.0                       # 計算する時間 [s] (Q のステップアップは t=50 s)
# 記録・解析する時間窓 [s]: ステップ前の準定常 (40〜44 s) と、ステップ後の過渡 (60〜64 s。振動が大きい)
const WINDOWS = [(40.0, 44.0), (60.0, 64.0)]
const RESAMPLE_DT = 0.002                # 再標本化の間隔 [s] (周波数の上限 = 250 Hz)
const FMAX = 100.0                       # スペクトルを調べる最大の周波数 [Hz]
# (名前, rtol, max_order)
const CONFIGS = [
    ("既定: rtol=1e-4, 次数2",  1.0e-4, 2),
    ("rtol=1e-5, 次数2",        1.0e-5, 2),
    ("rtol=1e-6, 次数2",        1.0e-6, 2),
    ("rtol=1e-4, 次数1",        1.0e-4, 1),
]
# --------------------------------------------------------------

const SRC = joinpath(@__DIR__, "..", "src")
for f in ("physics", "smoothing", "types", "proptables", "propaccess", "correlations",
          "heattransfer", "geometry_cc", "residuals_transport", "residuals_wall",
          "residuals_wick", "residuals_cc", "residuals_evapcase", "residuals_grv",
          "residuals_main", "initialcondition", "sparsity", "solverfuncs")
    include(joinpath(SRC, f * ".jl"))
end

fluid = "ammonia"
tables = init_fluid_lookups(fluid; T_min=230.0, T_max=360.0, P_max=4.5e6)
geo = LHPGeometry(
    WickGeometry(N=3, φ=0.65, ρ_wick=8900.0, Cp_wick=444.0),
    GrvGeometry(N=1),
    TLGeometry(N=20, L=2, D_i=0.00435, D_o=0.00635),
    CondGeometry(N=50, L=1.0, D_i=0.00435, D_o=0.00635),
    TLGeometry(N=20, L=2, D_i=0.00218, D_o=0.00318),
    CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037),
)
condition = LHPTransientCondition(
    fluid = fluid, T_cond = 15.0 + 273.15, T_amb = 25.0 + 273.15, T_init = 25.0 + 273.15,
    Q_apply_init = 20.0, Q_apply = 40.0, t_step = 50.0, τ_step = 1.0, ΔT_startup = 0.0,
    t_init = 0.0, t_end = T_RUN, saveat = T_RUN)

axes_ = getaxes(initialize_state(geo))     # 状態ベクトルの軸 (コールバックで、ComponentArray に包むのに使う)

# --- 窓の中の全ステップで記録する量 ---
const REC_NAMES = (:t, :u_ll_mean, :u_ll_first, :u_ll_last, :u_cond_mean, :u_vl_mean,
                   :dP_ll_mean, :dP_cond_mean, :n_vap, :T_ll_first)
function make_recorder(windows)
    rec = Dict(n => Float64[] for n in REC_NAMES)
    inwin(t) = any(w -> w[1] <= t <= w[2], windows)
    function record!(integrator)
        c = ComponentArray(copy(integrator.u), axes_)    # Sundials の内部バッファを直接触らないよう、コピーする
        P_cc = c.cc.P[1]
        push!(rec[:t], integrator.t)
        push!(rec[:u_ll_mean], mean(c.ll.u)); push!(rec[:u_ll_first], c.ll.u[1]); push!(rec[:u_ll_last], c.ll.u[end])
        push!(rec[:u_cond_mean], mean(c.cond.u)); push!(rec[:u_vl_mean], mean(c.vl.u))
        push!(rec[:dP_ll_mean], mean(c.ll.P) - P_cc); push!(rec[:dP_cond_mean], mean(c.cond.P) - P_cc)
        push!(rec[:n_vap], count(>(0.5), c.cond.x)); push!(rec[:T_ll_first], c.ll.T[1])
        # u を変えないので、積分器に、そう伝える。伝えないと、積分器は、IDA を再初期化し、
        # CheckInit で失敗する ("Initial conditions do not satisfy the DAE constraints")
        DiffEqBase.u_modified!(integrator, false)
        return nothing
    end
    return DiscreteCallback((u, t, integrator) -> inwin(t), record!; save_positions = (false, false)), rec
end

# --- 解析: 窓の中の時系列を、等間隔に再標本化して、スペクトルを求める ---
function interp1(ts, ys, t)
    k = searchsortedlast(ts, t)
    k >= length(ts) && return ys[end]
    k < 1 && return ys[1]
    w = (t - ts[k]) / (ts[k+1] - ts[k])
    return (1 - w) * ys[k] + w * ys[k+1]
end
function spectrum_peaks(tt, yy, w; npeaks=4)
    ts = collect(range(w[1], w[2], step = RESAMPLE_DT))
    y = [interp1(tt, yy, t) for t in ts]
    y .-= mean(y)
    freqs = collect(range(0.25, FMAX, step = 0.25))
    N = length(ts)
    amp = [abs(sum(y[j] * cis(-2π * f * (ts[j] - ts[1])) for j in 1:N)) * 2 / N for f in freqs]
    peaks = Tuple{Float64,Float64}[]
    for k in 2:length(freqs)-1
        (amp[k] > amp[k-1] && amp[k] >= amp[k+1]) && push!(peaks, (freqs[k], amp[k]))
    end
    sort!(peaks; by = x -> -x[2])
    return first(peaks, npeaks)
end

results = Dict{String,Any}()
for (name, rtol, order) in CONFIGS
    println("="^76)
    println("設定: ", name)
    println("="^76)
    cb, rec = make_recorder(WINDOWS)
    t_wall = @elapsed (sol, p, u0) = solve_transient_lhp!(geo, condition, tables;
        reltol = rtol, max_order = order, callback = cb, verbose = false)
    println("  retcode=", sol.retcode, ", 計算時間 ", round(t_wall, digits = 1), " s, 記録したステップ数 ", length(rec[:t]))
    out = []
    for w in WINDOWS
        idx = findall(t -> w[1] <= t <= w[2], rec[:t])
        isempty(idx) && (println("  窓 ", w, ": 記録なし"); continue)
        tt = rec[:t][idx]
        dts = diff(tt)
        println("  --- 窓 t=", w[1], "〜", w[2], " s (ステップ ", length(idx), ", 平均刻み ", round(mean(dts) * 1e6, digits = 0), " μs) ---")
        for (lab, key, unit) in (("液管の平均速度", :u_ll_mean, "m/s"), ("液管の最上流の速度", :u_ll_first, "m/s"),
                                  ("液管の最下流の速度", :u_ll_last, "m/s"), ("凝縮器の平均速度", :u_cond_mean, "m/s"),
                                  ("液管の P−P_cc (平均)", :dP_ll_mean, "Pa"), ("凝縮器の P−P_cc (平均)", :dP_cond_mean, "Pa"))
            y = rec[key][idx]
            half = length(y) ÷ 2
            sd1 = std(y[1:half]); sd2 = std(y[half+1:end])
            println("    ", rpad(lab, 24), " 平均 ", lpad(round(mean(y), sigdigits = 4), 10), " ", unit,
                    ", 標準偏差 ", lpad(round(std(y), sigdigits = 3), 9), " (前半 ", round(sd1, sigdigits = 3), " / 後半 ", round(sd2, sigdigits = 3),
                    "), 最大−最小 ", round(maximum(y) - minimum(y), sigdigits = 3))
        end
        pk = spectrum_peaks(tt, rec[:u_ll_mean][idx], w)
        println("    液管の平均速度のスペクトルの上位 (周波数 [Hz], 振幅 [m/s]): ",
                join(["($(round(f, digits = 2)), $(round(a, sigdigits = 3)))" for (f, a) in pk], ", "))
        pkp = spectrum_peaks(tt, rec[:dP_ll_mean][idx], w)
        println("    液管の P−P_cc のスペクトルの上位 (周波数 [Hz], 振幅 [Pa]):    ",
                join(["($(round(f, digits = 2)), $(round(a, sigdigits = 3)))" for (f, a) in pkp], ", "))
        println("    凝縮器の蒸気セル数 (x>0.5): ", Int(minimum(rec[:n_vap][idx])), "〜", Int(maximum(rec[:n_vap][idx])),
                ", 液管の最上流の温度 ", round(mean(rec[:T_ll_first][idx]), digits = 3), " K")
        push!(out, (w = w, sd_u = std(rec[:u_ll_mean][idx]), sd_P = std(rec[:dP_ll_mean][idx]),
                    pk = isempty(pk) ? (NaN, NaN) : pk[1], steps = length(idx)))
    end
    results[name] = out
end

println()
println("="^76)
println("まとめ: 液管の平均速度の標準偏差 [m/s] と、液管の P−P_cc の標準偏差 [Pa]、スペクトルの最大のピーク")
println("="^76)
for w in WINDOWS
    println("  窓 t=", w[1], "〜", w[2], " s:")
    println("    ", rpad("設定", 28), lpad("σ_u [m/s]", 12), lpad("σ_P [Pa]", 11), lpad("ピーク[Hz]", 11), lpad("ステップ", 10))
    for (name, _, _) in CONFIGS
        haskey(results, name) || continue
        r = findfirst(o -> o.w == w, results[name])
        r === nothing && continue
        o = results[name][r]
        println("    ", rpad(name, 28), lpad(round(o.sd_u, sigdigits = 3), 12), lpad(round(o.sd_P, sigdigits = 3), 11),
                lpad(round(o.pk[1], digits = 2), 11), lpad(o.steps, 10))
    end
end
println()
println("読み方:")
println("  ・rtol を 1e-4 → 1e-5 → 1e-6 と小さくして、σ_u, σ_P がほぼ比例して小さくなれば H1 (許容誤差レベルの雑音)")
println("  ・σ が rtol に依存せず、狭いピークが同じ周波数に出れば H2 (物理的な振動) または モデルの不連続")
println("  ・次数1 (後退オイラー、高周波の減衰が強い) で σ が大きく減るなら、高周波の成分が主 (数値的な影響が大きい)")
println("  ・窓の前半と後半の標準偏差が同じなら、定常的な振動 (減衰しない)。後半が小さければ減衰中")