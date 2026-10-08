# =============================================================================
# test_diagnose_fullscale.jl
#
# 実運転の規模 (transient.jl と同じ幾何: 蒸気管20, 凝縮器50, 液管20セル, 650変数) で、
#   1. 計算し、IDA の各ステップの時刻を、コールバックで記録する
#   2. 刻みの時間分布 (どの時間帯に、ステップが集中するか)
#   3. 物理の確認 (時刻ごとの表)
#   4. エネルギー収支
#   5. 残差評価 1 回のコスト (時間、メモリ確保、プロファイル)
# を調べる。
#
# 背景: transient.jl は 200 s の計算を 1674 s で完走した (受理 90.5 万ステップ、
# 平均刻み 221 μs、残差の評価 190 万回、1 回あたり 0.88 ms)。速くするには、
#   (a) 1 回の残差評価を速くする  → 5 のプロファイルが示す
#   (b) 評価回数を減らす          → 2 の刻みの分布が示す
#
# 実行時間の目安: T_END=100 s で、10〜15 分。解は先にファイルに保存するので、
# 分析で失敗しても、計算をやり直さなくてよい (SOL_FILE)。
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays,
      SparseArrays, LinearAlgebra, Random, Sundials, Statistics, Profile
import DiffEqBase      # (`using DiffEqBase: 名前...` は、列挙した名前だけを取り込み、モジュール名 DiffEqBase は束縛しない。
                       #  DiffEqBase.u_modified! のように、モジュール名を直接使うので、import も必要)
using DiffEqBase: DAEProblem, DAEFunction, solve, DiscreteCallback
# u_modified! は、SciMLBase の関数 (DiffEqBase が再エクスポートする)。コールバックで使うので、先に確認する
isdefined(DiffEqBase, :u_modified!) ||
    error("DiffEqBase.u_modified! が見つかりません。SciMLBase を、プロジェクトに追加して、" *
          "`using SciMLBase` し、コールバックの DiffEqBase.u_modified! を SciMLBase.u_modified! に置き換えてください")

# ---------------------------- 設定 ----------------------------
const T_END = 100.0        # Q のステップアップ (t=50 s) の後も含める
const WINDOW = 5.0         # 刻みの分布を集計する時間幅 [s]
const SOL_FILE = joinpath(@__DIR__, "diagnose_fullscale_sol.jls")
# --------------------------------------------------------------

const SRC = joinpath(@__DIR__, "..", "src")
for f in ("physics", "smoothing", "types", "proptables", "propaccess", "correlations",
          "heattransfer", "geometry_cc", "residuals_transport", "residuals_wall",
          "residuals_wick", "residuals_cc", "residuals_evapcase", "residuals_grv",
          "residuals_main", "initialcondition", "sparsity", "solverfuncs", "postprocess")
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
    t_init = 0.0, t_end = T_END, saveat = 1.0)

# =============================================================================
println("="^76)
println("1. 計算 (IDA の各ステップの時刻を記録)")
println("="^76)
step_times = Float64[]
# 【重要】 affect! の後、積分器は、既定で「u が書き換えられたかもしれない」と見なして、IDA を
# 再初期化し (reinit)、そのとき、現在の u, du で CheckInit を行う。IDA の現在の du は BDF の
# 推定値で、厳密には整合していないため、"Initial conditions do not satisfy the DAE constraints"
# で失敗する (前回の実行)。u を変えないコールバックでは、u_modified!(integrator, false) を呼ぶ。
# (呼ばないと、記録のたびに IDA が再初期化され、次数 1・小刻みから立ち上げ直すので、
#  計算が遅くなり、結果も変わる)
cb = DiscreteCallback((u, t, integrator) -> true,
                      integrator -> (push!(step_times, integrator.t);
                                     DiffEqBase.u_modified!(integrator, false); nothing),
                      save_positions = (false, false))
sol = nothing; p = nothing; u0 = nothing; t_wall = 0.0
try
    global t_wall = @elapsed begin
        global sol, p, u0 = solve_transient_lhp!(geo, condition, tables; callback = cb)
    end
catch err
    println("  コールバックつきの計算に失敗: ", sprint(showerror, err)[1:min(300, end)])
    println("  → コールバックなしで、やり直す (刻みの分布は得られない)")
    empty!(step_times)
    global t_wall = @elapsed begin
        global sol, p, u0 = solve_transient_lhp!(geo, condition, tables)
    end
end
println("  計算時間 ", round(t_wall, digits=1), " s (シミュレーション ", T_END, " s), retcode=", sol.retcode)
try
    serialize(SOL_FILE, (t = sol.t, u = sol.u))
    println("  解を保存: ", SOL_FILE)
catch err
    println("  解の保存に失敗: ", sprint(showerror, err)[1:min(200, end)])
end

# =============================================================================
println()
println("="^76)
println("2. 刻みの時間分布")
println("="^76)
try
    if isempty(step_times)
        println("  (ステップの時刻が記録されていません)")
    else
        ts = sort(step_times)
        dts = diff(ts)
        println("  記録したステップ数 ", length(ts), ", 平均刻み ", round(mean(dts) * 1e6, digits=1), " μs")
        println("  刻みのパーセンタイル [μs]: 5% ", round(quantile(dts, 0.05) * 1e6, digits=1),
                ", 50% ", round(quantile(dts, 0.5) * 1e6, digits=1),
                ", 95% ", round(quantile(dts, 0.95) * 1e6, digits=1),
                ", 最小 ", round(minimum(dts) * 1e6, digits=2))
        println("  時間幅 ", WINDOW, " s ごとの、ステップ数と平均刻み (ステップが集中する時間帯を探す):")
        println("    時間帯 [s]      ステップ数   平均刻み[μs]   全体に占める割合")
        edges = collect(0.0:WINDOW:T_END)
        total = length(ts)
        for k in 1:length(edges)-1
            n = count(t -> edges[k] <= t < edges[k+1], ts)
            avg = n > 0 ? WINDOW / n * 1e6 : NaN
            bar = repeat("#", round(Int, 40 * n / total))
            println("    ", lpad(edges[k], 5), "〜", lpad(edges[k+1], 5), "   ", lpad(n, 9), "   ", lpad(round(avg, digits=1), 10),
                    "      ", lpad(round(100 * n / total, digits=1), 5), " % ", bar)
        end
    end
catch err
    println("  失敗: ", sprint(showerror, err)[1:min(300, end)])
end

# =============================================================================
println()
println("="^76)
println("3. 物理の確認 (時刻ごと)")
println("="^76)
try
    println("  t[s] | Q[W]  | T_ec[K] | T_out[K] | T_cc[K] | P_cc[Pa]  | P_grv−P_cc[Pa] | ṁ[kg/s]  | ṁ h_fg[W] | 凝縮器の蒸気セル数 | CC液体積分率 | 凝縮器出口T[K]")
    for tsel in (0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 45.0, 49.0, 50.0, 51.0, 52.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0)
        tsel > sol.t[end] + 1e-9 && continue
        k = argmin(abs.(sol.t .- tsel))
        c = components_at(sol, p, k)
        P_grv = groove_pressure(c.grv)
        T_out = Tsat(P_grv, tables)
        ρ_l = tables.ρ_l_sat(c.cc.T[1])
        ṁ = c.cc.u[1] * ρ_l * p.p_cc.A
        h_fg = tables.h_v_sat(T_out) - tables.h_l_sat(T_out)
        n_vap = count(>(0.5), c.cond.x)
        sp = sat_props(c.cc.T[1], tables)
        f_l = liquid_volume_fraction(c.cc.ρm[1], c.cc.x[1], sp.ρ_l)
        println("  ", lpad(round(sol.t[k], digits=1), 5), " | ", lpad(round(applied_heat(sol.t[k], condition), digits=1), 5),
                " | ", round(c.T_wall.evap[1], digits=3), " | ", round(T_out, digits=3), " | ", round(c.cc.T[1], digits=3),
                " | ", round(c.cc.P[1], digits=0), " | ", lpad(round(P_grv - c.cc.P[1], digits=0), 14),
                " | ", round(ṁ, sigdigits=3), " | ", lpad(round(ṁ * h_fg, digits=2), 9),
                " | ", lpad(n_vap, 18), " | ", lpad(round(f_l, digits=4), 12), " | ", round(c.cond.T[end], digits=3))
    end
    println("  (ṁ h_fg ≈ Q_apply なら、熱が全て蒸発に使われている。凝縮器の蒸気セル数は、x>0.5 の凝縮器のセル数 (全50)。")
    println("   ステップアップ (t=50〜51 s) の後、T_ec と ṁ が、数秒〜数十秒かけて新しい値に落ち着くはず)")
catch err
    println("  失敗: ", sprint(showerror, err)[1:min(300, end)])
end

# =============================================================================
println()
println("="^76)
println("4. エネルギー収支")
println("="^76)
try
    energy_balance(sol, p)
    println("  (残差が入力の数 % 以内なら、エネルギーは閉じている)")
catch err
    println("  失敗: ", sprint(showerror, err)[1:min(300, end)])
end

# =============================================================================
println()
println("="^76)
println("5. 残差評価 1 回のコスト")
println("="^76)
try
    n = length(sol.u[1])
    u = copy(sol.u[end]); du = zeros(n); res = zeros(n)
    t_eval = sol.t[end]
    for _ in 1:20; f_dae_residual!(res, du, u, p, t_eval); end      # JIT のウォームアップ
    st = @timed for _ in 1:500; f_dae_residual!(res, du, u, p, t_eval); end
    alloc = @allocated f_dae_residual!(res, du, u, p, t_eval)
    println("  変数 ", n, " 個: 1 回あたり ", round(st.time / 500 * 1e3, digits=3), " ms (GC の割合 ",
            round(100 * st.gctime / st.time, digits=1), " %), メモリ確保 ", round(alloc / 1024, digits=1), " KiB/回")
    println("  → 1 変数あたり ", round(st.time / 500 / n * 1e6, digits=2), " μs")
    println()
    println("  プロファイル (3000 回の評価。右端 = サンプル数の多い順。Count は包含的 (子の関数を含む)):")
    Profile.clear()
    Profile.init(n = 10^7, delay = 0.0005)
    @profile for _ in 1:3000; f_dae_residual!(res, du, u, p, t_eval); end
    io = IOBuffer()
    Profile.print(IOContext(io, :displaysize => (100, 220)); format = :flat, sortedby = :count, mincount = 30, noisefloor = 2.0)
    lines = split(String(take!(io)), '\n')
    # 自分のコード (src の .jl) の行だけを抜き出して、サンプル数の多い順に表示
    mine = filter(l -> occursin(".jl", l), lines)       # (ファイル名は短縮されることがあるので、.jl だけで絞る)
    println("   Count  ファイル:行  関数")
    for l in mine[max(1, end - 39):end]
        println("  ", l)
    end
catch err
    println("  失敗: ", sprint(showerror, err)[1:min(300, end)])
end
println()
println("完了。この出力をそのまま共有してください。")