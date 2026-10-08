# =============================================================================
# test_grid_convergence.jl (空間格子の収束の確認)
#
# 目的: 起動時の動作を正確に再現するために、凝縮器・蒸気管・液管の格子の細かさが、結果を
# 左右するかを調べる。凝縮器 50 セル (2 cm) は、起動時に凝縮の前線が通過する場所で、
# 粗い可能性がある。
#
# 方法: 格子を、基準 → 2 倍 → 4 倍と細かくして、同じ条件 (起動 → t=50 s で 20 W → 40 W) で
# 計算し、次を比べる。
#   ・起動の指標: 蒸発量が最終値の 90% に達する時間、起動時の圧力のピーク
#   ・蒸発器ケースの温度、蒸発量、凝縮器の蒸気の長さ ∫x dz (前線の位置の、格子に依らない指標)
#   ・液管の速度の振動 (σ_u、周波数): 格子で変わるなら、数値的な成分を含む
#   ・3 つの格子があれば、収束の次数 p = log2 |d1/d2|  (d1 = 2倍−基準, d2 = 4倍−2倍)
#
# 結果はケースごとに OUT_DIR に保存する。RUN_CASES で、実行するケースを選び (既存のものは、
# FORCE_RERUN=false なら飛ばす)、解析は、保存済みの全ケースを使う。
# 解析だけ行うなら、RUN_CASES = String[] にする。
#
# 実行時間 (実測。T_RUN=70 s): 基準 390 s (650 変数)、2 倍 1087 s (1280 変数)。4 倍は、
# 変数が約 2600 個で、2 倍の 3 倍以上 (1 時間以上) かかる見込み。
# 一部だけ細かくしたケース (cond2x: 約 1000 変数、ll2x: 約 790 変数) は、数分〜10 分程度の見込み。
#
# 【基準・2 倍の結果 (実測)】 T_ec, ṁ, P_cc は、0.5% 以内で収束。L_vap (前線の位置) は、
# 基準が 2 倍より、10 s で 3.8%、45 s で 1.7% 短い。ステップ後の σ(u_ll) (60〜64 s) は、
# 1.1% の差で収束。一方、ステップ後の max|ΔP_ll| は 655 → 3605 Pa (5.5 倍) で、収束していない。
# (σ(u_ll) の 40〜44 s の窓は、rtol を変えても 2.8 倍ばらつく窓で、格子の影響とは断定できない。)
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays,
      SparseArrays, LinearAlgebra, Random, Sundials, Statistics
using DiffEqBase: DAEProblem, DAEFunction, solve

# ---------------------------- 設定 ----------------------------
const T_RUN = 70.0                # 計算する時間 [s]
const T_STEP_Q = 50.0             # Q のステップアップの時刻 [s]
const SAVE_DT = 0.05              # 保存間隔 [s] (saveat は、刻みに影響しない)
# (名前, 蒸気管のセル数, 凝縮器のセル数, 液管のセル数)。長さは、蒸気管 2 m, 凝縮器 1 m, 液管 2 m
const CASES = [
    ("base",  20, 50, 20),
    ("2x",    40, 100, 40),
    ("4x",    80, 200, 80),
    ("cond2x", 20, 100, 20),      # 凝縮器だけ 2 倍 (パルスの増加が、凝縮器の格子によるか)
    ("ll2x",   20, 50, 40),      # 液管だけ 2 倍 (パルスの増加が、液管の格子によるか)
]
const RUN_CASES = ["cond2x", "ll2x"]   # 今回、実行するケース (base, 2x は保存済み。4x は重いので、必要になったら追加)
const FORCE_RERUN = false         # true なら、保存済みのケースも、やり直す
const OUT_DIR = joinpath(@__DIR__, "grid_convergence")
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
condition = LHPTransientCondition(
    fluid = fluid, T_cond = 15.0 + 273.15, T_amb = 25.0 + 273.15, T_init = 25.0 + 273.15,
    Q_apply_init = 20.0, Q_apply = 40.0, t_step = T_STEP_Q, τ_step = 1.0, ΔT_startup = 0.0,
    t_init = 0.0, t_end = T_RUN, saveat = SAVE_DT)

make_geo(N_vl, N_cond, N_ll) = LHPGeometry(
    WickGeometry(N=3, φ=0.65, ρ_wick=8900.0, Cp_wick=444.0),
    GrvGeometry(N=1),
    TLGeometry(N=N_vl, L=2, D_i=0.00435, D_o=0.00635),
    CondGeometry(N=N_cond, L=1.0, D_i=0.00435, D_o=0.00635),
    TLGeometry(N=N_ll, L=2, D_i=0.00218, D_o=0.00318),
    CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037))

# --- 保存した解から、比べる量を取り出す ---
function extract(sol, p, tables)
    n = length(sol.t)
    T_ec = zeros(n); P_cc = zeros(n); P_grv = zeros(n); T_cc = zeros(n); mdot = zeros(n)
    L_vap = zeros(n); u_ll = zeros(n); dP_ll = zeros(n); u_cond = zeros(n); T_cout = zeros(n)
    n_liq = zeros(Int, n); n_dry = zeros(Int, n)      # 凝縮器の、液 (x<1e-3) と、乾き蒸気 (x>0.999) のセル数
    for k in 1:n
        c = components_at(sol, p, k)
        T_ec[k] = c.T_wall.evap[1]; P_cc[k] = c.cc.P[1]; T_cc[k] = c.cc.T[1]
        P_grv[k] = groove_pressure(c.grv)
        mdot[k] = c.cc.u[1] * tables.ρ_l_sat(c.cc.T[1]) * p.p_cc.A
        L_vap[k] = sum(c.cond.x) * p.p_cond.dl          # 凝縮器の蒸気の長さ ∫x dz [m]
        u_ll[k] = mean(c.ll.u); dP_ll[k] = mean(c.ll.P) - c.cc.P[1]
        u_cond[k] = mean(c.cond.u); T_cout[k] = c.cond.T[end]
        n_liq[k] = count(<(1e-3), c.cond.x); n_dry[k] = count(>(0.999), c.cond.x)
    end
    return (t = copy(sol.t), T_ec = T_ec, P_cc = P_cc, P_grv = P_grv, T_cc = T_cc, mdot = mdot,
            L_vap = L_vap, u_ll = u_ll, dP_ll = dP_ll, u_cond = u_cond, T_cout = T_cout,
            n_liq = n_liq, n_dry = n_dry)
end

# --- ケースの実行と保存 ---
mkpath(OUT_DIR)
case_file(label) = joinpath(OUT_DIR, "case_$(label).jls")
for (label, N_vl, N_cond, N_ll) in CASES
    label in RUN_CASES || continue
    if isfile(case_file(label)) && !FORCE_RERUN
        println("ケース ", label, ": 保存済み (", case_file(label), ") なので、飛ばします")
        continue
    end
    println("="^76); println("ケース ", label, ": 蒸気管 ", N_vl, ", 凝縮器 ", N_cond, ", 液管 ", N_ll, " セル"); println("="^76)
    geo = make_geo(N_vl, N_cond, N_ll)
    t_wall = @elapsed (sol, p, u0) = solve_transient_lhp!(geo, condition, tables; verbose = false)
    ser = extract(sol, p, tables)
    st = sol.stats
    res = merge(ser, (label = label, N_vl = N_vl, N_cond = N_cond, N_ll = N_ll, n_var = length(sol.u[1]),
                      t_wall = t_wall, retcode = string(sol.retcode), naccept = st.naccept, nreject = st.nreject,
                      njacs = st.njacs, nconvfail = st.nnonlinconvfail))
    serialize(case_file(label), res)
    println("  retcode=", res.retcode, ", 変数 ", res.n_var, " 個, 計算時間 ", round(t_wall, digits = 1), " s, 受理ステップ ", res.naccept,
            ", 棄却 ", res.nreject, ", 収束失敗 ", res.nconvfail, ", ヤコビアン ", res.njacs)
    println("  保存: ", case_file(label))
end

# =============================================================================
# 解析 (保存済みの全ケース)
# =============================================================================
results = []
for (label, _, _, _) in CASES
    isfile(case_file(label)) && push!(results, deserialize(case_file(label)))
end
labels = [r.label for r in results]
println()
println("="^76); println("解析: 保存済みのケース ", labels); println("="^76)
isempty(results) && error("保存済みのケースがありません。RUN_CASES を設定して、実行してください")

# 差を取る相手 (基準の格子): 最も細かい一様な格子。4x があれば 4x、なければ 2x、なければ base、
# それも無ければ最初のケース (一部だけ細かくしたケースを、基準にしないため)
ref_idx = something(findfirst(==("4x"), labels), findfirst(==("2x"), labels), findfirst(==("base"), labels), 1)
ref = results[ref_idx]
println("差を取る相手 (基準の格子): ", ref.label)

ix(r, tq) = argmin(abs.(r.t .- tq))
win_std(r, key, w) = std(getfield(r, key)[findall(τ -> w[1] <= τ <= w[2], r.t)])
win_max(r, key, w) = maximum(abs.(getfield(r, key)[findall(τ -> w[1] <= τ <= w[2], r.t)]))

function t_reach(r, key, frac, tref)       # 値が、tref での値の frac 倍に初めて達する時刻
    y = getfield(r, key); target = frac * y[ix(r, tref)]
    k = findfirst(>=(target), y)
    return k === nothing ? NaN : r.t[k]
end

# 窓の中の時系列の、スペクトルの上位のピーク (周波数 [Hz], 振幅)。0.25 Hz 刻み
function top_peaks(t, y, w; npeaks = 3)
    idx = findall(τ -> w[1] <= τ <= w[2], t)
    ts = t[idx]; yy = y[idx] .- mean(y[idx])
    freqs = collect(0.25:0.25:min(0.45 / SAVE_DT, 10.0))
    amp = [abs(sum(yy[j] * cis(-2π * f * (ts[j] - ts[1])) for j in eachindex(ts))) * 2 / length(ts) for f in freqs]
    pk = Tuple{Float64,Float64}[]
    for k in eachindex(freqs)
        lo = k > 1 ? amp[k-1] : -Inf; hi = k < length(freqs) ? amp[k+1] : -Inf
        (amp[k] > lo && amp[k] >= hi) && push!(pk, (freqs[k], amp[k]))
    end
    sort!(pk; by = x -> -x[2])
    return first(pk, min(npeaks, length(pk)))
end

# 液管の圧力 ΔP_ll の、中央値からのずれが大きいサンプル (互いに sep [s] 以上離れた上位 n 個)
function top_spikes(r, w; n = 5, sep = 0.5)
    idx = findall(τ -> w[1] <= τ <= w[2], r.t)
    med = median(r.dP_ll[idx])
    dev = abs.(r.dP_ll[idx] .- med)
    ord = sortperm(dev; rev = true)
    picked = Int[]
    for o in ord
        k = idx[o]
        all(abs(r.t[k] - r.t[q]) > sep for q in picked) || continue
        push!(picked, k)
        length(picked) >= n && break
    end
    return picked
end

const METRICS = [
    ("T_ec(5 s) [K]",                r -> r.T_ec[ix(r, 5.0)]),
    ("T_ec(20 s) [K]",               r -> r.T_ec[ix(r, 20.0)]),
    ("T_ec(45 s) [K]",               r -> r.T_ec[ix(r, 45.0)]),
    ("T_ec(60 s) [K]",               r -> r.T_ec[ix(r, 60.0)]),
    ("ṁ(2 s) [mg/s]",                r -> 1e6 * r.mdot[ix(r, 2.0)]),
    ("ṁ(5 s) [mg/s]",                r -> 1e6 * r.mdot[ix(r, 5.0)]),
    ("ṁ(45 s) [mg/s]",               r -> 1e6 * r.mdot[ix(r, 45.0)]),
    ("ṁ(60 s) [mg/s]",               r -> 1e6 * r.mdot[ix(r, 60.0)]),
    ("ṁ が 45 s の値の 90% に達する時刻 [s]", r -> t_reach(r, :mdot, 0.9, 45.0)),
    ("起動時の P_grv−P_cc の最大 (t<5 s) [Pa]", r -> maximum(r.P_grv[1:ix(r, 5.0)] .- r.P_cc[1:ix(r, 5.0)])),
    ("L_vap(10 s) [m]",              r -> r.L_vap[ix(r, 10.0)]),
    ("L_vap(45 s) [m]",              r -> r.L_vap[ix(r, 45.0)]),
    ("L_vap(69 s) [m]",              r -> r.L_vap[ix(r, 69.0)]),
    ("P_cc(45 s) [Pa]",              r -> r.P_cc[ix(r, 45.0)]),
    ("P_cc(69 s) [Pa]",              r -> r.P_cc[ix(r, 69.0)]),
    ("σ(u_ll) t=40〜44 s [mm/s] (窓が短く、ばらつく)", r -> 1e3 * win_std(r, :u_ll, (40.0, 44.0))),
    ("σ(u_ll) t=60〜64 s [mm/s]",    r -> 1e3 * win_std(r, :u_ll, (60.0, 64.0))),
    ("max|ΔP_ll| t=10〜50 s (ステップ前) [Pa]", r -> win_max(r, :dP_ll, (10.0, 50.0))),
    ("max|ΔP_ll| t=50〜70 s (ステップ後) [Pa]", r -> win_max(r, :dP_ll, (50.0, 70.0))),
]

println()
println("計算コスト:")
println("  ケース     変数    計算時間[s]   受理ステップ   棄却   収束失敗   ヤコビアン")
for r in results
    println("  ", rpad(r.label, 9), lpad(r.n_var, 6), lpad(round(r.t_wall, digits = 0), 12), lpad(r.naccept, 14),
            lpad(r.nreject, 8), lpad(r.nconvfail, 10), lpad(r.njacs, 12))
end

println()
println("比べる量 (ケースごと):")
println("  ", rpad("量", 52), join([lpad(l, 11) for l in labels], ""))
vals = [[f(r) for r in results] for (_, f) in METRICS]
for ((name, _), v) in zip(METRICS, vals)
    println("  ", rpad(name, 52), join([lpad(round(x, sigdigits = 5), 11) for x in v], ""))
end

# 基準の格子 (ref) との差: 各ケースについて、差 Δ と、相対差 [%]
println()
println("基準の格子 (", ref.label, ") との差: 各ケースの (値 − 基準の値) と、相対差 [%]")
others = [i for i in eachindex(results) if i != ref_idx]
println("  ", rpad("量", 52), join([lpad(results[i].label, 11) * lpad("[%]", 7) for i in others], ""))
for ((name, _), v) in zip(METRICS, vals)
    cells = String[]
    for i in others
        d = v[i] - v[ref_idx]
        rel = 100 * abs(d) / max(abs(v[ref_idx]), eps())
        push!(cells, lpad(round(d, sigdigits = 3), 11) * lpad(round(rel, digits = 1), 7))
    end
    println("  ", rpad(name, 52), join(cells, ""))
end

# 収束の次数 (base, 2x, 4x が、そろっているとき)
if all(l -> l in labels, ("base", "2x", "4x"))
    ib, i2, i4 = findfirst(==("base"), labels), findfirst(==("2x"), labels), findfirst(==("4x"), labels)
    println()
    println("収束の次数 p = log2|d1/d2| (d1 = 2x − base, d2 = 4x − 2x。p≈1 が 1 次風上の期待値。符号が逆転なら '-'):")
    for ((name, _), v) in zip(METRICS, vals)
        d1 = v[i2] - v[ib]; d2 = v[i4] - v[i2]
        p_ord = (d1 * d2 > 0 && abs(d2) > 0) ? log2(abs(d1 / d2)) : NaN
        println("  ", rpad(name, 52), lpad(isnan(p_ord) ? "-" : string(round(p_ord, digits = 2)), 8))
    end
end

# 時系列全体の差 (共通の時刻で、基準の格子との相対 RMS 差。変動幅に対する割合)
println()
println("時系列 (t=0〜", T_RUN, " s) の、基準の格子 (", ref.label, ") との、相対 RMS 差 (変動幅に対する割合 [%]):")
tc = collect(0.0:0.5:T_RUN)
println("  ", rpad("ケース", 9), join([lpad(s, 10) for s in ("T_ec", "ṁ", "L_vap", "P_cc")], ""))
for i in others
    r = results[i]; row = String[]
    for key in (:T_ec, :mdot, :L_vap, :P_cc)
        a = [getfield(r, key)[ix(r, τ)] for τ in tc]; b = [getfield(ref, key)[ix(ref, τ)] for τ in tc]
        rng = maximum(b) - minimum(b)
        push!(row, lpad(rng > 0 ? round(100 * sqrt(mean((a .- b) .^ 2)) / rng, digits = 2) : NaN, 10))
    end
    println("  ", rpad(r.label, 9), join(row, ""))
end

# 液管の速度のスペクトル: 上位 3 つのピーク
println()
println("液管の平均速度のスペクトル: 上位 3 つのピーク (周波数 [Hz], 振幅 [mm/s])。0.25 Hz 刻みで、上位が入れ替わりやすいので、並べて見る:")
for w in ((40.0, 44.0), (60.0, 64.0))
    println("  窓 t=", w[1], "〜", w[2], " s:")
    for r in results
        pk = top_peaks(r.t, r.u_ll, w)
        println("    ", rpad(r.label, 9), join(["($(f), $(round(1e3 * a, sigdigits = 3)))" for (f, a) in pk], "  "))
    end
end

# 液管の圧力パルス: 上位 5 個の時刻と大きさ、最大のパルス (ステップ後) の形
println()
println("="^76)
println("液管の圧力 ΔP_ll (液管のセル平均 − P_cc) の、突発的なパルス")
println("="^76)
println("  (ジューコフスキー則: ΔP = ρ c Δu。液の ρ c ≈ 8.4e5 kg/m²s → 1 mm/s の速度の跳びで 840 Pa)")
for r in results
    println()
    println("ケース ", r.label, ": 上位 5 個 (t=5〜", T_RUN, " s、互いに 0.5 s 以上離れたもの):")
    sp = top_spikes(r, (5.0, T_RUN))
    for k in sp
        println("    t = ", lpad(round(r.t[k], digits = 2), 6), " s   ΔP_ll = ", lpad(round(r.dP_ll[k], digits = 0), 7),
                " Pa   (u_ll = ", round(1e3 * r.u_ll[k], digits = 1), " mm/s)")
    end
    sp_post = top_spikes(r, (T_STEP_Q, T_RUN); n = 1)
    if !isempty(sp_post)
        k = sp_post[1]
        println("  ステップ後の最大のパルスの形 (t=", round(r.t[k], digits = 2), " s の ±0.5 s、0.1 s 刻みの ΔP_ll [Pa]):")
        println("    ", join([string(round(r.dP_ll[ix(r, τ)], digits = 0)) for τ in (r.t[k]-0.5):0.1:(r.t[k]+0.5)], ", "))
    end
    if hasproperty(r, :n_liq)
        chg = [r.t[i] for i in 2:length(r.t) if (r.n_liq[i] != r.n_liq[i-1] || r.n_dry[i] != r.n_dry[i-1]) && r.t[i] > 5.0]
        println("  凝縮器の前線のセルの状態が変わった時刻 (液セル数または乾きセル数が変化。t>5 s): ", length(chg), " 回")
        println("    ", join([string(round(c, digits = 2)) for c in chg[1:min(end, 30)]], ", "), length(chg) > 30 ? ", ..." : "")
        if !isempty(chg) && !isempty(sp)
            near = [minimum(abs.(chg .- r.t[k])) for k in sp]
            println("  上位 5 個のパルスから、最も近い状態変化までの時間 [s]: ", join([string(round(x, digits = 2)) for x in near], ", "))
        end
    else
        println("  (このケースの保存結果には、凝縮器の液セル数・乾きセル数が入っていません。新しく実行したケースだけ、照合できます)")
    end
end

println()
println("読み方:")
println("  ・T_ec, ṁ, L_vap, P_cc (巨視的な量) の相対差が、数 % 以下なら、その格子で収束している。")
println("    L_vap (凝縮器の蒸気の長さ) は、前線の位置の指標で、起動時の前線の進行の再現に関係する")
println("  ・ステップ後の σ(u_ll) (60〜64 s) が、格子で変わらなければ、ステップで励起される大振幅の振動は、")
println("    物理 (液柱の質量−ばね振動) の性質。σ(u_ll) の 40〜44 s の窓は、窓が短く、ばらつきが大きい (rtol を変えても 2.8 倍)")
println("  ・圧力パルス (max|ΔP_ll|) が、cond2x と ll2x の、どちらで大きくなるか (基準との差が大きい方) で、")
println("    パルスを決める格子 (凝縮器か液管か) が分かる。パルスの形が、狭く鋭いなら、セルの大きさで決まる事象")
println("  ・ケースが 3 つ (base, 2x, 4x) あれば、収束の次数が求まる")