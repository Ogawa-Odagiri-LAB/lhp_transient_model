# =============================================================================
# test_diagnose_domain.jl
#
# 低熱負荷 (例: Q_apply_init=5 W → Q_apply=10 W) で、凝縮器の h_cond (Shah) の
#   sqrt(g D ρ_v (ρ_l − ρ_v))  が負になって落ちる原因を切り分ける診断。
#
# 【前提の整理】
#   ・h_transport_wall は、セルの相にかかわらず、毎回、単相・Shah・Chen・Dougall-Rohsenow を
#     すべて計算して、tanh で滑らかにブレンドする (分岐を持たないため)。したがって、
#     単相のセルでも Shah の式は評価され、どのセルで落ちたかはスタックトレースに出ない。
#   ・sqrt の中身は、飽和密度 ρ_v_sat(T), ρ_l_sat(T) だけで決まる。乾き度 x が負でも、
#     この sqrt は落ちない (x は quality_for_correlation で [1e-4, 1−1e-4] に収められる)。
#     落ちるのは、(ρ_m, e) の flash が返す T が、物性テーブルの範囲 (230〜360 K) の外に出て、
#     飽和密度が線形外挿されて ρ_v<0 などになったときだけ。
#   ・(ρ_m, e) が表の外に出る経路の候補:
#       H1 圧縮液のセル: ∂P/∂ρ ≈ 1.1e6 Pa/(kg/m³) と硬いので、Newton の試行点で密度が
#          0.4% ずれるだけで P が P_max (4.5 MPa) を超え、表の外挿で T が異常になる
#          (凝縮器の大部分が液のとき、低熱負荷ほどこの状態になる)。
#       H2 ṁ<0 (蒸発器での凝縮) による逆流で、液が溜まりすぎる (ご提案の仮説)。
#       H3 Newton の試行点だけの一時的な飛び (収束後の解は正常)。
#       H4 本当の発散 (解が、物理的におかしな状態へ向かって進んでいる)。
#
# このスクリプトがすること:
#   1. 残差関数をラップして、全呼び出しのうち LOG_DT [s] ごとに、状態の要約 (ṁ、凝縮器・液管・CC
#      の T と P の範囲、凝縮器の蒸気セル数) を記録する。直近 RING 回の呼び出しの状態も保持する。
#   2. 最初の DomainError で、(a) 範囲外の flash 結果を持つセル (成分, 番号, ρm, Em, P, T, x,
#      飽和密度)、(b) そのセルの直近の呼び出しの履歴 (飛びか、ゆっくりした漂流か)、
#      (c) 直前の時系列 (ṁ が負になっていたか、P が上がっていたか) を表示する。
#   3. GUARD=true にすると、DomainError を捕まえて残差を NaN にして返し、IDA に刻みを縮めさせる。
#      最後まで走れば H3 (試行点だけの飛び)、走れなければ H4。捕まえた回数も表示する。
#   4. STARTUP_LIST に ΔT_startup を並べると、それぞれで実行して、ṁ<0 になる時刻・ṁ の最小値・
#      落ちる時刻を一覧にする (ΔT_startup が逆流を防ぐかの確認)。
#
# 使い方: ./test/ に置いて、プロジェクトのルートで  include("test/test_diagnose_domain.jl")
# =============================================================================

# 他の test_*.jl と同じ方式: ./src の各ファイルを、このスクリプトに直接 include する
# (モジュールの export に頼らない)。./test/ に置いて、プロジェクトのルートから実行する。
using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays,
      SparseArrays, LinearAlgebra, Random, Sundials, Statistics, Printf
import DiffEqBase
using DiffEqBase: DAEProblem, DAEFunction, solve

const SRC = joinpath(@__DIR__, "..", "src")
for f in ("physics", "smoothing", "types", "proptables", "propaccess", "correlations",
          "heattransfer", "geometry_cc", "geometry_flatcc", "residuals_transport", "residuals_wall",
          "residuals_wick", "residuals_cc", "residuals_evapcase", "residuals_grv", "residuals_main",
          "residuals_flatwick", "residuals_flatcc", "residuals_flatevapcase", "residuals_flatmain",
          "initialcondition", "sparsity", "solverfuncs", "postprocess")
    include(joinpath(SRC, f * ".jl"))
end
const M = @__MODULE__      # include した関数・型は、このスクリプトの名前空間 (Main) にある

# ---------------------------- 設定 ----------------------------
const GEOMETRY = :flat                 # :flat (平板型) または :cyl (円筒型)
const Q_INIT, Q_APPLY = 5.0, 10.0      # 印加熱 [W]
const T_STEP, TAU_STEP = 50.0, 1.0
const T_END = 200.0
const STARTUP_LIST = [0.0]             # 試す ΔT_startup [K] のリスト。例: [0.0, 0.5, 1.0, 2.0]
const GUARD = false                    # true: DomainError を捕まえて NaN を返す (H3 と H4 の見分け)
const T_LO, T_HI, P_HI = 230.0, 360.0, 4.5e6   # 物性テーブルの範囲 (init_fluid_lookups と同じ値)
const LOG_DT = 0.5                     # 状態の要約を記録する間隔 [s]
const RING = 40                        # 直近の残差評価の状態を保持する数
const N_PRINT = 30                     # 失敗時に表示する直前の時系列の行数
# --------------------------------------------------------------

tables = M.init_fluid_lookups("ammonia", T_min = T_LO, T_max = T_HI, P_max = P_HI)

function make_geo(kind::Symbol)
    # 運転と同じ幾何に合わせる (transient.jl / transient_flat.jl と同じ値にしてください)
    vl   = M.TLGeometry(N = 20, L = 1.17, D_i = 0.00176, D_o = 0.00318)
    cond = M.CondGeometry(N = 50, L = 0.56, D_i = 0.00176, D_o = 0.00318)
    ll   = M.TLGeometry(N = 20, L = 1.15, D_i = 0.00099, D_o = 0.00159)
    if kind === :flat
        return M.FlatLHPGeometry(
            M.FlatWickGeometry(N = 3, φ = 0.30, ρ_wick = 8900.0, Cp_wick = 444.0),
            M.SecondaryWickGeometry(φ_s = 0.70, ρ_swick = 8900.0, Cp_swick = 444.0),
            M.FlatGrvGeometry(N = 1), vl, cond, ll, M.FlatCCGeometry(N = 1))
    else
        return M.LHPGeometry(
            M.WickGeometry(N = 3, φ = 0.30, ρ_wick = 8900.0, Cp_wick = 444.0),
            M.GrvGeometry(N = 1), vl, cond, ll, M.CCGeometry(N = 1, L = 0.07, D_i = 0.035, D_o = 0.037))
    end
end

# ---------------- 診断用の状態と、ラップしたパラメータ ----------------
mutable struct DiagState
    ncalls::Int
    next_log::Float64
    log::Vector{Any}
    ring_t::Vector{Float64}
    ring_u::Vector{Vector{Float64}}
    nerr::Int
    reported::Bool
    t_first_err::Float64
end

struct DiagParams{P} <: M.AbstractLHPParams
    inner::P
    st::DiagState
end
Base.getproperty(d::DiagParams, s::Symbol) =
    (s === :inner || s === :st) ? getfield(d, s) : getproperty(getfield(d, :inner), s)

"""全流体セルの、(ρm, Em) から表の flash で求めた P, T, x。h_transport_wall と同じ経路。"""
function fluid_cells(pin, x)
    g = pin.cond.g
    out = NamedTuple[]
    for (name, pc) in ((:grv, pin.p_grv), (:vl, pin.p_vl), (:cond, pin.p_cond),
                       (:ll, pin.p_ll), (:cc, pin.p_cc))
        cm = getproperty(x, name)
        for i in 1:pc.N
            z = M.cell_elevation(pc, i)
            em = M.internal_energy(cm.Em[i], cm.u[i], z, g)
            P, T, xq = try
                M.flash_ρe(cm.ρm[i], em, pin.tables)
            catch
                (NaN, NaN, NaN)
            end
            push!(out, (comp = name, i = i, ρ = cm.ρm[i], u = cm.u[i], Em = cm.Em[i], em = em,
                        P = P, T = T, x = xq))
        end
    end
    return out
end

isbad(c) = !isfinite(c.T) || !isfinite(c.P) || c.T < T_LO || c.T > T_HI || c.P > P_HI || c.P < 1.0e4

function snapshot(pin, u, t)
    x = ComponentArray(u, pin.axes)
    cells = fluid_cells(pin, x)
    sel(name) = filter(c -> c.comp === name, cells)
    rng(name, f) = (v = [getfield(c, f) for c in sel(name)]; (minimum(v), maximum(v)))
    ρ_l = pin.tables.ρ_l_sat(x.cc.T[1])
    return (t = t, Q = M.applied_heat(t, pin.cond),
            mdot = x.cc.u[1] * ρ_l * pin.p_cc.A,
            T_ec = x.T_wall.evap[1], T_cc = x.cc.T[1], P_cc = x.cc.P[1],
            cond_T = rng(:cond, :T), cond_P = rng(:cond, :P),
            ll_T = rng(:ll, :T), ll_P = rng(:ll, :P),
            nvap_cond = count(>(0.5), x.cond.x), nvap_vl = count(>(0.5), x.vl.x),
            n_bad = count(isbad, cells))
end

function print_log_rows(rows)
    @printf("%7s %5s %10s %8s %8s %9s | %7s %7s %9s %9s | %7s %9s | %4s %4s\n",
            "t[s]", "Q[W]", "mdot[mg/s]", "T_ec", "T_cc", "P_cc[kPa]",
            "cdT_min", "cdT_max", "cdP_min", "cdP_max", "llT_min", "llP_max", "vap_c", "vap_v")
    for r in rows
        @printf("%7.2f %5.1f %10.4f %8.3f %8.3f %9.1f | %7.2f %7.2f %9.1f %9.1f | %7.2f %9.1f | %4d %4d\n",
                r.t, r.Q, r.mdot * 1e6, r.T_ec, r.T_cc, r.P_cc / 1e3,
                r.cond_T[1], r.cond_T[2], r.cond_P[1] / 1e3, r.cond_P[2] / 1e3,
                r.ll_T[1], r.ll_P[2] / 1e3, r.nvap_cond, r.nvap_vl)
    end
end

function report_failure(pin, st::DiagState, u, t, err)
    println("\n", "="^78)
    println("DomainError を捕まえました: t = ", t, " s (残差評価 ", st.ncalls, " 回目)")
    println("  ", sprint(showerror, err))
    x = ComponentArray(u, pin.axes)
    cells = fluid_cells(pin, x)
    bads = filter(isbad, cells)
    println("\n(a) 範囲外の flash 結果を持つセル (T: ", T_LO, "〜", T_HI, " K, P ≤ ", P_HI / 1e6, " MPa):  ", length(bads), " セル")
    for c in first(bads, 8)
        ρl = pin.tables.ρ_l_sat(c.T); ρv = pin.tables.ρ_v_sat(c.T)
        @printf("   %-5s #%-3d ρm=%9.3f u=%9.2e Em=%11.1f em=%11.1f -> P=%11.1f Pa T=%8.2f K x=%4.2f | ρl_sat=%8.2f ρv_sat=%8.3f\n",
                string(c.comp), c.i, c.ρ, c.u, c.Em, c.em, c.P, c.T, c.x, ρl, ρv)
    end
    if isempty(bads)
        println("   (この u では範囲外のセルが無い。flash は範囲内でも、Shah の入力 (G, T_wall など) が原因の可能性)")
    else
        c0 = first(bads)
        println("\n(b) 最初のセル (", c0.comp, " #", c0.i, ") の、直近 ", RING, " 回の残差評価での履歴 (古い順):")
        println("      call        t[s]        ρm        Em          P[Pa]      T[K]   x")
        j0 = max(1, st.ncalls - RING + 1)
        for j in j0:st.ncalls
            k = mod1(j, RING)
            xr = ComponentArray(copy(st.ring_u[k]), pin.axes)
            cc = first(filter(c -> c.comp === c0.comp && c.i === c0.i, fluid_cells(pin, xr)))
            @printf("   %8d %10.4f %9.3f %11.1f %12.1f %8.2f %4.2f\n", j, st.ring_t[k], cc.ρ, cc.Em, cc.P, cc.T, cc.x)
        end
        println("   → 1 回の評価で ρm や P が飛んでいれば H1/H3 (Newton の試行点)、")
        println("     何十回もかけて少しずつ動いていれば H2/H4 (状態が本当にそちらへ進んでいる)")
    end
    println("\n(c) 失敗の直前の時系列 (", LOG_DT, " s ごと。mdot の単位は mg/s):")
    print_log_rows(last(st.log, N_PRINT))
    println("   → mdot が負になっていれば H2 (逆流)、cdP_max が P_max (", P_HI / 1e3, " kPa) に近づいていれば H1 (液の圧縮)")
    println("="^78, "\n")
end

function M.f_dae_residual!(res, du, u, p::DiagParams, t)
    st = p.st
    st.ncalls += 1
    k = mod1(st.ncalls, RING)
    st.ring_t[k] = t
    copyto!(st.ring_u[k], u)
    if t >= st.next_log
        push!(st.log, snapshot(p.inner, u, t))
        st.next_log = t + LOG_DT
    end
    try
        M.f_dae_residual!(res, du, u, p.inner, t)
    catch e
        e isa DomainError || rethrow()
        st.nerr += 1
        st.nerr == 1 && (st.t_first_err = t)
        if !st.reported
            st.reported = true
            report_failure(p.inner, st, u, t, e)
        end
        GUARD || rethrow()
        fill!(res, NaN)
    end
    return nothing
end

# ---------------------------- 1 ケースの実行 ----------------------------
function run_case(ΔT_startup)
    geo = make_geo(GEOMETRY)
    condition = M.LHPTransientCondition(
        fluid = "ammonia", T_cond = 15.0 + 273.15, T_amb = 25.0 + 273.15, T_init = 25.0 + 273.15,
        Q_apply_init = Q_INIT, Q_apply = Q_APPLY, t_step = T_STEP, τ_step = TAU_STEP,
        ΔT_startup = ΔT_startup, t_init = 0.0, t_end = T_END, saveat = 1.0)

    u0, du0, p = M.initial_condition(geo, condition, tables; verbose = false)
    n = length(getdata(u0))
    st = DiagState(0, condition.t_init, Any[], zeros(RING), [zeros(n) for _ in 1:RING], 0, false, NaN)
    pd = DiagParams(p, st)

    diffmask = Bool.(getdata(M.initialize_differential_vars(geo)))
    scales = M.variable_scales(geo, condition)
    abstol = M.absolute_tolerances(geo, condition; atol_frac = 1.0e-9, atol_frac_u = 1.0e-7)
    prob = M.build_dae_problem(pd, copy(getdata(u0)), du0, diffmask,
                               (condition.t_init, condition.t_end), scales;
                               linear_solver = :KLU, verbose = false,
                               guard_domain_errors = false)   # ガードは、このスクリプトの GUARD で切り替える
    solver = IDA(; linear_solver = :KLU, max_order = 2)
    empty!(st.log); st.next_log = condition.t_init; st.ncalls = 0    # パターン検出の呼び出しは除く
    st.nerr = 0; st.reported = false

    sol = try
        solve(prob, solver; abstol = abstol, reltol = 1.0e-4,
              tstops = M.applied_heat_breakpoints(condition), saveat = 1.0, maxiters = 10^7)
    catch e
        e isa DomainError ? nothing : rethrow()
    end

    rows = st.log
    mdots = [r.mdot for r in rows]
    i_min = isempty(mdots) ? 0 : argmin(mdots)
    i_neg = isempty(mdots) ? 0 : something(findfirst(<(-1.0e-9), mdots), 0)
    return (ΔT_startup = ΔT_startup,
            retcode = sol === nothing ? "例外で中断" : string(sol.retcode),
            t_final = sol === nothing ? (isempty(rows) ? NaN : rows[end].t) : sol.t[end],
            n_err = st.nerr, t_first_err = st.t_first_err,
            mdot_min = i_min == 0 ? NaN : mdots[i_min], t_mdot_min = i_min == 0 ? NaN : rows[i_min].t,
            t_first_neg = i_neg == 0 ? NaN : rows[i_neg].t,
            n_bad_max = isempty(rows) ? 0 : maximum(r.n_bad for r in rows))
end

results = [run_case(d) for d in STARTUP_LIST]

println("\n要約 (GEOMETRY=", GEOMETRY, ", Q ", Q_INIT, " → ", Q_APPLY, " W, GUARD=", GUARD, ")")
@printf("%10s %-22s %8s %6s %10s %12s %10s %10s %7s\n", "ΔT_start", "retcode", "t_end", "#err",
        "t_1st_err", "mdot_min[mg/s]", "t(mdot_min)", "t(mdot<0)", "#bad")
for r in results
    @printf("%10.2f %-22s %8.2f %6d %10.2f %12.4f %10.2f %10.2f %7d\n", r.ΔT_startup, r.retcode, r.t_final,
            r.n_err, r.t_first_err, r.mdot_min * 1e6, r.t_mdot_min, r.t_first_neg, r.n_bad_max)
end
println("  #bad: 記録した時刻のうち、範囲外の flash 結果を持つセル数の最大値 (GUARD=true のとき、落ちずに範囲外へ出た回数の目安)")