# =============================================================================
# test_solver_sweep.jl
#
# IDA の設定 (BDF の最大次数、絶対許容誤差、相対許容誤差) を変えて、同じ計算を
# 繰り返し、計算時間とステップ数、解の差を比較する。
#
# 背景: エネルギー式を保存形に直した後も、IDA の刻みが小さい。これまでの結果:
#   ・次数 5 の刻みは、精度でなく、安定性で決まる (rtol を 1e-4 → 1e-6 と 100 倍厳しく
#     しても、15 s で 48〜57 万ステップのまま。平均 約 30 μs)。BDF の 3〜5 次は、減衰の
#     弱い振動 (音波) に対して、安定領域が狭い。
#   ・次数 2 (A 安定) は、解をほぼ変えずに、ステップを大きく減らす (前回: 約 13 分の 1)。
#     残る問題は、Newton 法の収束失敗 (ステップの 4 割) と、そのたびのヤコビアンの再計算。
#   ・KLU (疎行列) は、Dense と一致し、2.8 倍速い。
#
# 【saveat の影響の疑い】 前回のスイープ (T_END=5, saveat=0.05) の次数 5 の基準は、
# 棄却 63 062、ヤコビアン 98 169 で、同じ設定の 15 s の計算 (saveat=0.5: 棄却 6 908、
# ヤコビアン 17 778) より、はるかに多かった (15 s の中に 5 s が含まれるので、あり得ない)。
# 保存時刻ごとに IDA が停止・再開して、BDF の履歴が失われ、次数 1・小さい刻みから
# 立ち上げ直している可能性がある (未確認)。このスクリプトは、設定ごとに saveat を
# 指定でき、同じ設定で saveat だけ変えた組を入れて、直接確かめる。
#
# 【T_END=15 s】 蒸気の前線が凝縮器に届く t≈7〜11 s を含める。起動直後の数秒に、刻みの
# 縮小が集中するので、短い区間では、準定常での 1 秒あたりのコストが分からない。
#
# 基準 (1番目) = 次数 5、rtol 1e-6、atol_u 1e-9 (test_solver.jl と同じ。KLU で約 124 s)。
# 各設定の、基準との解の差、時間、ステップ数、収束失敗を比べる。
# 実行時間の目安: 1設定あたり 2〜10 分。全体で 30 分前後。他の計算を同時に動かさないこと。
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays,
      SparseArrays, LinearAlgebra, Random, Sundials, Statistics
using DiffEqBase: DAEProblem, DAEFunction, solve

# ---------------------------- 設定 ----------------------------
const T_END = 15.0
const LINEAR_SOLVER = :KLU        # :Dense または :KLU
# (名前, reltol, atol_frac, atol_frac_u, max_order, max_nonlinear_iters, saveat)
const CONFIGS = [
    ("基準: 次数5, rtol=1e-6, atol_u=1e-9",            1.0e-6, 1.0e-9, 1.0e-9, 5, 4,  0.5),
    ("次数2, rtol=1e-4, atol_u=1e-7",                  1.0e-4, 1.0e-9, 1.0e-7, 2, 4,  0.5),
    ("次数2, 同上, saveat=0.05 (saveat の影響)",         1.0e-4, 1.0e-9, 1.0e-7, 2, 4,  0.05),
    ("次数2, Newton 反復の上限 8",                      1.0e-4, 1.0e-9, 1.0e-7, 2, 8,  0.5),
    ("次数2, Newton 反復の上限 15",                     1.0e-4, 1.0e-9, 1.0e-7, 2, 15, 0.5),
    ("次数3, 上限 8",                                   1.0e-4, 1.0e-9, 1.0e-7, 3, 8,  0.5),
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
println("テーブル構築中 (初回は数分かかる場合があります。キャッシュされます)...")
tables = init_fluid_lookups(fluid; T_min=230.0, T_max=360.0, P_min=1.0e5, P_max=4.5e6,
                            Ne=400, Δδ=1.0, N_T_sat=200)

wick = WickGeometry(N=3, φ=0.65, ρ_wick=8900.0, Cp_wick=444.0)
geo = LHPGeometry(wick,
                  GrvGeometry(N=1),
                  TLGeometry(N=5, L=1.0, D_i=0.00435, D_o=0.00635),
                  CondGeometry(N=5, L=1.0, D_i=0.00435, D_o=0.00635),
                  TLGeometry(N=5, L=1.0, D_i=0.00218, D_o=0.00318),
                  CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037))
make_condition(saveat) = LHPTransientCondition(Q_apply=20.0, t_init=0.0, t_end=T_END, saveat=saveat, fluid=fluid)

# 注目する量 (最終時刻の値)
function summarize(sol, p)
    c_end = components_at(sol, p, length(sol.t))
    ρ_l = tables.ρ_l_sat(c_end.cc.T[1])
    ṁ = c_end.cc.u[1] * ρ_l * p.p_cc.A
    return (T_ec = c_end.T_wall.evap[1], P_cc = c_end.cc.P[1], T_cc = c_end.cc.T[1], ṁ = ṁ,
            P_grv = c_end.grv.P[1], T_cond1 = c_end.T_wall.cond[1])
end

stat_names = ("nf", "njacs", "naccept", "nreject", "nnonlinconvfail")
function get_stats(sol)
    st = sol.stats
    names = fieldnames(typeof(st))
    g(n) = n in names ? getfield(st, n) : -1
    return Dict(n => g(Symbol(n)) for n in stat_names)
end

results = []
for (iconf, (name, rt, af, afu, mo, mni, sv)) in enumerate(CONFIGS)
    println("="^76)
    println("設定: ", name)
    println("="^76)
    local sol, p, t_wall
    try
        t_wall = @elapsed (sol, p, _) = solve_transient_lhp!(geo, make_condition(sv), tables;
            linear_solver=LINEAR_SOLVER, reltol=rt, atol_frac=af, atol_frac_u=afu,
            max_order=mo, max_nonlinear_iters=mni, verbose=(iconf == 1))   # 1番目だけ統計を全表示
        ok = occursin("Success", string(sol.retcode))
        s = summarize(sol, p)
        st = get_stats(sol)
        push!(results, (name = name, ok = ok, t = t_wall, s = s, st = st))
        println("  retcode=", sol.retcode, ", 計算時間 ", round(t_wall, digits=1), " s, ",
                "受理ステップ ", st["naccept"], " (平均刻み ", round(T_END / max(st["naccept"], 1) * 1e6, digits=1), " μs), ",
                "棄却 ", st["nreject"], ", ヤコビアン ", st["njacs"], ", 収束失敗 ", st["nnonlinconvfail"])
        println("  最終: T_ec=", round(s.T_ec, digits=4), " K, P_cc=", round(s.P_cc, digits=1),
                " Pa, T_cc=", round(s.T_cc, digits=4), " K, ṁ=", round(s.ṁ, sigdigits=4), " kg/s")
    catch err
        println("  失敗: ", sprint(showerror, err)[1:min(300, end)])
        push!(results, (name = name, ok = false, t = NaN, s = nothing, st = nothing))
    end
end

println()
println("="^76)
println("まとめ (基準 = 1番目との比較)")
println("="^76)
ref = results[1]
println(rpad("設定", 46), lpad("時間[s]", 9), lpad("ステップ", 10), lpad("棄却", 8), lpad("収束失敗", 9),
        lpad("ヤコビアン", 11), lpad("ΔT_ec[K]", 10), lpad("ΔP_cc[Pa]", 10), lpad("Δṁ/ṁ", 9))
for r in results
    if r.s === nothing || ref.s === nothing
        println(rpad(r.name, 46), "  (失敗)")
        continue
    end
    println(rpad(r.name, 46), lpad(round(r.t, digits=1), 9), lpad(r.st["naccept"], 10), lpad(r.st["nreject"], 8),
            lpad(r.st["nnonlinconvfail"], 9), lpad(r.st["njacs"], 11),
            lpad(round(r.s.T_ec - ref.s.T_ec, digits=4), 10), lpad(round(r.s.P_cc - ref.s.P_cc, digits=1), 10),
            lpad(round((r.s.ṁ - ref.s.ṁ) / ref.s.ṁ, sigdigits=2), 9))
end
println("\n(解の差が小さい (ΔT_ec < 0.01 K, ΔP_cc < 100 Pa, Δṁ/ṁ < 1%) まま、時間とステップ数が減る設定を採用する)")
println("(2番目と3番目は、saveat だけが違う。ステップ・棄却・ヤコビアンが大きく違えば、saveat が刻みを歪めている)")