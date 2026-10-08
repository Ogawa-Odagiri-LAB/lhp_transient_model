# =============================================================================
# test_benchmark_residual.jl
#
# 残差評価 1 回のコストと、メモリ確保を、積分なしで、素早く測る (1〜2 分)。
# 長い時間積分を待たずに、性能改善 (型の安定化など) の効果を確認する。
#
# 背景: test_diagnose_fullscale.jl のプロファイルで、残差評価 1 回 (650 変数) は 1.136 ms、
# メモリ確保 589.8 KiB/回。時間の 85% は、セルごとの h_transport_wall と
# friction_force_transport (どちらも flash_ρe, sat_props, transport_props を呼ぶ)。
# NamedTuple の生成だけで 12% の自己時間があった。FluidTables の補間のフィールドが
# `Any` 型で、動的ディスパッチと箱詰めが起きていた疑い。FluidTables を型パラメータで
# 具体化して、確認する。
#
# 出力:
#   1. 残差評価 1 回のコスト (時間、メモリ確保、GC の割合)
#   2. 個々の物性関数 1 回あたりの時間とメモリ確保。
#      メモリ確保が 0 バイトなら、型が安定している (ヒープ確保なし)。
#   3. 型推論の確認 (戻り値の型が、具体的か)
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays,
      SparseArrays, LinearAlgebra, Random, Sundials, InteractiveUtils
using DiffEqBase: DAEProblem, DAEFunction, solve

const SRC = joinpath(@__DIR__, "..", "src")
for f in ("physics", "smoothing", "types", "proptables", "propaccess", "correlations",
          "heattransfer", "geometry_cc", "residuals_transport", "residuals_wall",
          "residuals_wick", "residuals_cc", "residuals_evapcase", "residuals_grv",
          "residuals_main", "initialcondition")
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
    t_init = 0.0, t_end = 100.0, saveat = 1.0)

println("FluidTables の型パラメータ数: ", length(typeof(tables).parameters),
        " (6 なら型パラメータ化されている。0 なら旧版)")
println("  P_of_ρe のフィールドの型は具体的か: ", isconcretetype(fieldtype(typeof(tables), :P_of_ρe)),
        ",  μ_l_sat: ", isconcretetype(fieldtype(typeof(tables), :μ_l_sat)),
        ",  e_range: ", isconcretetype(fieldtype(typeof(tables), :e_range)))

# 初期条件 (整合初期条件) の状態で、残差を評価する
u0, du0, p = initial_condition(geo, condition, tables; verbose = false)
n = length(getdata(u0))
u = copy(getdata(u0)); du = copy(du0); res = zeros(n)

println()
println("="^72)
println("1. 残差評価 1 回のコスト (変数 ", n, " 個、初期状態)")
println("="^72)
for _ in 1:20; f_dae_residual!(res, du, u, p, 0.0); end            # JIT のウォームアップ
st = @timed for _ in 1:1000; f_dae_residual!(res, du, u, p, 0.0); end
alloc = @allocated f_dae_residual!(res, du, u, p, 0.0)
println("  1 回あたり ", round(st.time / 1000 * 1e3, digits = 3), " ms (GC の割合 ",
        round(100 * st.gctime / st.time, digits = 1), " %), メモリ確保 ", round(alloc / 1024, digits = 1), " KiB/回")
println("  (改善前: 1.136 ms、589.8 KiB/回。ただし、これは、t=100 s の状態で測った値で、初期状態とは、")
println("   二相のセル数などが違うので、厳密な比較ではない。同じ状態での、改善前後の比較は、")
println("   proptables.jl を旧版に戻して、このスクリプトを実行すると、できる)")

println()
println("="^72)
println("2. 物性関数 1 回あたり (T=298 K、液・蒸気・二相の代表的な状態)")
println("="^72)
function bench(name, f, nrep = 200_000)
    f()                              # ウォームアップ
    st = @timed for _ in 1:nrep; f(); end
    a = @allocated f()
    println("  ", rpad(name, 34), lpad(round(st.time / nrep * 1e9, digits = 0), 8), " ns/回,  メモリ確保 ", lpad(a, 6), " バイト/回")
end
T0 = 298.0; ρ_liq = 604.0; e_liq = 4.5e5; ρ_vap = 8.0; e_vap = 1.5e6; ρ_2p = 305.0; e_2p = 4.75e5
bench("sat_props(T)", () -> sat_props(T0, tables))
bench("flash_ρe (液)", () -> flash_ρe(ρ_liq, e_liq, tables))
bench("flash_ρe (蒸気)", () -> flash_ρe(ρ_vap, e_vap, tables))
bench("flash_ρe (二相)", () -> flash_ρe(ρ_2p, e_2p, tables))
bench("transport_props (液)", () -> transport_props(ρ_liq, e_liq, tables))
bench("liquid_props(T, P)", () -> liquid_props(T0, 1.0e6, tables))
bench("Tsat(P)", () -> Tsat(1.0e6, tables))
bench("h_transport_wall (液)", () -> h_transport_wall(ρ_liq, 0.01, e_liq, 0.0, 296.0, 0.00435, ρ_liq * 0.01, 9.8, tables))
bench("h_transport_wall (蒸気)", () -> h_transport_wall(ρ_vap, 0.18, e_vap, 0.0, 296.0, 0.00435, ρ_vap * 0.18, 9.8, tables))
bench("friction_force_transport (液)", () -> friction_force_transport(ρ_liq, 0.01, e_liq, 0.0, 0.00435, ρ_liq * 0.01, 9.8, tables))
bench("friction_force_transport (蒸気)", () -> friction_force_transport(ρ_vap, 0.18, e_vap, 0.0, 0.00435, ρ_vap * 0.18, 9.8, tables))

println()
println("="^72)
println("3. 型推論 (戻り値の型が具体的か。Any があると、動的ディスパッチの疑い)")
println("="^72)
for (name, f, argtypes) in (("sat_props", sat_props, (Float64, typeof(tables))),
                            ("flash_ρe", flash_ρe, (Float64, Float64, typeof(tables))),
                            ("transport_props", transport_props, (Float64, Float64, typeof(tables))))
    rt = Base.return_types(f, argtypes)
    println("  ", rpad(name, 18), " → ", join(string.(rt), " | "))
end
println()
println("読み方: メモリ確保が 0 バイト、かつ、戻り値の型が NamedTuple{..., Tuple{Float64, ...}} なら、型が安定している。")
println("        Any や Tuple{Any, ...} が出たら、動的ディスパッチが残っている。")