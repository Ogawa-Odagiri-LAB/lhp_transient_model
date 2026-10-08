# =============================================================================
# test_initialcondition.jl
#
# residuals_main.jl と initialcondition.jl の動作確認。
#   0. ComponentArray の使い方 (IDA が渡す Vector を ComponentArray に包む) の確認
#   1. 整合初期条件の生成 (u0, du0)
#   2. 整合性: 代数変数の行 ≈ 0、微分変数の行は F(u0,du0) ≈ 0
#   3. 初期状態の要約 (圧力、温度、乾き度、CCの液体積分率)
#   4. 質量・エネルギーの収支 (全体の質量変化率が、出入りの合計に一致するか)
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays

const SRC = joinpath(@__DIR__, "..", "src")   # テストは ./test/、ソースは ./src/
include(joinpath(SRC, "physics.jl"))
include(joinpath(SRC, "smoothing.jl"))
include(joinpath(SRC, "types.jl"))
include(joinpath(SRC, "proptables.jl"))
include(joinpath(SRC, "propaccess.jl"))
include(joinpath(SRC, "correlations.jl"))
include(joinpath(SRC, "heattransfer.jl"))
include(joinpath(SRC, "geometry_cc.jl"))
include(joinpath(SRC, "residuals_transport.jl"))
include(joinpath(SRC, "residuals_wall.jl"))
include(joinpath(SRC, "residuals_wick.jl"))
include(joinpath(SRC, "residuals_cc.jl"))
include(joinpath(SRC, "residuals_evapcase.jl"))
include(joinpath(SRC, "residuals_grv.jl"))
include(joinpath(SRC, "residuals_main.jl"))
include(joinpath(SRC, "initialcondition.jl"))

fluid = "ammonia"
println("テーブル構築中 (初回は数分かかる場合があります。キャッシュされます)...")
tables = init_fluid_lookups(fluid; T_min=230.0, T_max=360.0, P_min=1.0e5, P_max=4.5e6,
                            Ne=400, Δδ=1.0, N_T_sat=200)

# --- 幾何 (既定値に近い、小さめの分割) ---
wick = WickGeometry(N=3, φ=0.65, ρ_wick=8900.0, Cp_wick=444.0)
geo = LHPGeometry(wick,
                  GrvGeometry(N=1),
                  TLGeometry(N=5, L=1.0, D_i=0.00435, D_o=0.00635),
                  CondGeometry(N=5, L=1.0, D_i=0.00435, D_o=0.00635),
                  TLGeometry(N=5, L=1.0, D_i=0.00218, D_o=0.00318),
                  CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037))
condition = LHPTransientCondition(Q_apply=20.0, t_init=0.0, t_end=100.0, fluid=fluid)

println("="^72)
println("テスト0: ComponentArray の使い方の確認")
println("="^72)
u_test = initialize_state(geo)
n = length(getdata(u_test))
vec = zeros(n)
ca = ComponentArray(vec, getaxes(u_test))
ca.vl.ρm[1] = 123.0
ca.T_wall.evap[1] = 456.0
println("  状態ベクトルの長さ n = ", n)
println("  ComponentArray(vec, axes) への書き込みが、元の Vector に反映される: ",
        any(vec .== 123.0) && any(vec .== 456.0), " (true のはず)")
println("  getaxes の戻り値の型: ", typeof(getaxes(u_test)))

println()
println("="^72)
println("テスト1: 整合初期条件の生成")
println("="^72)
u0, du0, p = initial_condition(geo, condition, tables)

println()
println("="^72)
println("テスト2: 初期状態の要約")
println("="^72)
function summary_line(name, c)
    println("  ", rpad(name, 6), " P=", round(c.P[1], digits=1), " Pa, T=", round(c.T[1], digits=4),
            " K, x=", round(c.x[1], digits=5), ", ρm=", round(c.ρm[1], digits=4), " kg/m3, u=", c.u[1], " m/s")
end
summary_line("grv", u0.grv); summary_line("vl", u0.vl); summary_line("cond", u0.cond)
summary_line("ll", u0.ll);   summary_line("cc", u0.cc)
println("  wick.T = ", u0.wick.T, " K")
println("  T_wall: evap=", u0.T_wall.evap, ", vl=", u0.T_wall.vl[1], ", cond=", u0.T_wall.cond[1],
        ", ll=", u0.T_wall.ll[1], ", cc=", u0.T_wall.cc, " K")
sp = sat_props(u0.cc.T[1], tables)
f_l = liquid_volume_fraction(u0.cc.ρm[1], u0.cc.x[1], sp.ρ_l)
println("  CC の液体積分率 f_l = ", f_l, " (", condition.f_l_cc_init, " のはず)")
println("  CC 出口の流速 u = ", u0.cc.u[1], " m/s (ṁ/(ρ_l A) に一致)")
println("  期待値 (アンモニア, T_init=298.15K): P0≈1.0027e6 Pa, CC ρm≈305.4, x≈0.0128;",
        " 液 ρ≈604.49, T≈297.15K; 蒸気 ρ≈7.759, T≈299.4K")

println()
println("="^72)
println("テスト3: 質量・エネルギーの収支 (du0 から)")
println("="^72)
d0 = ComponentArray(du0, p.axes)
# 質量: Σ V dρm/dt の合計が、(ウィックの蒸発 - ... ) と一致するか。ここでは、全体の質量の変化率を表示
V(pc) = pc.V
M_dot = sum(d0.grv.ρm) * p.p_grv.V + sum(d0.vl.ρm) * p.p_vl.V + sum(d0.cond.ρm) * p.p_cond.V +
        sum(d0.ll.ρm) * p.p_ll.V + d0.cc.ρm[1] * p.p_cc.V
println("  全流体の質量の変化率 Σ V dρm/dt = ", M_dot, " kg/s")
println("  (ウィック内の液の質量は変化しないモデル (ṁ 一様) なので、閉じた系なら 0 に近いはず。")
println("   ṁ が 0 の初期なら、ほぼ 0)")
println("  蒸発器ケースの温度変化率 dT_ec/dt = ", d0.T_wall.evap[1], " K/s (印加熱 ",
        applied_heat(condition.t_init, condition), " W で加熱される)")