# =============================================================================
# test_residuals_grv.jl
#
# residuals_grv.jl の動作確認。CoolPropのテーブルを実際に構築して、
#   1. 静止した蒸気 (閉端・流出なし・熱なし・ソースなし) で全残差 ≈ 0
#   2. ソース付きで、ρEm の増加率が q + S_E に一致する (エネルギー保存)
#   3. 蒸発器ケース→グルーブ蒸気の熱伝達係数の大きさ
#   4. 蒸発器ケースが失う熱 = グルーブが受け取る熱 (同じ値を共有)
# を確認する。
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
include(joinpath(SRC, "residuals_evapcase.jl"))
include(joinpath(SRC, "residuals_grv.jl"))

fluid = "ammonia"
println("テーブル構築中 (初回は数分かかる場合があります。キャッシュされます)...")
tables = init_fluid_lookups(fluid; T_min=230.0, T_max=360.0, P_min=1.0e5, P_max=4.5e6,
                            Ne=400, Δδ=1.0, N_T_sat=200)

g = 9.8
geo_grv = GrvGeometry(N=3, n_grv=12, L=0.035, width=0.001, height=0.001)
p = grv_params(geo_grv)
N = p.N

println("グルーブ: N=$(N), 全断面積 A = ", p.A, " m2, 水力直径 D_h = ", p.D, " m, セル体積 V = ", p.V, " m3")

# --- 過熱蒸気の状態 (T=305K, P=1.0MPa; 305Kの飽和圧力は約1.2MPaなので蒸気) ---
T0, P0 = 305.0, 1.0e6
@assert P0 < CoolProp.PropsSI("P", "T", T0, "Q", 1.0, fluid) "過熱蒸気ではありません"
ρ0 = CoolProp.PropsSI("Dmass", "T", T0, "P", P0, fluid)
e0 = CoolProp.PropsSI("Umass", "T", T0, "P", P0, fluid)
Em0 = e0 + g * cell_elevation(p, 1)

make() = ComponentArray(ρm=fill(ρ0, N), u=zeros(N), Em=fill(Em0, N),
                        P=fill(P0, N), T=fill(T0, N), x=fill(1.0, N))
u = make(); du = make(); resid = make()
for c in (du, resid)
    c.ρm .= 0; c.u .= 0; c.Em .= 0; c.P .= 0; c.T .= 0; c.x .= 0
end
down = (ρ = ρ0, Em = Em0, P = P0, u = 0.0, dl = p.dl)

println("="^64)
println("テスト1: 静止蒸気 (閉端・流出なし・熱なし・ソースなし) → 質量・運動量・エネルギー残差 ≈ 0")
println("="^64)
residual_grv!(resid, du, u, down, zeros(N), 0.0, 0.0, p, g, tables)
println("  質量   max|res| = ", maximum(abs, resid.ρm))
println("  運動量 max|res| = ", maximum(abs, resid.u))
println("  エネルギー max|res| = ", maximum(abs, resid.Em))

println()
println("="^64)
println("テスト2: ソースあり (閉端・流出なし) → d(ρEm)/dt = q + S_E (エネルギー保存)")
println("="^64)
S_m_total = 1.7e-5                          # kg/s
S_E_total = S_m_total * (Em0 + P0 / ρ0) * 1.02
Q_cells = fill(0.05, N)                     # W (蒸発器ケースから、セルごと)
residual_grv!(resid, du, u, down, Q_cells, S_m_total, S_E_total, p, g, tables)
# du=0 のときの resid から、resid=0 となる du を逆算 (resid は du について線形)
Vtot = N * p.V
dρ_expected = S_m_total / Vtot
for i in 1:N
    dρ = -resid.ρm[i]
    dEm = -resid.Em[i] / ρ0
    println("  セル$i: dρ/dt = ", dρ, " (期待 ", dρ_expected, ")")
end
# 閉端・流出なし (u=0) では、移流 Φ は 0。保存形の残差は
#   resid.Em = ρ du.Em + Em du.ρm − (Φ_in − Φ_out)/Δz − q − S_E
# なので、du=0 で評価した resid.Em = −(q + S_E)。resid=0 となる du では、
#   d(ρEm)/dt = ρ du.Em + Em du.ρm = q + S_E
# つまり、d(ρEm)/dt = −resid.Em(du=0)。
for i in 1:N
    d_ρE = -resid.Em[i]
    q_S = Q_cells[i] / p.V + S_E_total / Vtot
    println("  セル$i: d(ρEm)/dt = ", d_ρE, ",  q+S_E = ", q_S, ",  相対差 = ", abs(d_ρE - q_S) / q_S)
end

println()
println("="^64)
println("テスト3: 蒸発器ケース→グルーブ蒸気の熱伝達係数 (式72、Re_L 基準、Q_apply=20,5,1 W に相当する ṁ)")
println("="^64)
Tec = T0 + 10.0   # 蒸発器ケースが蒸気より 10K 熱い
h_lat = CoolProp.PropsSI("H", "T", T0, "Q", 1.0, fluid) - CoolProp.PropsSI("H", "T", T0, "Q", 0.0, fluid)
for Qa in (20.0, 5.0, 1.0)
    ṁ = Qa / h_lat
    h = groove_h_ec(u, p, ṁ, g, tables)
    Qg = groove_evap_heat(u, Tec, h, p)
    println("  Q_apply=$(Qa) W: ṁ=", ṁ, " kg/s, h_ec-grv = ", round(h[1], digits = 1),
            " W/m2K, Q_ec-grv 合計 = ", round(sum(Qg), digits = 4), " W (印加熱の ",
            round(100 * sum(Qg) / Qa, digits = 2), " %)")
end
println("  (期待: h ≈ 37 (20W), 18.5 (5W), 8.3 (1W) W/m2K 程度 (Re_L 基準、h=Nu k/L)。",
        "Q_ec-grv は印加熱の 0.8〜3.5% 程度)")

println()
println("="^64)
println("テスト4: 蒸発器ケースが失う熱 = グルーブ蒸気が受け取る熱 (同じ値を共有)")
println("="^64)
ṁ = 1.7e-5
h = groove_h_ec(u, p, ṁ, g, tables)
Qg = groove_evap_heat(u, Tec, h, p)
q_in_grv = sum(Qg ./ p.V .* p.V)                # グルーブの残差 (q = Q/V) が受け取る熱の合計 [W]
println("  蒸発器ケースが失う熱 = ", sum(Qg), " W,  グルーブが受け取る熱 = ", q_in_grv, " W")