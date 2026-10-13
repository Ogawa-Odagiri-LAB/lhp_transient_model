# =============================================================================
# residuals_grv.jl (グルーブ)
#
# グルーブは、ウィックの蒸発界面から蒸気が湧き出し、蒸発器ケースから
# 過熱され、蒸気管へ流れ出す、軸方向の流路。n_grv 本の溝を、全溝の断面積の
# 合計を持つ1本の管として扱い、輸送管の核 (residual_transport!) を、
# 質量・エネルギーのソースと、蒸発器ケースからの熱を与えて呼ぶ。
#
#   質量:       ∂ρm/∂t + ∂(ρm u)/∂z = S_m                       (S_m: ウィックの蒸発)
#   エネルギー: ρ ∂Em/∂t + ρ u ∂φ/∂z − (P/ρ)∂ρm/∂t = q + S_E − φ S_m
#                q   : 蒸発器ケースからグルーブ蒸気への熱 (式45)
#                S_E : 蒸発界面を通過するエネルギー (wick_groove_source)
#
# 【幾何】 A = n_grv A_grv_cross (全溝の合計)、水力直径 D_h = 2wh/(w+h)。
# 溝の断面は矩形 (幅 w、深さ h)。矩形ダクトの層流摩擦係数は、円管の 64/Re
# ではなく、アスペクト比に依存する (正方形で約 57/Re)。ここでは円管の式を、
# 水力直径で使っている (誤差は約12%)。グルーブの摩擦圧損は小さい
# (蒸気 ṁ=1.7e-5 kg/s で、溝1本の全長 35 mm で約 2 Pa) ので、影響は軽微。
#
# 【境界】 入口 (溝の端) は閉端 (u=0)。ゴーストは、セル1と同じ値にして、
# 入口側の勾配をゼロにする。出口は蒸気管の最初のセル
# (transport_ghost_down(u_vl, p_vl, p_grv) で作る)。
#
# 【蒸発器ケースからグルーブ蒸気への熱 Q_ec-grv (式45)】
#   Q_ec-grv,i = h_ec-grv,i (A_ec-grv / N) (T_ec − T_i)
# 過熱面は蒸発器ケースとグルーブ蒸気が接する上面だけ (A_ec_grv = width L n_grv)。
#
#   h_ec-grv は、ドキュメントの式72 (水平平板の強制対流) に従う:
#     Re_L = G L / μ_v,   G = (ṁ_evap/n_grv)/A_grv_cross
#     Nu_L = 0.664 Re_L^(1/2) Pr^(1/3)    (Re_L < 5e5)
#          = 0.0296 Re_L^(4/5) Pr^(1/3)   (それ以外)
#     h = Nu_L k_v / L                    (L: 溝の長さ。Re と h の代表長さを L に統一)
#   ṁ_evap はウィックの蒸発量 (全体)。溝1本あたりに直すため n_grv で割る。
#   物性 (μ_v, k_v, Pr) はグルーブのセルの状態から取る (ρ は Re の中で相殺)。
#
# 【旧コードとの違い】
#   - 旧 calc_h_ec_gr は、Re を溝の深さ (height) 基準、h を Nu k/L (溝の長さ L で
#     割る) としており、代表長さが混在していた (h が約1/6 になる)。
#     本実装は、ドキュメントに従い、Re も h も溝の長さ L に統一する。
#   - 旧は Re≤2300 で層流/乱流を切り替え、乱流側は 0.0371 Pr Re^0.8 だった。
#     本実装は式72 (遷移 Re=5e5、乱流側 0.0296 Re^0.8 Pr^(1/3))。
#   - 旧 steadyfuncs.jl は溝1本あたりの流量、過渡側は 12本まとめて計算していた。
#     本実装は、全体の ṁ_evap を n_grv で割って溝1本あたりにする。
#
# 【判断事項】
#   1. ウィックの蒸発界面の圧力 P_grv は、グルーブの圧力の平均
#      (groove_pressure)。ウィックは半径方向の1次元で、軸方向の位置を
#      区別しないため。
# =============================================================================

# -----------------------------------------------------------------------------
# 1. 幾何
# -----------------------------------------------------------------------------

"""
    grv_params(geo::Union{GrvGeometry,FlatGrvGeometry}) -> NamedTuple

グルーブの幾何量。pipe_params と同じフィールド (N, dl, D, A, V, θ, elev0) を
持つので、residual_transport! や transport_ghost_up/down にそのまま渡せる。
A, V は全溝の合計。D は矩形溝の水力直径。
"""
grv_params(geo::Union{GrvGeometry,FlatGrvGeometry}) = (
    N = geo.N, dl = geo.dl, L = geo.L,
    D = 2 * geo.width * geo.height / (geo.width + geo.height),
    A = geo.n_grv * geo.A_grv_cross,
    V = geo.n_grv * geo.A_grv_cross * geo.dl,
    θ = geo.θ, elev0 = geo.elev0,
    n_grv = geo.n_grv, A_one = geo.A_grv_cross, A_ec_grv = geo.A_ec_grv,
)

"""
    groove_pressure(u_grv) -> P_grv

ウィックの蒸発界面温度 T_sat(P_grv) に使う、グルーブの代表圧力
(セルの圧力の平均)。N=1 なら、そのセルの圧力。
"""
groove_pressure(u_grv) = sum(u_grv.P) / length(u_grv.P)

# -----------------------------------------------------------------------------
# 2. 蒸発器ケースからグルーブ蒸気への熱
# -----------------------------------------------------------------------------

"""
    groove_h_ec(u, p, ṁ_evap, g, tables; δ_m=1e-9, ε_Re=2e4) -> h [W/m²K] (セルごと)

蒸発器ケース→グルーブ蒸気の熱伝達係数 (式72、水平平板の強制対流)。
  Re_L = G L/μ_v,  G = (|ṁ_evap|/n_grv)/A_grv_cross
  h = Nu_L k_v / L   (Nu_L は Nu_flatplate_smooth: 遷移 Re=5e5)
μ_v, k_v, Pr はグルーブのセルの状態 (ρm, Em, u) から transport_props で取る。
ṁ_evap の符号は soft_abs で処理する (δ_m は ṁ の典型値 1e-5 kg/s より十分小さく)。
入口は閉端 (u=0)。
"""
function groove_h_ec(u, p, ṁ_evap, g, tables::FluidTables; δ_m=1.0e-9, ε_Re=2.0e4)
    h = Vector{Float64}(undef, p.N)
    G = soft_abs(ṁ_evap, δ_m) / (p.n_grv * p.A_one)
    for i in 1:p.N
        u_c = 0.5 * ((i == 1 ? 0.0 : u.u[i-1]) + u.u[i])
        em = internal_energy(u.Em[i], u_c, cell_elevation(p, i), g)
        tp = transport_props(u.ρm[i], em, tables)
        Re = G * p.L / tp.μ
        h[i] = h_ec_grv(Re, tp.Pr, p.L, tp.k; ε_Re = ε_Re)
    end
    return h
end

"""
    groove_evap_heat(u, T_ec, h, p) -> Q [W] (セルごと)

蒸発器ケースからグルーブの各セルの蒸気へ渡す熱 (式45)。蒸発器ケースが
熱ければ正。蒸発器ケースの残差 (residual_evapcase) では sum(Q) を失う熱と
して、グルーブの残差 (residual_grv!) ではこの値を熱源として使う
(同じ値を共有するので、エネルギーが厳密に一致する)。
"""
groove_evap_heat(u, T_ec, h, p) = [h[i] * (p.A_ec_grv / p.N) * (T_ec - u.T[i]) for i in 1:p.N]

# -----------------------------------------------------------------------------
# 3. グルーブの残差
# -----------------------------------------------------------------------------

"""
    residual_grv!(resid, du, u, down, Q_ec_grv, S_m_total, S_E_total, p, g, tables; kwargs...)

グルーブ (N セル) の流体残差を resid に書き込む。resid, du, u は
グルーブの ComponentVector (ρm,u,Em,P,T,x)。

  down       : 下流 (蒸気管の最初のセル) のゴースト値。
               transport_ghost_down(u_vl, p_vl, p) で作る
  Q_ec_grv   : groove_evap_heat の戻り値 [W] (セルごとのベクトル)
  S_m_total  : ウィックから蒸発して入る質量流量 [kg/s] (wick_groove_source の第1戻り値)
  S_E_total  : ウィックから入るエネルギー [W] (wick_groove_source の第2戻り値)
  p          : grv_params(geo.grv)

S_m, S_E は、全セルに一様に分配する (蒸発は溝の全長で一様に起きるとする)。
kwargs は residual_transport! へ渡される (ε_u, δ_u, ε_Re, ε_phase)。
"""
function residual_grv!(resid, du, u, down, Q_ec_grv, S_m_total, S_E_total, p, g,
    tables::FluidTables; kwargs...)

    N = p.N
    # 入口は閉端: u=0、ゴーストはセル1と同じ値 (入口側の勾配をゼロにする)
    up = (ρ = u.ρm[1], Em = u.Em[1], P = u.P[1], u = 0.0, dl = p.dl)

    V_total = N * p.V
    S_m = fill(S_m_total / V_total, N)
    S_E = fill(S_E_total / V_total, N)
    q = Q_ec_grv ./ p.V

    return residual_transport!(resid, du, u, up, down, q, p, g, tables;
                               S_m = S_m, S_E = S_E, kwargs...)
end