# =============================================================================
# residuals_flatcc.jl (平板型のCC)
#
# CC は、曲面を上、平面を下にした半円筒 (D_i = 27 mm)。気液は重力で分離し、蒸気は
# 曲面の最上部に中心角 θ_vap の円冠として溜まる (geometry_flatcc.jl)。
# CC流体の残差 (residual_cc!) は円筒型と同じものをそのまま使う (貯水池、N=1)。
# 円筒型との違いは次の3点。
#   1. 流体と接する「壁」は CCケースではなく secondary wick。secondary wick の温度
#      T_swick からの自然対流 (式70) が、CCケース → 流体 の自然対流の代わりになる。
#   2. 接触面積は、曲面 + 両端面 (S_inner) + 平面側 (D_i · dl)。平面側は常に液に
#      接する (蒸気は上)。両端面は、蒸気の円冠の分だけ蒸気側とする。
#   3. CCケース (T_wall.cc) は、流体ではなく secondary wick と熱伝導で熱をやり取りする
#      (Q_case_swick、residuals_flatwick.jl)。残差は residual_cc_wall! をそのまま使う。
#
# 【判断事項】
#   - 自然対流の代表長さ L は円筒型と同じ内径 D_i。平面側の液にも同じ Nu(Ra) を使う
#     (secondary wick が平面の下にあるので、実際は加熱面が下の対流になる)。
#   - 液管からの流入は secondary wick を通らず、CCに直接入る。
# =============================================================================

"""
    pipe_params(geo::FlatCCGeometry) は residuals_transport.jl 側で共通化済み。
    flatcc_params(geo::FlatCCGeometry) -> NamedTuple

CC流体の幾何量 (pipe_params に、半径 R と、自然対流の面積に使う量を加えたもの)。
N=1 のみ対応。S_total は、曲面 + 両端面 + 平面側の全面積。
"""
function flatcc_params(geo::FlatCCGeometry)
    geo.N == 1 || error("CCは N=1 (貯水池) のみ対応しています。geo.cc.N=$(geo.N)")
    return merge(pipe_params(geo),
                 (R = geo.R, A_flat = geo.D_i * geo.dl, S_total = geo.S_inner + geo.D_i * geo.dl))
end

"""
    flatcc_wall_params(geo::FlatCCGeometry) -> NamedTuple

CCケース (半円筒の壁、内径 D_swick、外径 D_o) の幾何量。V_w, A_cross は半円筒なので
円筒の半分。residuals_wall.jl の wall_ghost / wall_axial_conductance にそのまま渡せる。
"""
flatcc_wall_params(geo::FlatCCGeometry) = (
    N = geo.N, dl = geo.dl, D_i = geo.D_swick, D_o = geo.D_o,
    V_w = 0.5 * π * (geo.D_o^2 - geo.D_swick^2) / 4 * geo.dl,
    A_cross = 0.5 * π * (geo.D_o^2 - geo.D_swick^2) / 4,
)

"""
    flatcc_natural_convection(u_cc, T_swick, p, g, tables; δ_T=1e-3) -> NamedTuple

secondary wick (温度 T_swick) から、CC内の液と蒸気への熱を、気液別の自然対流で求める。

    Q_l, Q_v  … secondary wick → 液、→ 蒸気の熱 [W] (secondary wick が熱ければ正)
    Q_fluid   … Q_l + Q_v (流体が受け取る熱。CC流体と secondary wick の両方の残差で共有)
    geom, Ra_l, Ra_v, h_l, h_v

面積: 蒸気 A_vap = R θ_vap dl (曲面) + R²(θ_vap − sin θ_vap) (両端面の円冠2つ)、
      液   A_liq = S_total − A_vap (曲面の残り + 両端面の残り + 平面側)。
"""
function flatcc_natural_convection(u_cc, T_swick, p, g, tables::FluidTables; δ_T=1e-3)
    T = u_cc.T[1]
    sp = sat_props(T, tables)
    f_v = vapor_volume_fraction(u_cc.ρm[1], u_cc.x[1], sp.ρ_v)
    geom = cc_vapor_geometry(f_v, p.R, p.dl)
    θ = geom.θ_vap
    A_vap = p.R * θ * p.dl + p.R^2 * (θ - sin(θ))
    A_liq = p.S_total - A_vap

    ΔT_abs = soft_abs(T_swick - T, δ_T)
    L = p.D

    ν_l = sp.μ_l / sp.ρ_l
    α_l = sp.k_l / (sp.ρ_l * sp.cp_l)
    Ra_l = g * sp.β_l * ΔT_abs * L^3 / (ν_l * α_l)
    h_l = Nu_cc_smooth(Ra_l) * sp.k_l / L

    ν_v = sp.μ_v / sp.ρ_v
    α_v = sp.k_v / (sp.ρ_v * sp.cp_v)
    Ra_v = g * sp.β_v * ΔT_abs * L^3 / (ν_v * α_v)
    h_v = Nu_cc_smooth(Ra_v) * sp.k_v / L

    Q_l = h_l * A_liq * (T_swick - T)
    Q_v = h_v * A_vap * (T_swick - T)
    return (Q_l = Q_l, Q_v = Q_v, Q_fluid = Q_l + Q_v, geom = geom,
            Ra_l = Ra_l, Ra_v = Ra_v, h_l = h_l, h_v = h_v)
end