# =============================================================================
# residuals_flatevapcase.jl (平板型の蒸発器ケース)
#
# 蒸発器ケースは、円筒型と同じく集中容量の壁 (T_wall.evap)。残差は residual_evapcase
# (式41) をそのまま使い、幾何と熱伝導の経路だけ平板用にする。
#
#   V_w = W_case L_case H_case − W L H      (ケースのブロックから primary wick の領域を除く)
#
# 隣の壁への熱伝導 (どちらも residuals_wall.jl の wall_axial_conductance を使うので、
# 相手の壁の残差が受け取る熱と厳密に一致する):
#   ・蒸気管の最初の壁 (軸方向): 面積 A_cross_ax = W_case H_case − W H、
#     経路 L_case/2 (ケース側)。ec = (A_cross = A_cross_ax, dl = L_case)。
#   ・CCケース (上方向): ケースの壁の接合部の断面積 A_cross_h = W_case L_case − W L。
#     ケース側の経路 H_case/2、CCケース側の経路 D_o/4 (半円筒の壁の中心までの目安)。
#     ec_h = (A_cross = A_cross_h, dl = H_case)、cw_h = (A_cross = A_cross_h, dl = D_o/2)。
#     [判断事項: 接合部の面積と経路長は、断面図からの目安。実測・CADで見直せる。]
# primary wick の側面への入熱は無視する (ec → wick は界面 Q_ec のみ)。
# =============================================================================

"""
    flat_evap_case_params(geo::FlatGrvGeometry) -> NamedTuple

蒸発器ケースの幾何量 (V_w, A_cross, dl は、円筒型の evap_case_params と同じ意味)。
"""
flat_evap_case_params(geo::FlatGrvGeometry) = (
    D_i = 0.0, D_o = 0.0,
    dl = geo.L_case,
    A_cross = geo.W_case * geo.H_case - geo.W * geo.H,
    V_w = geo.W_case * geo.L_case * geo.H_case - geo.W * geo.L * geo.H,
)

"""
    flat_ec_cc_params(geo_grv::FlatGrvGeometry, geo_cc::FlatCCGeometry) -> (ec_h, cw_h)

蒸発器ケース ↔ CCケースの熱伝導の経路 (上の説明)。
Q_ec_cc = Q_ec_wall(T_ec, ec_h, T_ccw, cw_h) で、wall_axial_conductance を使う。
"""
function flat_ec_cc_params(geo_grv::FlatGrvGeometry, geo_cc::FlatCCGeometry)
    A_h = geo_grv.W_case * geo_grv.L_case - geo_grv.W * geo_grv.L
    return (A_cross = A_h, dl = geo_grv.H_case), (A_cross = A_h, dl = geo_cc.D_o / 2)
end