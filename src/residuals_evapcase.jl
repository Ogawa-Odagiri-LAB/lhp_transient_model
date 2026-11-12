# =============================================================================
# residuals_evapcase.jl (蒸発器ケース)
#
# 設計ドキュメント 式41-45 に対応する。蒸発器ケースは、軸方向に分割しない
# 集中容量の壁 (T_wall.evap、スカラー) で、微分変数のみ (代数変数なし)。
#
#   C_ec dT_ec/dt = Q_apply − (Q_ec-vl + Q_ec-cc + Q_ec-wick + Q_ec-grv)   (式41)
#
#   Q_apply   … 印加熱量 (ヒーター) [W]
#   Q_ec-vl   … 蒸気管の最初の壁セルへの軸方向熱伝導 [W]
#   Q_ec-cc   … CCケースへの軸方向熱伝導 [W]
#   Q_ec-wick … ウィックの蒸発界面への熱 (式44)。wick_interface が計算した値 iface.Q_ec
#   Q_ec-grv  … グルーブ蒸気への熱 (式45)。groove_evap_heat (residuals_grv.jl) が
#               セルごとに計算した値の合計。グルーブ側の熱源にもなるので、
#               ここでは呼び出し側が計算した値を受け取る
#
# 【共有する量】 Q_ec-wick (= wick_interface の iface.Q_ec)、Q_ec-vl、
# Q_ec-cc、Q_ec-grv は、いずれも別の残差 (ウィック・蒸気管の壁・CCの壁・
# グルーブ) でも使うので、全体の残差関数 (residuals_main.jl) が状態から
# 1回だけ計算して、各残差に配る。蒸発器ケース側で別々に計算し直すと、
# 失う熱と受け取る熱がずれて、エネルギーが保存しなくなる。
#
# 【Q_ec-vl, Q_ec-cc (式42, 43) の扱い】
# ドキュメントの式42, 43 は、径方向の円筒殻の熱伝導の形
#     Q = k 2πL / ln(Do/Di) (T_ec − T_vl)
# だが、本実装は既存コードの実装 (軸方向の熱伝導) に従う:
#     Q = G_ave (T_ec − T_wall,neighbor)
#     G_ave = [1/G_ec + 1/G_neighbor]⁻¹,   G = k π(Do²−Di²)/4 / (l/2)
# つまり、各壁の環状断面積 A=π(Do²−Di²)/4 を通る軸方向熱伝導で、
# 各壁セルの「半長さ (l/2)」ずつを直列に結合する。蒸発器ケースは長さ L
# (geo.grv.L)、蒸気管の最初のセルは dl、CCケースは長さ L (geo.cc.L)。
# これは residuals_wall.jl の wall_axial_conductance と同じ計算なので、
# 蒸気管・CCの壁の残差 (residual_wall_line!, residual_cc_wall!) が受け取る熱と、
# ここで蒸発器ケースが失う熱が、厳密に一致する (同じ関数・同じ引数で
# 計算するため)。ドキュメントの式42, 43 は、この形に直す必要がある。
#
# 【旧コードからの修正】
#   - 熱容量 C_ec = V ρ c の V に π が抜けていた (旧: (D_o²−D_i²)×0.25×L)。
#     本実装は V = π(D_o²−D_i²)/4 × L。
#   - 旧 calc_Q_ec_cc は T_cc に T_wall.evap を代入しており、常に Q_ec-cc=0 だった。
#     本実装は CCケースの壁温度 (T_wall.cc) を使う。
#
# 【蒸発器ケースの幾何】 GrvGeometry の D_i, D_o, L (ケースの内径・外径・長さ)。
# =============================================================================

# -----------------------------------------------------------------------------
# 1. 幾何・熱容量
# -----------------------------------------------------------------------------

"""
    evap_case_params(geo::GrvGeometry) -> NamedTuple

蒸発器ケース (集中容量の壁) の幾何量。residuals_wall.jl の wall_ghost /
wall_axial_conductance に渡せるよう、A_cross (環状断面積) と dl (軸方向長さ) を
持たせている (蒸気管の最初の壁セルの上流側ゴーストとして使う)。
"""
evap_case_params(geo::GrvGeometry) = (
    D_i = geo.D_i, D_o = geo.D_o, dl = geo.L,
    A_cross = π * (geo.D_o^2 - geo.D_i^2) / 4,
    V_w = π * (geo.D_o^2 - geo.D_i^2) / 4 * geo.L,
)

"""
    evap_case_heat_capacity(T, ec) -> C_ec [J/K]

蒸発器ケースの熱容量 (SUS304、physics.jl の ρ_sus304, c_sus304)。
"""
evap_case_heat_capacity(T, ec) = ec.V_w * ρ_sus304(T) * c_sus304(T)

# -----------------------------------------------------------------------------
# 2. 各熱流 (蒸発器ケースが失う向きを正)
# -----------------------------------------------------------------------------

"""
    Q_ec_wall(T_ec, ec, T_nb, wp_nb) -> Q [W]

蒸発器ケースから、隣接する壁 (蒸気管の最初のセル、またはCCケース) へ
軸方向熱伝導で出ていく熱 (蒸発器ケースが熱ければ正)。
wp_nb は隣接する壁の (A_cross, dl) を持つ NamedTuple
(蒸気管なら wall_params(geo.vl)、CCなら cc_wall_params(geo.cc))。

wall_axial_conductance は対称 (引数の順序によらない) なので、
隣接する壁の残差側が計算する「受け取る熱」と、符号を除いて一致する。
"""
Q_ec_wall(T_ec, ec, T_nb, wp_nb) = wall_axial_conductance(T_ec, ec, T_nb, wp_nb) * (T_ec - T_nb)

# -----------------------------------------------------------------------------
# 3. 蒸発器ケースの残差
# -----------------------------------------------------------------------------

"""
    residual_evapcase(du_Tec, T_ec, ec, Q_apply, Q_ec_vl, Q_ec_cc, Q_ec_wick, Q_ec_grv) -> 残差

蒸発器ケースの残差 (式41)。他の残差関数 (residual_*!) と違って、値を返す
(! なし)。蒸発器ケースは集中容量で1要素しかないため、配列を受けて書き込む形に
する必要がないので。呼び出し側で resid.T_wall.evap[1] に代入する
(T_wall.evap は types.jl で長さ1の配列)。

    residual = C_ec dT_ec/dt − (Q_apply − Q_ec_vl − Q_ec_cc − Q_ec_wick − Q_ec_grv)

  du_Tec    : du.T_wall.evap[1] (dT_ec/dt)
  T_ec      : u.T_wall.evap[1]

  Q_apply   : 印加熱量 [W]
  Q_ec_vl   : Q_ec_wall(T_ec, ec, T_vl1, wall_params(geo.vl)) の値
  Q_ec_cc   : Q_ec_wall(T_ec, ec, T_ccw, cc_wall_params(geo.cc)) の値
  Q_ec_wick : wick_interface の iface.Q_ec (式44)
  Q_ec_grv  : グルーブ蒸気への熱 (式45)

"""
function residual_evapcase(du_Tec, T_ec, ec, Q_apply, Q_ec_vl, Q_ec_cc, Q_ec_wick, Q_ec_grv; C_extra = 0.0)
    C = evap_case_heat_capacity(T_ec, ec) + C_extra
    return C * du_Tec - (Q_apply - Q_ec_vl - Q_ec_cc - Q_ec_wick - Q_ec_grv)
end