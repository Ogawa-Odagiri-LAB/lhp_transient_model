# =============================================================================
# residuals_cc.jl (CC: コンペンセーションチャンバー)
#
# CC内は貯水池のようなもので、流れがほとんどない (20 W のアンモニアで
# ṁ≈1.7e-5 kg/s、内径35 mm → u≈3e-5 m/s)。輸送管と同じ式
# (質量・運動量・エネルギー) をそのまま使うと、次の問題がある。
#   1. 出口 (CC→ウィック) の面の運動量式に意味がない。ウィックの圧力損失は
#      モデルに無く、この面の流量は蒸発量 ṁ (式48) で決まる。
#   2. 流速が風上ブレンドの幅 ε_u (液管では1e-4 m/s) より小さく、流れの向きが
#      判定できない (w≈0.5〜0.64)。エネルギー式の移流項の36〜50%が、
#      下流側 (ウィック) の値を使う前進差分になってしまう。
#      ε_u をCC用に小さくしても、流量が変わるたびに調整が必要になる。
# そこで、CCは貯水池 (0次元の集中定数) として扱い、式1, 3, 4〜6 を
# CCの体積全体で積分した形を使う。運動量式 (式2) は使わない。
#
#   質量:       V dρm/dt = ṁ_in − ṁ_out
#   エネルギー: ρ V dEm/dt = ṁ_in(φ_in − Em) − ṁ_out(φ_out − Em) + Q
#   状態拘束:   P, T, x = (ρm, em) から、em = Em − g z  (運動エネルギーは無視)
#   出口の流速: u_out = ṁ_out/(ρ_l A)   (代数拘束。運動量式の代わり)
#
# エネルギー式は、CCの全エネルギー ρ V Em の保存則
#   d(ρ V Em)/dt = ṁ_in φ_in − ṁ_out φ_out + Q
# から、質量の式を使って整理したもの (保存形)。
#
# 【CCの状態変数】 輸送管と同じ6変数 (ρm,u,Em,P,T,x) をそのまま使う
# (types.jl)。ただし u は代数変数 (出口面の速度)。CC は N=1 のみ対応
# (式50のCC全体の体積・質量を前提とするため)。
#
# 【入口 (液管→CC)】 液管の質量残差と同じ式 (face_upwind) で ṁ_in を計算する
# (cc_inlet)。液管側の ε_u と同じ値を渡すこと。
#
# 【出口 (CC→ウィック)】 流出するのは液だけ (ウィックはCC内の液につかって
# いる)。エンタルピーは混合物の φ ではなく、液のエンタルピー
# h_l(T,P) + g z を使う。混合物の φ を使うと、蒸気のエネルギーまで
# ウィックへ持ち出してしまう。
#
# 【CCケース (壁) からの熱】 気液別の自然対流 (式70)。
#   Ra = g β |T_wall − T| L³/(ν α_th)、Nu = Nu_cc_smooth(Ra)、h = Nu k/L
#   代表長さ L はCCの内径 (Daney の水平円筒の定義)。接触面積は
#   geometry_cc.jl の θ_liq から求める (液: A_liq、蒸気: A_vap)。
#
# 【判断事項 (ドキュメントに明記がなく、仮定を置いたもの)】
#   1. 上記の貯水池 (0次元) としての扱い。運動量式を使わない。
#   2. 自然対流の代表長さ L は、気液どちらも内径 D_i。
#      (乱流側は Nu ∝ Ra^(1/3) で h が L に依存しないが、層流側は依存する。)
#   3. CCケースの周囲への熱 (Q_amb) は、既存コードと同じく既定は 0 (断熱)。
#      ドキュメントに CC の Qamb の定義がない。
#   4. CCから流出する液のエンタルピーに含まれる位置エネルギー ṁ g z は、
#      ウィック側 (位置エネルギーを含まない) とは勘定が合わない (ṁ g z ≈ 1e-4 W
#      程度で無視できる大きさ)。
# =============================================================================

# -----------------------------------------------------------------------------
# 1. 幾何
# -----------------------------------------------------------------------------

"""
    cc_params(geo::CCGeometry) -> NamedTuple

CCの流体側の幾何量 (pipe_params に、半径 R を加えたもの)。N=1 のみ対応。
"""
function cc_params(geo::CCGeometry)
    geo.N == 1 || error("CCは N=1 (貯水池) のみ対応しています。geo.cc.N=$(geo.N)")
    return merge(pipe_params(geo), (R = geo.R,))
end

"""
    cc_wall_params(geo::CCGeometry) -> NamedTuple

CCケース (壁) の幾何量。residuals_wall.jl の wall_ghost / wall_axial_conductance
に渡せるよう、A_cross (環状断面積) と dl (軸方向長さ) を持たせている。
液管の最終セルの壁とCCケースの壁との軸方向熱伝導は、既存コードと同じく
CCケースの長さ全体 (dl=L) の半分ずつを直列に結合する形になる。
"""
cc_wall_params(geo::CCGeometry) = (
    N = geo.N, dl = geo.dl, D_i = geo.D_i, D_o = geo.D_o,
    V_w = π * (geo.D_o^2 - geo.D_i^2) / 4 * geo.dl,
    A_cross = π * (geo.D_o^2 - geo.D_i^2) / 4,
)

# -----------------------------------------------------------------------------
# 2. 入口 (液管 → CC)
# -----------------------------------------------------------------------------

"""
    cc_inlet(u_ll, p_ll, u_cc; ε_u=1e-4) -> (ṁ, φ)

液管の最終面 (= CCの入口面) を通って流入する質量流量 ṁ [kg/s] と、
その風上エンタルピー φ_in [J/kg] (全エンタルピー Em + P/ρm の風上値)。

液管の質量残差 (residual_transport! の面 N) と同じ式 face_upwind を使い、
同じ ε_u を使うことで、液管から出る質量と CC に入る質量が厳密に一致する。
u_ll, u_cc はそれぞれ液管、CC の状態 (ComponentVector)。
"""
function cc_inlet(u_ll, p_ll, u_cc; ε_u=1e-4)
    uf = u_ll.u[end]
    ρf, w = face_upwind(uf, u_ll.ρm[end], u_cc.ρm[1], ε_u)
    φ_ll = u_ll.Em[end] + u_ll.P[end] / u_ll.ρm[end]
    φ_cc = u_cc.Em[1] + u_cc.P[1] / u_cc.ρm[1]
    φ_in = w * φ_ll + (1 - w) * φ_cc
    return (ṁ = ρf * uf * p_ll.A, φ = φ_in)
end

# -----------------------------------------------------------------------------
# 3. 自然対流による CCケース → 流体の熱 (気液別)
# -----------------------------------------------------------------------------

"""
    cc_natural_convection(u_cc, Tw, p, g, tables; δ_T=1e-3) -> NamedTuple

CCケース (壁温度 Tw) から、CC内の液と蒸気への熱を、気液別の自然対流で
求める。戻り値:

    Q_l, Q_v   … 壁→液、壁→蒸気の熱 [W] (壁が熱ければ正)
    Q_fluid    … Q_l + Q_v (流体が壁から受け取る熱)
    geom       … cc_liquid_geometry の戻り値 (θ_liq, A_liq, A_vap など)
    Ra_l, Ra_v, h_l, h_v

温度 T は状態変数 u_cc.T[1] (飽和温度) を使い、物性は sat_props から取る。
ΔT=Tw−T の符号は soft_abs で処理し、Q には符号付きの ΔT を使う
(壁が冷たければ負)。
"""
function cc_natural_convection(u_cc, Tw, p, g, tables::FluidTables; δ_T=1e-3)
    T = u_cc.T[1]
    sp = sat_props(T, tables)
    f_l = liquid_volume_fraction(u_cc.ρm[1], u_cc.x[1], sp.ρ_l)
    geom = cc_liquid_geometry(f_l, p.R, p.dl)

    ΔT_abs = soft_abs(Tw - T, δ_T)
    L = p.D

    ν_l = sp.μ_l / sp.ρ_l
    α_l = sp.k_l / (sp.ρ_l * sp.cp_l)
    Ra_l = g * sp.β_l * ΔT_abs * L^3 / (ν_l * α_l)
    h_l = Nu_cc_smooth(Ra_l) * sp.k_l / L

    ν_v = sp.μ_v / sp.ρ_v
    α_v = sp.k_v / (sp.ρ_v * sp.cp_v)
    Ra_v = g * sp.β_v * ΔT_abs * L^3 / (ν_v * α_v)
    h_v = Nu_cc_smooth(Ra_v) * sp.k_v / L

    Q_l = h_l * geom.A_liq * (Tw - T)
    Q_v = h_v * geom.A_vap * (Tw - T)
    return (Q_l = Q_l, Q_v = Q_v, Q_fluid = Q_l + Q_v, geom = geom,
            Ra_l = Ra_l, Ra_v = Ra_v, h_l = h_l, h_v = h_v)
end

# -----------------------------------------------------------------------------
# 4. CC流体の残差
# -----------------------------------------------------------------------------

"""
    residual_cc!(resid, du, u, inlet, ṁ_out, Q_fluid, Q_leak, p, g, tables)

CC (貯水池、N=1) の流体残差を resid に書き込む。resid, du, u は CC の
ComponentVector (ρm,u,Em,P,T,x)。

  inlet   : cc_inlet の戻り値 (ṁ_in [kg/s], φ_in [J/kg])
  ṁ_out   : ウィックへ流出する質量流量 [kg/s] (wick_interface の ṁ)
  Q_fluid : CCケースから流体が受け取る熱 [W] (cc_natural_convection の Q_fluid)
  Q_leak  : ウィックから CC の液へ内表面を通って入る熱 [W]
            (residual_wick! の戻り値)
  p       : cc_params(geo.cc)

残差の単位は、輸送管と揃えて体積あたり (質量: kg/m³/s、エネルギー: W/m³)。
"""
function residual_cc!(resid, du, u, inlet, ṁ_out, Q_fluid, Q_leak, p, g, tables::FluidTables; ΔP_ne = 0.0)
    V = p.V
    ρ, Em, P, T = u.ρm[1], u.Em[1], u.P[1], u.T[1]
    elev = cell_elevation(p, 1)

    # --- 質量: V dρ/dt = ṁ_in − ṁ_out ---
    resid.ρm[1] = du.ρm[1] - (inlet.ṁ - ṁ_out) / V

    # --- エネルギー (保存形) ---
    # 流出するのは液だけ: エンタルピーは液の h_l(T,P) + g z
    φ_out = liquid_props(T, P, tables).h + g * elev
    Q_total = Q_fluid + Q_leak
    resid.Em[1] = ρ * du.Em[1] -
                  (inlet.ṁ * (inlet.φ - Em) - ṁ_out * (φ_out - Em) + Q_total) / V

    # --- 出口の流速 (運動量式の代わりの代数拘束) ---
    ρ_l = tables.ρ_l_sat(T)
    resid.u[1] = u.u[1] - ṁ_out / (ρ_l * p.A)

    # --- 状態拘束 (式4-6): 運動エネルギーは無視 ---
    em = Em - g * elev
    st = flash_ρe(ρ, em, tables)
    resid.P[1] = P - st.P - ΔP_ne      # ΔP_ne: 有限速度の凝縮による蒸気側の過剰圧力 (既定 0)
    resid.T[1] = T - st.T
    resid.x[1] = u.x[1] - st.x
    return nothing
end

# -----------------------------------------------------------------------------
# 5. CCケース (壁) の残差
# -----------------------------------------------------------------------------

"""
    cc_wall_heat_capacity(T, cw) -> C [J/K]

CCケースの熱容量 (SUS304、physics.jl の ρ_sus304, c_sus304)。
V_w = π(D_o²−D_i²)/4 · L (円筒側面のみ。端板は含まない)。
"""
cc_wall_heat_capacity(T, cw) = cw.V_w * ρ_sus304(T) * c_sus304(T)

"""
    Q_cc_ll_gain(T_ccw, cw, T_ll, wp_ll) -> Q [W]

液管の最終セルの壁から、CCケースが軸方向熱伝導で受け取る熱
(液管が熱ければ正)。residuals_wall.jl の wall_axial_conductance を使うので、
液管側の残差 (wall_ghost で CC を下流側ゴーストにしたもの) と対称で、
熱の出入りが厳密に一致する。
"""
Q_cc_ll_gain(T_ccw, cw, T_ll, wp_ll) = wall_axial_conductance(T_ll, wp_ll, T_ccw, cw) * (T_ll - T_ccw)

"""
    residual_cc_wall!(resid_Tw, du_Tw, Tw, Q_fluid, Q_ec_cc, Q_cc_ll, cw; Q_amb=0.0)

CCケース (壁、N=1) の残差。

    C dTw/dt = Q_ec-cc + Q_cc-ll − Q_fluid + Q_amb

  Q_fluid : CCケースから流体へ渡す熱 [W]。cc_natural_convection の Q_fluid
            (流体側の残差と同じ値を共有するので、エネルギーが厳密に一致する)
  Q_ec_cc : 蒸発器ケースから CCケースへ入る熱 [W] (式43、residuals_evapcase.jl)
  Q_cc_ll : 液管の壁から CCケースへ入る熱 [W] (Q_cc_ll_gain)
  Q_amb   : 周囲との熱。既定は 0 (断熱。判断事項3)
"""
function residual_cc_wall!(resid_Tw, du_Tw, Tw, Q_fluid, Q_ec_cc, Q_cc_ll, cw; Q_amb=0.0, C_extra=0.0)
    resid_Tw[1] = (cc_wall_heat_capacity(Tw[1], cw) + C_extra) * du_Tw[1] -
                  (Q_ec_cc + Q_cc_ll - Q_fluid + Q_amb)
    return nothing
end