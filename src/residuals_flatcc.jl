# =============================================================================
# residuals_flatcc.jl (平板型のCC: 気相 1 セル + 液相 N_CC_LAYER 層)
#
# CC は、曲面を上、平面 (primary wick 側) を下にした半円筒 (D_i = 27 mm)。気液は重力で分離し、蒸気は
# 曲面の最上部に中心角 θ_vap の円冠として溜まる (geometry_flatcc.jl)。
#
# 【旧版からの変更】 旧版は CC を気液平衡の 1 セル (貯水池) とし、液の流入による蒸気側の過剰圧力 ΔP_ne を、
# 界面の液側の 1 次元熱伝導層 (cc_layer_cond) の温度から別に与えていた。新版は次のとおり。
#   ・気相 1 セル (飽和蒸気、状態 mv) + 液相 N_CC_LAYER = 16 層 (状態 Tl[k]、層 1 = 気液界面側、層 N = 底)。
#   ・液管からの流入は液相の最下層 (層 N) に入る。ウィック (secondary wick) への流出も層 N から出る。
#   ・液の層は、液全体の体積 V_L,nom の一定の割合 f_k (等比 CC_LAYER_RATIO、界面側ほど細かい) で、液面が動くと層の境界も動く
#     (液は、層の間を上向きの質量流量 F_k で押し上げられる)。層の間は熱伝導だけ (成層の混合は無し)。
#   ・気相と層 1 の間に、界面の運動論的な凝縮 (ドキュメント 2.9 式 97〜99):
#         ṁ_cond λ = q_v,i − q_l,i = h_i (T_v − T_i),   h_i = (2σ/(2−σ)) (λ²/(T_v v_lv)) √(M/(2π R T_v)) (1 − P_v v_lv/(2λ))
#     T_v = T_sat(P) (気相の圧力 P に対する飽和温度)、T_i = 層 1 の温度、λ = h_v,sat(T_v) − h_l,sat(T_i)。
#     ΔP_ne = P − P_sat(T_1) > 0 のとき凝縮 (ṁ_cond > 0)、< 0 のとき蒸発。凝縮の質量とエンタルピーは層 1 に入る。
#     h_i は大きく (アンモニア 20 °C で約 3e7 W/m²K、界面全体で約 1.7e4 W/K)、界面の運動論の抵抗はほぼ 0 なので、
#     P ≈ P_sat(T_1) に張り付く。ΔP_ne の大きさは、層 1 の温度がどこまで上がるか (液側の熱伝導) で決まる。
#
# 【状態変数】 (types.jl の ccl / cc)
#   ccl.Tl[1:N]  液の各層の温度 [K] (微分変数)
#   ccl.X        名目の空き容積 X = V − V_L,nom [m³]。V_L,nom = 液の質量 / 飽和液密度 で見た液の体積 (微分変数)。
#                気相があるときは X ≈ 気相の体積。X < 0 は液が満水で圧縮されている状態。
#   ccl.mv       気相の質量 [kg] (微分変数)
#   cc.P         CC の圧力 (代数変数。下の閉じた式で決まる。気相の蒸気密度 = mv/V_v と、液の圧縮性)
#   cc.ρm, u, Em, T, x   旧版と同じ名前の要約 (代数変数)。診断・保存・後処理が従来どおり使えるようにしたもの。
#                T = T_sat(P) (気相の飽和温度)、x = mv/M、ρm = M/V、Em = 全内部エネルギー/M + g z、
#                u = ṁ_out/(ρ_l A) (旧版と同じ出口の流速の代数拘束)。
#
# 【圧力 P の式 (気相が消えても同じ式)】
#     mv/ρ_v,sat(P) = V_v(P) = V − V_L(P)
#     V_L(P) = V_L,nom ρb_sat / ρb(P),  ρb(P) = Σ f_k ρ_sat(T_k) (1 + (P − P_sat(T_k))/K_liq)
#   気相があるときは、左辺 (蒸気の密度) が圧力を決める (液の圧縮はほぼ効かない)。気相が尽きる (mv → 0) と、右辺の液の圧縮
#   (体積弾性率 K_liq) だけが圧力を決める (旧版が、単相液を表の flash で扱っていたのと同じ硬さ)。モードの切り替え (イベント) は無い。
#   mv が 0 に近づくときの凝縮量は、ṁ_cond に mv/(mv + m_ε) を掛けて 0 にする (mv が負にならないように)。
#
# 【液の質量と体積】 層 k の質量は m_k = f_k V_L,nom ρ_sat(T_k)。液の質量の式 dM_L/dt = ṁ_in − ṁ_out + ṁ_cond から
#   ρb dX/dt = −(ṁ_in − ṁ_out + ṁ_cond) + V_L,nom dρb/dt、層の間の質量流量 F_k は、層の質量の式 dm_k/dt = F_k − F_{k−1} から下から順に出す。
#
# 【層のエネルギー (温度の式)】 非保存形、圧力の仕事は層の質量の式の P/ρ の項と打ち消すので現れない:
#   m_k c_k dT_k/dt = F_k c_k (T*_k − T_k) − F_{k−1} c_k (T*_{k−1} − T_k) + (層間の熱伝導) + Q_k
#   T*_k = 層 k と k+1 の間の風上の温度 (流れの向きは、du を含まない近似の Ftil で決める。残差を du について線形に保つため)。
#   Q_k = secondary wick からの自然対流 (層ごと) 、層 1: + 蒸気側の壁からの熱 + ṁ_cond λ、層 N: + ṁ_in (φ_in − φ_N) + Q_cc_fluid (診断用)
#
# 【判断事項・仮定】
#   ・気相は飽和蒸気 (T_v = T_sat(P))。気相の熱容量 (約 0.05 J/K) と、断熱圧縮で蒸気が過熱する分 (約 0.03 J) は無視する。
#   ・層は、液の断面積を一様 (A_slab = V_L,nom/h_liq) として扱う (層の厚さ = 層の割合 × 液深さ h_liq)。
#     界面の面積 (凝縮の h_i A_if) だけは、半円筒の弦の幅 × 長さ (cc_vapor_geometry)。
#   ・液の圧縮性は、体積弾性率 K_liq (T_init の飽和液、表の ρ(P, e) から) 一定。
#   ・壁面凝縮 (cc_wall_cond = true): 液面の凝縮に比例して壁面 (A_vap) にも凝縮する (比 = 熱浸透率 × 面積の比 × 温度差の比 g_T (cc_wall_dT = true のとき)。フィット係数なし)。壁面の潜熱は secondary wick へ、凝縮水は層 1 へ。false なら液面だけ。
#   ・自然対流の代表長さ L は円筒型と同じ内径 D_i。平面側 (A_flat) は最下層、曲面・端面の液は層の体積の割合で配分する。
# =============================================================================

# -----------------------------------------------------------------------------
# 0. 定数・パラメータ
# -----------------------------------------------------------------------------

"""
    CC_LAYER_FRAC

液の層の体積の割合 f_k (k = 1: 気液界面側、N = 底)。等比 CC_LAYER_RATIO、合計 1。
"""
const CC_LAYER_FRAC = let r = CC_LAYER_RATIO, N = N_CC_LAYER
    [(r - 1) * r^(k - 1) / (r^N - 1) for k in 1:N]
end

const CC_R_GAS   = 8.314462618    # 一般気体定数 [J/mol K]
const CC_K_SCALE = 1.0e9          # 圧力の式 (閉じた式) を Pa の次元にそろえる係数 [Pa] (液の体積弾性率程度)
const CC_M_EPS   = 1.0e-8         # 気相の質量が尽きるときの凝縮量の抑制の尺度 [kg] (通常の気相の質量 約 2e-5 kg の 0.05 %)
const CC_ε_F     = 1.0e-8         # 層の間の流れの向きの判定の幅 [kg/s]
const CC_ε_MC    = 1.0e-8         # 凝縮・蒸発の向きの判定の幅 [kg/s]

"""
    liquid_bulk_modulus(T, tables) -> K [Pa]

飽和液 (温度 T) の体積弾性率の目安 K = ρ (∂P/∂ρ)_e。表の P(ρ, e) を、飽和液の点から ρ を 0.2 % 増やして差分する。
build_params で T_init について 1 回だけ計算する。
"""
function liquid_bulk_modulus(T, tables::FluidTables)
    ρ = tables.ρ_l_sat(T)
    e = tables.h_l_sat(T) - tables.Psat_T(T) / ρ
    ε = 2.0e-3
    P0 = flash_ρe(ρ, e, tables).P
    P1 = flash_ρe(ρ * (1 + ε), e, tables).P
    return (P1 - P0) / ε
end

"""
    pipe_params(geo::FlatCCGeometry) は residuals_transport.jl 側で共通化済み。
    flatcc_params(geo::FlatCCGeometry; K_liq, M_mol) -> NamedTuple

CC流体の幾何量 (pipe_params に、半径 R と、自然対流の面積に使う量を加えたもの)。N=1 のみ対応。
S_total は、曲面 + 両端面 + 平面側の全面積。K_liq は液の体積弾性率 [Pa]、M_mol は作動流体のモル質量 [kg/mol]
(build_params が渡す。省略すると NaN)。
"""
function flatcc_params(geo::FlatCCGeometry; K_liq = NaN, M_mol = NaN)
    geo.N == 1 || error("CCは N=1 (貯水池) のみ対応しています。geo.cc.N=$(geo.N)")
    return merge(pipe_params(geo),
                 (R = geo.R, A_flat = geo.D_i * geo.dl, S_total = geo.S_inner + geo.D_i * geo.dl,
                  K_liq = K_liq, M_mol = M_mol))
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

# -----------------------------------------------------------------------------
# 1. 底の層の状態・入口
# -----------------------------------------------------------------------------

"""
    cc_bottom_state(T_N, P, p, g, tables) -> (ρ, e, h, φ, Em)

液相の最下層 (液管からの流入・ウィックへの流出の位置) の状態。φ = h + g z (旧版の φ = Em + P/ρ と同じ位置エネルギーの勘定)、
Em = e + g z。ρ は飽和液密度 (液管の最終面の質量流束・ゴースト値と cc_inlet で同じ値を使う)。
"""
function cc_bottom_state(T_N, P, p, g, tables::FluidTables)
    lp = liquid_props(T_N, P, tables)
    z = cell_elevation(p, 1)
    return (ρ = lp.ρ, e = lp.e, h = lp.h, φ = lp.h + g * z, Em = lp.e + g * z)
end

"""
    cc_inlet(u_ll, p_ll, ρ_cc, φ_cc; ε_u=1e-4) -> (ṁ, φ)

平板型 (新しい CC): 液管の最終面 (= CC の底の層への入口面) を通って流入する質量流量 ṁ [kg/s] と、風上のエンタルピー φ_in。
旧版の cc_inlet (CC の状態 u_cc を渡す版) と同じ式 (face_upwind、同じ ε_u) で、ρ_cc, φ_cc は、液相の最下層の密度とエンタルピー (cc_bottom_state)。
"""
function cc_inlet(u_ll, p_ll, ρ_cc, φ_cc; ε_u=1e-4)
    uf = u_ll.u[end]
    ρf, w = face_upwind(uf, u_ll.ρm[end], ρ_cc, ε_u)
    φ_ll = u_ll.Em[end] + u_ll.P[end] / u_ll.ρm[end]
    φ_in = w * φ_ll + (1 - w) * φ_cc
    return (ṁ = ρf * uf * p_ll.A, φ = φ_in)
end

# -----------------------------------------------------------------------------
# 2. 幾何・圧力の閉じた式
# -----------------------------------------------------------------------------

"""
    cc_layer_geometry(X, p) -> NamedTuple

液の層の幾何。f_v = X/V から cc_vapor_geometry で液深さ h_liq を出し、層の厚さ dz_k = f_k h_liq、断面積 A_slab = V_L,nom/h_liq (一様)。
A_if は気液界面の面積 (弦の幅 × 長さ)。
"""
function cc_layer_geometry(X, p)
    geom  = cc_vapor_geometry(X / p.V, p.R, p.dl)
    h_liq = geom.h_liq
    V_nom = p.V - X
    return (geom = geom, h_liq = h_liq, V_nom = V_nom, A_slab = V_nom / h_liq,
            dz = CC_LAYER_FRAC .* h_liq,
            A_if = 2 * sqrt(max(p.R^2 - h_liq^2, 1.0e-8)) * p.dl)
end

"""
    cc_closure(P, Tl, X, p, tables) -> (V_L, V_v)

圧力 P のときの液の体積 V_L (圧縮性込み) と、残りの気相の体積 V_v = V − V_L。
"""
function cc_closure(P, Tl, X, p, tables::FluidTables)
    ρbs = 0.0
    ρb  = 0.0
    for k in 1:N_CC_LAYER
        ρs = tables.ρ_l_sat(Tl[k])
        ρbs += CC_LAYER_FRAC[k] * ρs
        ρb  += CC_LAYER_FRAC[k] * ρs * (1 + (P - tables.Psat_T(Tl[k])) / p.K_liq)
    end
    V_L = (p.V - X) * ρbs / ρb
    return (V_L = V_L, V_v = p.V - V_L)
end

"""
    cc_solve_P(Tl, X, mv, p, tables) -> P [Pa]

圧力の式 mv/ρ_v,sat(P) = V_v(P) を P について解く (二分法)。初期条件の作成と診断用。式の左辺 − 右辺は P について単調に減る。
"""
function cc_solve_P(Tl, X, mv, p, tables::FluidTables)
    Tm = sum(CC_LAYER_FRAC .* Tl)
    P0 = tables.Psat_T(Tm)
    fP(P) = mv / tables.ρ_v_sat(Tsat(P, tables)) - cc_closure(P, Tl, X, p, tables).V_v
    lo, hi = 0.8 * P0, 1.2 * P0
    n = 0
    while fP(hi) > 0 && n < 60
        hi *= 1.3; n += 1
    end
    n = 0
    while fP(lo) < 0 && n < 20
        lo *= 0.85; n += 1
    end
    (fP(lo) >= 0 && fP(hi) <= 0) || error("CC の圧力の式の根が [", lo, ", ", hi, "] Pa に見つかりません")
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        fP(mid) > 0 ? (lo = mid) : (hi = mid)
        (hi - lo) < 1.0e-9 * hi && break
    end
    return 0.5 * (lo + hi)
end

# -----------------------------------------------------------------------------
# 3. 気液界面の凝縮 (ドキュメント 2.9 式 97〜99)
# -----------------------------------------------------------------------------

"""
    cc_condensation(P, T1, X, mv, p, tables) -> NamedTuple

気相 (圧力 P、T_v = T_sat(P)) と層 1 (温度 T1) の間の凝縮の質量流量 ṁ [kg/s] (凝縮が正)。
  ṁ_raw = h_i A_if (T_v − T1)/λ ,   λ = h_v,sat(T_v) − h_l,sat(T1)
  h_i = (2σ/(2−σ)) (λ²/(T_v v_lv)) √(M/(2π R T_v)) (1 − P v_lv/(2λ)),   σ = [1 − (v_l/v_v)^(1/3)] exp(−1/(2 (v_v/v_l)^(1/3) − 1))
(σ の式 (99) の指数の分母は、σ が 0〜1 になる 2 (v_v/v_l)^(1/3) − 1 を採用。アンモニア 20 °C で σ ≈ 0.69。)
気相が尽きるとき (mv → 0) は、凝縮を mv/(mv + m_ε) で抑える (蒸発は抑えない)。
"""
function cc_condensation(P, T1, X, mv, p, tables::FluidTables; swp = nothing, T_sw = NaN, Tbar = NaN)
    T_v = Tsat(P, tables)
    sp  = sat_props(T_v, tables)
    v_l = 1 / sp.ρ_l
    v_v = 1 / sp.ρ_v
    v_lv = v_v - v_l
    λ   = sp.h_v - tables.h_l_sat(T1)
    r3  = (v_v / v_l)^(1 / 3)
    σa  = (1 - 1 / r3) * exp(-1 / (2 * r3 - 1))
    h_i = 2 * σa / (2 - σa) * λ^2 / (T_v * v_lv) * sqrt(p.M_mol / (2π * CC_R_GAS * T_v)) * (1 - P * v_lv / (2λ))
    lg = cc_layer_geometry(X, p)
    A_if = lg.A_if
    ṁ_raw = h_i * A_if * (T_v - T1) / λ
    mvp = 0.5 * (mv + sqrt(mv^2 + 1.0e-14^2))
    s_v = mvp / (mvp + CC_M_EPS)
    w   = smooth_step(ṁ_raw, CC_ε_MC)
    # 壁面凝縮 (swp = secondary wick のパラメータ。nothing なら無し): 蒸気が接する壁 (A_vap) にも凝縮する。
    # 液面と壁面はどちらも半無限体として蒸気に接するので、凝縮量の比は熱浸透率 e = √(kρc) × 面積 の比 (フィット係数なし)。
    # 壁面の凝縮は液面の凝縮 (ṁ_raw > 0 の分) に比例させ、蒸発には付けない。凝縮水は液面 (層 1) へ流れ落ちる。
    # T_sw, Tbar が与えられたとき (NaN でないとき) は、駆動する温度差の比 g_T = (T_v − T_sw)/(T_v − Tbar) をかける
    # (T_sw = secondary wick の温度、Tbar = 液の体積平均温度 = それぞれの遠方の温度)。半無限体の熱流束 e (T_v − T_遠方)/√(πt) の、温度差の部分。
    # 分母の下限 0.05 K (起動前で T_v ≈ Tbar のとき。このとき ṁ_raw ≈ 0 なので影響しない)、上限 g_max = 2 は数値の安全のためで、通常は効かない。
    g_T = 1.0
    r_w = 0.0
    ṁ_w = 0.0
    if swp !== nothing
        e_l = sqrt(sp.k_l * sp.ρ_l * sp.cp_l)
        e_w = sqrt(swp.k_eff * ((1 - swp.por) * swp.ρ_w * swp.C_w + swp.por * sp.ρ_l * sp.cp_l))
        r_w = e_w * lg.geom.A_vap / (e_l * A_if)
        if !isnan(T_sw) && !isnan(Tbar)
            den = soft_clamp(T_v - Tbar, 0.05, 10.0, 0.01)
            g_T = soft_clamp((T_v - T_sw) / den, 0.0, 2.0, 0.02)
        end
        ṁ_w = r_w * g_T * ṁ_raw * w * s_v
    end
    return (ṁ = ṁ_raw * (w * s_v + (1 - w)), ṁ_raw = ṁ_raw, h_i = h_i, λ = λ, T_v = T_v, A_if = A_if, σ = σa,
            ṁ_w = ṁ_w, r_w = r_w, g_T = g_T, Q_w = ṁ_w * λ)
end

# -----------------------------------------------------------------------------
# 4. secondary wick → CC流体の自然対流 (層ごと)
# -----------------------------------------------------------------------------

"""
    flatcc_natural_convection(ccl, P, T_A, T_B, p, g, tables; δ_T=1e-3, boiling=false, r_c=5e-6, ε_onset=0.1)

secondary wick の 2 節点 (A: 平面側 + 両端面、温度 T_A。B: 半円環 = 曲面、温度 T_B) から、CC内の液 (各層) と蒸気への熱を、
気液別の自然対流で求める (式70)。

    Ql        … 各層が受け取る熱 [W] (長さ N_CC_LAYER。secondary wick が熱ければ正)
    Q_l, Q_v  … 液全体、蒸気が受け取る熱 [W]
    Q_fluid   … Q_l + Q_v (流体が受け取る熱)
    Q_A, Q_B  … 節点 A、B が流体へ渡す熱 (swick の残差で、それぞれの節点から引く。Q_A + Q_B = Q_fluid)
    geom, Ra_l, Ra_v, h_l (液の面積平均), h_v, w_on

面積: 蒸気 = 曲面の R θ_vap dl (B) + 両端面の円冠 R²(θ_vap − sin θ_vap) (A)。
      液   = 曲面の濡れ R (π − θ_vap) dl (B) + 両端面の濡れ π R² − 円冠 (A) + 平面側 A_flat (A)。
液の面積は、平面側 A_flat を最下層に、残りを層の体積の割合 f_k で各層に配る。各層の ΔT = T_節点 − T_k、物性は層の平均温度の飽和液。
"""
function flatcc_natural_convection(ccl, P, T_A, T_B, p, g, tables::FluidTables; δ_T=1e-3,
    boiling=false, r_c=5.0e-6, ε_onset=0.1)
    N = N_CC_LAYER
    f = CC_LAYER_FRAC
    Tl = ccl.Tl
    T_v = Tsat(P, tables)
    geom = cc_vapor_geometry(ccl.X[1] / p.V, p.R, p.dl)
    θ = geom.θ_vap
    A_vap_c = p.R * θ * p.dl                       # 曲面の蒸気側 (B)
    A_vap_e = p.R^2 * (θ - sin(θ))                 # 両端面の円冠 (A)
    A_vap = A_vap_c + A_vap_e
    A_liq = p.S_total - A_vap
    A_cw = p.R * (π - θ) * p.dl                    # 曲面の濡れ (B)
    A_ew = A_liq - p.A_flat - A_cw                 # 両端面の濡れ (A)
    L = p.D

    Tm = sum(f .* Tl)
    sp = sat_props(Tm, tables)
    ν_l = sp.μ_l / sp.ρ_l
    α_l = sp.k_l / (sp.ρ_l * sp.cp_l)
    sv = sat_props(T_v, tables)

    # 核沸騰 (ドキュメント 2.8.2): CC での沸騰は自然対流が支配的なので、プール沸騰の h を使う。
    # 節点の表面温度が T_sat(P_cc) + ΔT_onset を超えたら沸騰が起きる (式82、r_c は secondary wick の孔径)。
    ΔT_onset = 2 * sv.σ * T_v / (sv.ρ_v * sv.h_fg * r_c)
    function boil(Tw)
        boiling || return (0.0, 0.0)
        hc = h_chen_smooth(Tw, T_v, Psat(Tw, tables), P, tables.Pcrit, 0.0, 1.0e-3, p.D, sv, sv.h_fg)
        return (hc, smooth_step(Tw - T_v - ΔT_onset, ε_onset))
    end
    hcA, wA = boil(T_A)
    hcB, wB = boil(T_B)

    Ql = zeros(N)
    Q_Al = 0.0
    Q_Bl = 0.0
    hA  = 0.0
    Ra_sum = 0.0
    for k in 1:N
        ΔT_A = T_A - Tl[k]
        ΔT_B = T_B - Tl[k]
        Ra_A = g * sp.β_l * soft_abs(ΔT_A, δ_T) * L^3 / (ν_l * α_l)
        Ra_B = g * sp.β_l * soft_abs(ΔT_B, δ_T) * L^3 / (ν_l * α_l)
        h_A = (1 - wA) * Nu_cc_smooth(Ra_A) * sp.k_l / L + wA * hcA
        h_B = (1 - wB) * Nu_cc_smooth(Ra_B) * sp.k_l / L + wB * hcB
        A_Ak = A_ew * f[k] + (k == N ? p.A_flat : 0.0)
        A_Bk = A_cw * f[k]
        QA_k = h_A * A_Ak * ΔT_A
        QB_k = h_B * A_Bk * ΔT_B
        Ql[k] = QA_k + QB_k
        Q_Al += QA_k
        Q_Bl += QB_k
        hA += h_A * A_Ak + h_B * A_Bk
        Ra_sum += f[k] * 0.5 * (Ra_A + Ra_B)
    end

    ν_v = sv.μ_v / sv.ρ_v
    α_v = sv.k_v / (sv.ρ_v * sv.cp_v)
    function h_vap(Tw)
        Ra = g * sv.β_v * soft_abs(Tw - T_v, δ_T) * L^3 / (ν_v * α_v)
        return Nu_cc_smooth(Ra) * sv.k_v / L, Ra
    end
    h_vA, Ra_vA = h_vap(T_A)
    h_vB, Ra_vB = h_vap(T_B)
    Q_vA = h_vA * A_vap_e * (T_A - T_v)
    Q_vB = h_vB * A_vap_c * (T_B - T_v)
    Q_l = Q_Al + Q_Bl
    Q_v = Q_vA + Q_vB

    return (Ql = Ql, Q_l = Q_l, Q_v = Q_v, Q_fluid = Q_l + Q_v, Q_A = Q_Al + Q_vA, Q_B = Q_Bl + Q_vB, geom = geom,
            Ra_l = Ra_sum, Ra_v = Ra_vB, h_l = hA / A_liq, h_v = h_vB, w_on = 0.5 * (wA + wB))
end

# -----------------------------------------------------------------------------
# 5. 液相の残差 (層の温度 Tl、空き容積 X、気相の質量 mv)
# -----------------------------------------------------------------------------

"""
    residual_ccl!(r, d, x, P, inlet, ṁ_out, nc, Q_extra, p, g, tables) -> NamedTuple

CC の液相 (各層の温度 Tl、名目の空き容積 X) と気相の質量 mv の残差を r に書き込む。r, d, x は ccl の ComponentVector (Tl, X, mv)。

  P       : CC の圧力 (cc.P、代数変数)
  inlet   : cc_inlet (液管の最終面の ṁ_in [kg/s] と風上のエンタルピー φ_in)。液相の最下層に入る
  ṁ_out   : 最下層からウィック (secondary wick) へ出る質量流量 [kg/s] (wick_interface の ṁ)
  nc      : flatcc_natural_convection の戻り値 (Ql: 各層、Q_v: 蒸気側の壁から。蒸気側は層 1 に入れる)
  Q_extra : 最下層へ入る診断用の熱 [W] (Q_cc_fluid)

戻り値: ṁ_cond (凝縮が正)、cd (cc_condensation の戻り値)、F (層の間の質量流量 [kg/s]。F[k] = 層 k+1 → 層 k の上向き、F[N] = ṁ_in − ṁ_out)。
残差は du について線形 (流れの向きの判定には du を含まない近似の Ftil を使う)。
"""
function residual_ccl!(r, d, x, P, inlet, ṁ_out, nc, Q_extra, p, g, tables::FluidTables; swp = nothing, T_sw = NaN)
    N = N_CC_LAYER
    f = CC_LAYER_FRAC
    Tl = x.Tl
    X  = x.X[1]
    mv = x.mv[1]
    V_nom = p.V - X

    ρs  = [tables.ρ_l_sat(Tl[k]) for k in 1:N]
    cps = [tables.cp_l_sat(Tl[k]) for k in 1:N]
    δT  = 0.05
    dρdT = [(tables.ρ_l_sat(Tl[k] + δT) - tables.ρ_l_sat(Tl[k] - δT)) / (2δT) for k in 1:N]
    ρb = sum(f .* ρs)
    dρbdt = sum(f .* dρdT .* d.Tl)

    # --- 界面の凝縮と、液・気相の質量 ---
    cd = cc_condensation(P, Tl[1], X, mv, p, tables; swp = swp, T_sw = T_sw, Tbar = sum(f .* Tl))
    ṁ_cond = cd.ṁ                        # 液面での凝縮 (潜熱は層 1 へ)
    ṁ_tot  = cd.ṁ + cd.ṁ_w               # 気相の質量の減少 = 液面 + 壁面 (壁面の凝縮水は液面 (層 1) へ流れ落ちる。潜熱は壁 = secondary wick へ)
    Δṁ = inlet.ṁ - ṁ_out                 # 最下層への正味の流入
    dMLdt = Δṁ + ṁ_tot                   # 液の質量の変化率
    r.X[1]  = ρb * d.X[1] + (dMLdt - V_nom * dρbdt)       # [kg/s]
    r.mv[1] = d.mv[1] + ṁ_tot                           # [kg/s]

    # --- 層の間の質量流量 F (上向き正)。F[k]: 層 k+1 → 層 k (k = 1..N−1)、F[N]: 底の境界の正味の流入 ---
    dmdt = [f[k] * (-ρs[k] * d.X[1] + V_nom * dρdT[k] * d.Tl[k]) for k in 1:N]
    F = zeros(N)
    F[N] = Δṁ
    for k in N:-1:2
        F[k-1] = F[k] - dmdt[k]
    end
    # 流れの向きの判定用 (du を含まない近似): Ftil[k] = Δṁ − (Δṁ + ṁ_cond) Σ_{j>k} f_j
    Tstar = zeros(N - 1)
    below = 0.0
    for k in N-1:-1:1
        below += f[k+1]
        Ftil = Δṁ - dMLdt * below
        w = smooth_step(Ftil, CC_ε_F)
        Tstar[k] = w * Tl[k+1] + (1 - w) * Tl[k]
    end

    # --- 層の間の熱伝導 (層 k → k+1) ---
    geo = cc_layer_geometry(X, p)
    Qc = zeros(N - 1)
    for k in 1:N-1
        kl = tables.k_l_sat(0.5 * (Tl[k] + Tl[k+1]))
        Qc[k] = kl * geo.A_slab / (0.5 * (geo.dz[k] + geo.dz[k+1])) * (Tl[k] - Tl[k+1])
    end

    # --- 層の温度の式 ---
    bt = cc_bottom_state(Tl[N], P, p, g, tables)
    for k in 1:N
        adv = 0.0
        k <= N - 1 && (adv += F[k] * cps[k] * (Tstar[k] - Tl[k]))        # 下の面 (層 k+1 との間)
        k >= 2     && (adv -= F[k-1] * cps[k] * (Tstar[k-1] - Tl[k]))    # 上の面 (層 k−1 との間)
        cnd = (k >= 2 ? Qc[k-1] : 0.0) - (k <= N - 1 ? Qc[k] : 0.0)
        Q = nc.Ql[k]
        k == 1 && (Q += nc.Q_v + ṁ_cond * cd.λ)                         # 蒸気側の壁から、凝縮の潜熱 (h_v − h_l,sat(T1))
        k == N && (Q += inlet.ṁ * (inlet.φ - bt.φ) + Q_extra)           # 液管からの流入 (φ_in − φ_N)、診断用の熱
        r.Tl[k] = f[k] * V_nom * ρs[k] * cps[k] * d.Tl[k] - (adv + cnd + Q)
    end
    return (ṁ_cond = ṁ_cond, cd = cd, F = F, Q_wall = cd.Q_w)
end

"""
    residual_cc_summary!(r, x_cc, x_ccl, ṁ_out, p, g, tables)

CC の圧力 cc.P (閉じた式) と、旧版と同じ名前の要約 cc.ρm, u, Em, T, x の残差 (全て代数拘束)。r, x_cc は cc の ComponentVector。
  P : (CC_K_SCALE/V) (mv/ρ_v,sat(P) − V_v(P)) = 0
  T = T_sat(P),  x = mv/M,  ρm = M/V,  Em = (mv u_v + Σ m_k e_k)/M + g z,  u = ṁ_out/(ρ_l,sat(T) A)
"""
function residual_cc_summary!(r, x_cc, x_ccl, ṁ_out, p, g, tables::FluidTables)
    P  = x_cc.P[1]
    Tl = x_ccl.Tl
    X  = x_ccl.X[1]
    mv = x_ccl.mv[1]
    T_v = Tsat(P, tables)
    ρ_v = tables.ρ_v_sat(T_v)
    cl  = cc_closure(P, Tl, X, p, tables)
    r.P[1] = (CC_K_SCALE / p.V) * (mv / ρ_v - cl.V_v)

    sm = cc_summary_values(P, Tl, X, mv, p, g, tables)
    r.ρm[1] = x_cc.ρm[1] - sm.ρm
    r.x[1]  = x_cc.x[1] - sm.x
    r.T[1]  = x_cc.T[1] - sm.T
    r.Em[1] = x_cc.Em[1] - sm.Em
    r.u[1]  = x_cc.u[1] - ṁ_out / (tables.ρ_l_sat(x_cc.T[1]) * p.A)
    return nothing
end

"""
    cc_summary_values(P, Tl, X, mv, p, g, tables) -> (ρm, x, T, Em)

CC の要約 (旧版と同じ名前の変数 cc.ρm, x, T, Em の値)。残差と初期条件で同じ式を使う。
  T = T_sat(P),  M = mv + Σ m_k (m_k = f_k V_L,nom ρ_sat(T_k)),  ρm = M/V,  x = mv/M,  Em = (mv u_v + Σ m_k e_k)/M + g z
"""
function cc_summary_values(P, Tl, X, mv, p, g, tables::FluidTables)
    f = CC_LAYER_FRAC
    V_nom = p.V - X
    T_v = Tsat(P, tables)
    ρ_v = tables.ρ_v_sat(T_v)
    M_L = 0.0
    U_L = 0.0
    for k in 1:N_CC_LAYER
        ρk = tables.ρ_l_sat(Tl[k])
        mk = f[k] * V_nom * ρk
        M_L += mk
        U_L += mk * (tables.h_l_sat(Tl[k]) - tables.Psat_T(Tl[k]) / ρk)
    end
    u_v = tables.h_v_sat(T_v) - P / ρ_v
    M = mv + M_L
    return (ρm = M / p.V, x = mv / M, T = T_v, Em = (mv * u_v + U_L) / M + g * cell_elevation(p, 1))
end

"""
    cc_state_diag(c, p, tables) -> NamedTuple

components_at の結果 c から、新しい CC の診断量: P, T_v (気相の飽和温度)、T1 (層 1 = 界面側)、TN (層 N = 底)、f_v (気相の体積分率)、
mv [kg]、dP_ne = P − P_sat(T1) [Pa]、ṁ_cond [kg/s] (凝縮が正)、h_i [W/m²K]、A_if [m²]、h_liq [m]。
"""
function cc_state_diag(c, p, tables::FluidTables)
    pc = p.p_cc
    P  = c.cc.P[1]
    Tl = c.ccl.Tl
    X  = c.ccl.X[1]
    mv = c.ccl.mv[1]
    cd = cc_condensation(P, Tl[1], X, mv, pc, tables; swp = (p.cond.cc_wall_cond ? p.swick_p : nothing),
                         T_sw = (p.cond.cc_wall_dT ? c.swick.T[2] : NaN), Tbar = sum(CC_LAYER_FRAC .* Tl))
    geo = cc_layer_geometry(X, pc)
    cl = cc_closure(P, Tl, X, pc, tables)
    return (P = P, T_v = cd.T_v, T1 = Tl[1], TN = Tl[end], f_v = cl.V_v / pc.V, mv = mv,
            dP_ne = P - tables.Psat_T(Tl[1]), ṁ_cond = cd.ṁ, h_i = cd.h_i, A_if = cd.A_if, h_liq = geo.h_liq,
            λ = cd.λ, σ = cd.σ, ṁ_w = cd.ṁ_w, Q_wall = cd.Q_w, r_w = cd.r_w, g_T = cd.g_T)
end