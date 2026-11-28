# =============================================================================
# heattransfer.jl (Level 3: 領域ブレンド済み相関 / Level 4: コンポーネント別
# ディスパッチ)
#
# correlations.jl (Level1) の個別相関式を、smoothing.jl (Level0) の
# tanh/べき乗ブレンドで統合し (Level3)、最終的に残差関数から直接呼べる
#「ある壁面・ある流体セルの熱伝達係数/摩擦力」を返す関数 (Level4) を
# 提供する。
#
# 【要確認事項 (ドキュメントとの照合で見つかったもの)】
#   1. 式57: 原ドキュメントには hconv = hI (Regime I); hI+hLT (Regime II)
#      と書かれているが、原典 (Shah 2009/2013) の Regime II は h_I + h_Nu
#      (式56のNusselt型重力項を加算)。「hLT」は「hNu」の誤記。
#      【確認済み (2026-11-02)】本実装は h_I + h_Nu で、原典どおり。
#   2. 式62 S(Retp) の Re_tp がドキュメント中で定義されていない。
#      Chenの式の標準的な定義 (Re_tp = Re_l * F^1.25) を採用している。
#      Carey原典 (参考文献[11], 式61の修正でも参照したもの) での
#      確認を推奨。
#   3. 【確認済み】式64 Dougall-RohsenowのRev。原典 (Dougall & Rohsenow,
#      MIT report 9079-26, 式26) で確認した結果、括弧内は
#      throughput Reynolds数 Re_t = (ρ_v D/μ_v)(Q_l+Q_v)/A_P であり、
#      Re_v = G*D/μ_v (全質量流束ベース) に [x+(1-x)ρ_v/ρ_l] を掛けた
#      形と等価。本実装 (Re_v = G*D/μ_v) は原典と一致する。
#      物性はすべて飽和温度で評価する (sat_props を使う本実装と一致)。
#      なお原典の Re_f (式26a) は、乾き度を壁面熱流束と距離で表した
#      別表現で、局所の乾き度が既知の本モデルでは不要。
# =============================================================================

# =============================================================================
# Level 3: 領域ブレンド済み相関
# =============================================================================

# -----------------------------------------------------------------------------
# 3.1 無次元数の安全な計算
# -----------------------------------------------------------------------------

"""
    quality_for_correlation(x; δ=1e-4)

二相の相関式に渡す乾き度を、定義できる範囲 [δ, 1−δ] に、滑らかに収める。

Shah, Chen, Dougall-Rohsenow の式や、Lockhart-Martinelli の摩擦は、乾き度 x に
ついて (1−x)^0.8, (1/x−1)^0.8, ((1−x)/x)^0.9 などの小数べきを含み、
x≤0 または x≥1 では、底が 0 または負になって定義できない (Julia では
DomainError)。flash の乾き度 x は、単相では 0 または 1 ちょうどで、
補間や摂動で、わずかに [0,1] の外に出ることもある。

soft_floor (下限だけを守る) では、不十分だった: soft_floor(x, δ) は、常に x より
わずかに大きい値 (≈ x + δ²/(x−δ)) を返すため、x=1 のとき 1.00000001 となり、
1/x − 1 = −1e-8 < 0 で、DomainError になる。soft_clamp は、両側を守り、
x∈[δ, 1−δ] の範囲では、ほぼ恒等 (x=0.3,0.5,0.9 での差は 1e-8 程度)。
"""
quality_for_correlation(x; δ=1.0e-4) = soft_clamp(x, δ, 1 - δ, δ)

"""
    reynolds_number_safe(ρ, u, D, μ; δ_u=1e-4)

Re = ρ|u|D/μ を、流れ方向反転時にも滑らかになるよう soft_abs を
使って計算する。
"""
reynolds_number_safe(ρ, u, D, μ; δ_u=1e-4) = ρ * soft_abs(u, δ_u) * D / μ

# -----------------------------------------------------------------------------
# 3.2 単相強制対流・摩擦係数 (式27, 式50)
# -----------------------------------------------------------------------------

"""
    dittus_boelter_exponent_smooth(T_wall, T_bulk; ε=0.5)

Dittus-Boelter式のPr指数 (加熱時0.4 / 冷却時0.3, 式50の説明文) を、
壁と流体の温度差の符号でなめらかに切り替える。
"""
dittus_boelter_exponent_smooth(T_wall, T_bulk; ε=0.5) = blend_tanh(0.4, 0.3, T_wall - T_bulk, ε)

"""
    Nu_pipe_smooth(Re, Pr, n; ε_Re=200.0)

単相強制対流のNu (式50) を、層流⇔乱流でtanhブレンドする。
"""
function Nu_pipe_smooth(Re, Pr, n; ε_Re=200.0)
    Nu_lam = Nu_pipe_laminar()
    Nu_turb = Nu_pipe_turbulent(max(Re, 1.0), Pr, n)
    return blend_tanh(Nu_turb, Nu_lam, Re - 2300.0, ε_Re)
end

"""
    friction_factor_smooth(Re; ε_Re=200.0)

摩擦係数 (式27) を層流⇔乱流でtanhブレンドする。
"""
function friction_factor_smooth(Re; ε_Re=200.0)
    f_lam = friction_factor_laminar(max(Re, 1.0))
    f_turb = friction_factor_turbulent(max(Re, 1.0))
    return blend_tanh(f_turb, f_lam, Re - 2300.0, ε_Re)
end

"""
    C_lockhart_martinelli_smooth(Re_gas, Re_liq; ε_Re=200.0)

Lockhart-Martinelli型の二相圧力損失係数C。気相・液相それぞれの
Reynolds数に対する2つの独立なtanh重みで、C_LM_table (correlations.jl)
の4値をバイリニア的にブレンドする (以前検算済み)。
"""
function C_lockhart_martinelli_smooth(Re_gas, Re_liq; ε_Re=200.0)
    w_gas = smooth_step(Re_gas - 2300.0, ε_Re)
    w_liq = smooth_step(Re_liq - 2300.0, ε_Re)
    return w_gas * w_liq * C_LM_table.gas_turb_liq_turb +
           w_gas * (1 - w_liq) * C_LM_table.gas_turb_liq_lam +
           (1 - w_gas) * w_liq * C_LM_table.gas_lam_liq_turb +
           (1 - w_gas) * (1 - w_liq) * C_LM_table.gas_lam_liq_lam
end

# -----------------------------------------------------------------------------
# 3.3 凝縮 (Shah) (式51-57)
# -----------------------------------------------------------------------------

"""
    h_condense_shah_smooth(x, G, D, g, P, Pcrit, props; ε_J=0.1)

Shahの凝縮熱伝達係数 (式51-57) を、Regime I/IIでtanhブレンドする。
props は sat_props(T, tables) の戻り値 (μ_l,μ_v,ρ_l,ρ_v,k_l,Pr_l 等
を含むNamedTuple)。

【確認済み (2026-11-02)】式57のRegime IIは h_I + h_Nu で、原典どおり
(ファイル冒頭のコメント参照)。
"""
function h_condense_shah_smooth(x, G, D, g, P, Pcrit, props; ε_J=0.1)
    P_r = soft_floor(P / Pcrit, 1.0e-4)   # P_r^0.4 などの底。表の外挿で P≤0 になっても定義域エラーにしない
    x_safe = quality_for_correlation(x)

    ρ_l, ρ_v = props.ρ_l, props.ρ_v
    μ_l, μ_v = props.μ_l, props.μ_v
    k_l, Pr_l = props.k_l, props.Pr_l

    Jv = Jv_shah(x_safe, G, g, D, ρ_v, ρ_l)
    Z = Z_shah(x_safe, P_r)
    Jv_b = Jv_boundary_shah(Z)

    Re_LT = G * D / μ_l
    Re_LS = G * (1 - x_safe) * D / μ_l

    h_LT = h_LT_shah(k_l, D, Re_LT, Pr_l)
    n = n_shah(P_r)
    h_I = h_I_shah(h_LT, μ_l, μ_v, n, x_safe, P_r)
    h_Nu = h_Nu_shah(Re_LS, ρ_l, ρ_v, g, k_l, μ_l)

    h_regimeII = h_I + h_Nu

    return blend_tanh(h_I, h_regimeII, Jv - Jv_b, ε_J)
end

# -----------------------------------------------------------------------------
# 3.4 沸騰 (Chen / Dougall-Rohsenow) (式58-64)
# -----------------------------------------------------------------------------

"""
    h_chen_smooth(T_wall, T_sat, P_sat_Twall, P, Pcrit, x, G, D, props, h_fg; ε_invXtt=0.02)

Chenの沸騰熱伝達係数 (式58-63) をまとめて計算する。
F(Xtt)のRegime境界 (1/Xtt=0.1) をtanhブレンドする
(式61は以前の会話で原典に基づき修正済みの形を使用)。

【要確認2】S(Re_tp)のRe_tpは Re_tp = Re_l * F^1.25 (Chenの標準的な
定義) として実装している。
"""
function h_chen_smooth(T_wall, T_sat, P_sat_Twall, P, Pcrit, x, G, D, props, h_fg; ε_invXtt=0.02)
    P_r = soft_floor(P / Pcrit, 1.0e-4)   # P_r^0.4 などの底。表の外挿で P≤0 になっても定義域エラーにしない
    x_safe = quality_for_correlation(x)

    ρ_l, ρ_v = props.ρ_l, props.ρ_v
    μ_l = props.μ_l
    k_l, cp_l, Pr_l, σ = props.k_l, props.cp_l, props.Pr_l, props.σ

    Xtt = Xtt_martinelli(x_safe, ρ_v, ρ_l, μ_l, props.μ_v)
    inv_Xtt = 1 / Xtt

    F = blend_tanh(F_chen_highX(inv_Xtt), F_chen_lowX(), inv_Xtt - 0.1, ε_invXtt)

    Re_l = G * (1 - x_safe) * D / μ_l
    h_mac = h_mac_chen(k_l, D, Re_l, Pr_l, F)

    ΔT_sup = soft_floor(T_wall - T_sat, 0.0)
    ΔP_sat = soft_floor(P_sat_Twall - P, 0.0)
    Re_tp = Re_l * F^1.25
    S = S_chen(Re_tp)

    h_mic = h_mic_chen(k_l, cp_l, ρ_l, σ, μ_l, h_fg, ρ_v, ΔT_sup, ΔP_sat, S)

    return h_mac + h_mic
end

"""
    h_dougall_rohsenow_smooth(x, G, D, props)

ドライアウト後の熱伝達係数 (式64)。

Re_v = G*D/μ_v (全質量流束ベース)。原典 (Dougall & Rohsenow, 式26) の
throughput Reynolds数と一致することを確認済み。
"""
function h_dougall_rohsenow_smooth(x, G, D, props)
    ρ_l, ρ_v, μ_v, k_v, Pr_v = props.ρ_l, props.ρ_v, props.μ_v, props.k_v, props.Pr_v
    Re_v = G * D / μ_v
    x_safe = quality_for_correlation(x)
    return h_pd_dougall_rohsenow(k_v, D, Re_v, x_safe, ρ_v, ρ_l, Pr_v)
end

"""
    h_boiling_smooth(x, h_chen, h_pd; ε_x=0.05, x_center=0.8)

沸騰熱伝達係数を、Chenの式 (0<x<0.7) とDougall-Rohsenow (0.9<x<1.0)
でtanhブレンドする (0.7<x<0.9の遷移域をカバー)。原文の3次スプライン
の代わりに採用 (以前の議論: オーバーシュート回避、Jacobianとの相性)。
"""
h_boiling_smooth(x, h_chen, h_pd; ε_x=0.05, x_center=0.8) = blend_tanh(h_pd, h_chen, x - x_center, ε_x)

# -----------------------------------------------------------------------------
# 3.5 CC自然対流 (式65)
# -----------------------------------------------------------------------------

"""
    Nu_cc_smooth(Ra; ε_cond=..., ε_turb=..., Ra_c1=1e4, Ra_c2=None)

CC内部の自然対流Nu (式65) を、伝導→層流→乱流の3段階でべき乗
ブレンドする。Ra→0でNu→1 (伝導限界) に正しく収束する設計
(以前の議論: 層流Nu単体ではRa→0でNu→0になってしまう問題への対策)。

層流・乱流の交差点Ra_c2は、係数比から解析的に求まる
((0.54/0.15)^12 ≈ 4.7e6 程度) が、ここでは引数として明示的に
渡せるようにしている。
"""
function Nu_cc_smooth(Ra; n1=6, n2=6)
    Ra_safe = soft_floor(Ra, 1.0)   # Ra<0 (ΔT<0側) は呼び出し側でsoft_absにより
                                      # 事前に非負化されている想定。念のための保険。
    Nu_cond = Nu_cc_conduction()
    Nu_lam = Nu_cc_laminar(Ra_safe)
    Nu_turb = Nu_cc_turbulent(Ra_safe)

    step1 = blend_power(Nu_cond, Nu_lam, n1)
    return blend_power(step1, Nu_turb, n2)
end

"""
    Nu_plate_down_smooth(Ra; n=6)

加熱面が下向きの水平平板の自然対流Nu (Nu_plate_down = 0.52 Ra^(1/5)) に、伝導の下限 Nu = 1
(Nu_cc_conduction) をべき乗ブレンドしたもの。Nu_cc_smooth と同じく、Ra→0 で Nu→0 になって熱の経路が
消えるのを防ぐための数値的な下限 (2 つの式の交点は Ra = (1/0.52)^5 ≈ 26 で、相関式の適用範囲
Ra ≳ 10^4 ではブレンドの影響は無視できる)。
"""
function Nu_plate_down_smooth(Ra; n=6)
    Ra_safe = soft_floor(Ra, 1.0)
    return blend_power(Nu_cc_conduction(), Nu_plate_down(Ra_safe), n)
end

# -----------------------------------------------------------------------------
# 3.6 グルーブ過熱・周囲環境 (式66-67)
# -----------------------------------------------------------------------------

"""
    Nu_flatplate_smooth(Re, Pr; ε_Re=2e4)

グルーブ蒸気への過熱、平板強制対流のNu (式67) を層流⇔乱流で
tanhブレンドする (遷移基準Re=5×10^5)。
"""
function Nu_flatplate_smooth(Re, Pr; ε_Re=2.0e4)
    Nu_lam = Nu_flatplate_laminar(max(Re, 1.0), Pr)
    Nu_turb = Nu_flatplate_turbulent(max(Re, 1.0), Pr)
    return blend_tanh(Nu_turb, Nu_lam, Re - 5.0e5, ε_Re)
end

"""
    h_ambient_smooth(T_wall, T_amb, D_out; δ_T=0.01)

周囲環境への自然熱伝達係数 (式66, Holman)。ΔTの符号に関わらず
なめらかに計算できるよう soft_abs を適用する。
"""
function h_ambient_smooth(T_wall, T_amb, D_out; δ_T=0.01)
    ΔT_abs = soft_abs(T_wall - T_amb, δ_T)
    return h_nc_holman(ΔT_abs, D_out)
end

# =============================================================================
# Level 4: コンポーネント別ディスパッチ
# =============================================================================

"""
    h_transport_wall(ρ, u, Em, elev, T_wall, D, G, g, tables::FluidTables; kwargs...)

輸送管・グルーブ・CC (いずれも6変数のTransportState) の、流体-壁間
熱伝達係数を1つの値として返す。単相 (液/気)⇔二相 (凝縮⇔沸騰) の
すべての切り替えを、疑似乾き度・tanhブレンドで滑らかに統合する。

局所の(ρ,e)における輸送物性 (μ,k,cp) は transport_props から直接
取得するため、「単相液側/単相気側」で別々の計算をする必要がなく
(その場の(ρ,e)が自動的に正しい相の物性を与えるため)、疑似乾き度に
よる単相⇔二相のブレンドだけで済む。
"""
function h_transport_wall(ρ, u, Em, elev, T_wall, D, G, g, tables::FluidTables;
    ε_phase=0.05, ε_mode=0.1, ε_Re=200.0, ε_J=0.1, ε_invXtt=0.02, ε_x=0.05, f_cond=1.0, xps_mode=:P)

    # 質量流束 G は、ρ|u| で、状態変数 ρm が Newton 反復の途中で負になると負になりうる。
    # Re_LS^(-1/3) などの小数べきの底になるので、下限を守る (通常の G は 1e-2 以上)。
    G = soft_floor(G, 1.0e-6)

    em = internal_energy(Em, u, elev, g)
    P, T, x = flash_ρe(ρ, em, tables)
    # 熱伝達のブレンド用: 圧力の飽和温度基準 (過冷却液は x_pseudo<0、過熱蒸気は >1)。
    # (T 基準の pseudo_quality は、単相液でも ≈0 になり、二相の h が半分掛かっていた)
    x_pseudo = xps_mode === :P ? pseudo_quality_P(em, ρ, P, tables) : pseudo_quality(em, ρ, P, T, tables)   # :T は従来 (切り分け用)

    # --- 単相側: その場の(ρ,e)から直接、正しい相の物性が得られる ---
    μ, k, cp, Pr = transport_props(ρ, em, tables)
    Re = reynolds_number_safe(ρ, u, D, μ)
    n_exp = dittus_boelter_exponent_smooth(T_wall, T)
    h_singlephase = h_from_Nu(Nu_pipe_smooth(Re, Pr, n_exp; ε_Re=ε_Re), k, D)

    # --- 二相側: 凝縮⇔沸騰(Chen⇔Dougall-Rohsenow)をブレンド ---
    props = sat_props(T, tables)
    T_sat = Tsat(P, tables)
    P_sat_Twall = Psat(T_wall, tables)

    h_cond = f_cond * h_condense_shah_smooth(x, G, D, g, P, tables.Pcrit, props; ε_J=ε_J)   # f_cond: 凝縮の h の係数 (既定 1)
    h_chen = h_chen_smooth(T_wall, T_sat, P_sat_Twall, P, tables.Pcrit, x, G, D, props, props.h_fg; ε_invXtt=ε_invXtt)
    h_pd = h_dougall_rohsenow_smooth(x, G, D, props)
    h_boil = h_boiling_smooth(x, h_chen, h_pd; ε_x=ε_x)

    # T_wall<T_sat (s=T_sat-T_wall>0) で凝縮側へ
    h_twophase = blend_tanh(h_cond, h_boil, T_sat - T_wall, ε_mode)

    # --- 単相⇔二相の接続 (疑似乾き度による3領域ブレンド) ---
    w_high = smooth_step(x_pseudo - 1.0, ε_phase)   # →1: 単相蒸気側
    w_low = smooth_step(-x_pseudo, ε_phase)         # →1: 単相液側
    w_twophase = (1 - w_high) * (1 - w_low)

    return w_twophase * h_twophase + (1 - w_twophase) * h_singlephase
end

"""
    friction_force_transport(ρ, u, Em, elev, D, G, g, tables::FluidTables; kwargs...)

輸送管・グルーブ・CCの摩擦力の大きさ |Fw| [Pa/m] (式24-28) を1つの値として
返す (向き sign(u) は含まない。呼び出し側で掛ける)。
単相⇔二相 (Lockhart-Martinelli) の切り替えを滑らかに統合する。
"""
function friction_force_transport(ρ, u, Em, elev, D, G, g, tables::FluidTables;
    ε_phase=0.05, ε_Re=200.0)

    # 質量流束 G は、ρ|u| で、状態変数 ρm が Newton 反復の途中で負になると負になりうる。
    # Re_LS^(-1/3) などの小数べきの底になるので、下限を守る (通常の G は 1e-2 以上)。
    G = soft_floor(G, 1.0e-6)

    em = internal_energy(Em, u, elev, g)
    P, T, x = flash_ρe(ρ, em, tables)
    x_pseudo = pseudo_quality(em, ρ, P, T, tables)

    # --- 単相側 ---
    μ, k, cp, Pr = transport_props(ρ, em, tables)
    Re_sp = reynolds_number_safe(ρ, u, D, μ)
    f_sp = friction_factor_smooth(Re_sp; ε_Re=ε_Re)
    τw_sp = wall_shear_stress(f_sp, ρ, u)
    ΔP_sp = pressure_drop_single_phase(D, τw_sp)

    # --- 二相側 (Lockhart-Martinelli) ---
    props = sat_props(T, tables)
    x_safe = quality_for_correlation(x)
    G_l = G * (1 - x_safe)
    G_v = G * x_safe
    u_l = G_l / props.ρ_l
    u_v = G_v / props.ρ_v

    Re_l = reynolds_number_safe(props.ρ_l, u_l, D, props.μ_l)
    Re_v = reynolds_number_safe(props.ρ_v, u_v, D, props.μ_v)

    f_l = friction_factor_smooth(Re_l; ε_Re=ε_Re)
    f_v = friction_factor_smooth(Re_v; ε_Re=ε_Re)
    τw_l = wall_shear_stress(f_l, props.ρ_l, u_l)
    τw_v = wall_shear_stress(f_v, props.ρ_v, u_v)
    ΔP_l = pressure_drop_single_phase(D, τw_l)
    ΔP_v = pressure_drop_single_phase(D, τw_v)

    C = C_lockhart_martinelli_smooth(Re_v, Re_l; ε_Re=ε_Re)
    ΔP_2p = pressure_drop_twophase(ΔP_l, ΔP_v, C)

    # --- 単相⇔二相の接続 ---
    w_high = smooth_step(x_pseudo - 1.0, ε_phase)
    w_low = smooth_step(-x_pseudo, ε_phase)
    w_twophase = (1 - w_high) * (1 - w_low)

    return w_twophase * ΔP_2p + (1 - w_twophase) * ΔP_sp
end

"""
    h_evaporation()

蒸発熱伝達係数 (Kaya[14]を参考に設定した定数値、以前の議論参照)。
"""
h_evaporation() = 10000.0

"""
    h_ec_grv(Re, Pr, D, k; ε_Re=2e4)

蒸発器ケースからグルーブ蒸気への過熱熱伝達係数 (式45, 67)。
"""
function h_ec_grv(Re, Pr, D, k; ε_Re=2.0e4)
    return h_from_Nu(Nu_flatplate_smooth(Re, Pr; ε_Re=ε_Re), k, D)
end