# =============================================================================
# residuals_flatwick.jl (平板型の primary wick と secondary wick)
#
# 【構成】 CC流体 ─(ṁ)→ secondary wick (1ノード) ─(ṁ)→ primary wick (N層) → 蒸発界面
#   ・液管の流入は secondary wick を通らず、CCに直接入る (residual_cc! の入口そのまま)。
#   ・CCから出る液 ṁ は、secondary wick を通って primary wick の第1層に入る。
#     (secondary wick が「流体の通るノード」であり、同時にCCに対する管壁の役割も持つ)
#   ・primary wick は、高さ方向 (secondary wick 側 → 蒸発界面) の1次元 N 層。
#     蒸発器ケースから primary wick の側面への入熱は無視する。
#   ・secondary wick は常に完全に液で満たされている (気相なし、状態は E, T の2変数)。
#
# 【再利用】 primary wick の残差は円筒型の residual_wick! (residuals_wick.jl) を
# そのまま使う。違いは次の2点だけで、どちらも引数で表現できる。
#   1. 内側のノード0の温度が T_CC ではなく T_swick (secondary wick の温度)。
#      residual_wick! の第7引数 T_cc に T_swick を渡す。
#      戻り値 Q_leak は「primary wick → secondary wick の熱 G_0 (T_1 − T_swick)」になる。
#   2. 幾何 (G, V) は、平板の高さ方向の1次元に作り直す (flatwick_params)。
# secondary wick の残差は residual_swick! (本ファイル)。
#
# 【エネルギー収支】 secondary wick (温度 T_s):
#   dE_s/dt = ṁ(⟨h⟩_in − ⟨h⟩_out)           … CCから入り、primary wick へ出る液
#           + Q_cs                             … CCケース → secondary wick (熱伝導)
#           − Q_nc                             … secondary wick → CC流体 (自然対流)
#           + Q_leak                           … primary wick → secondary wick (熱伝導)
#   ⟨h⟩_in  = wm h(T_CC) + (1−wm) h(T_s),  ⟨h⟩_out = wm h(T_s) + (1−wm) h(T_1)
#   (wm = smooth_step(ṁ, ε_m)。⟨h⟩_out は primary wick 第1層の入口の風上値
#    hf[1] と同じ式なので、secondary wick が失う分と第1層が得る分が一致する。)
#   E_s = V_s [(1−φ_s) ρ_swick C_swick T_s + φ_s ρ_l_ref e_l(T_s)]    (primary wick の式46と同じ形)
#
# 【形状の読み取り (断面図)】 secondary wick は、CC流体 (半円 D_i) を囲む半円環
#   (半径方向 (D_swick−D_i)/2) と、平面側・両端面の厚さ T_swick の層から成る。
#   外寸は直径 D_swick、長さ dl + 2 T_swick の半円筒で、V_s = 外寸の体積 − CC流体の体積。
# =============================================================================

"""
    flatwick_params(geo::FlatWickGeometry, cc::FlatCCGeometry, sgeo::SecondaryWickGeometry, k_s, k_l)
        -> NamedTuple

primary wick (高さ方向 N 層) の幾何量。residual_wick! が使うフィールド
(N, V, G, por, ρ_w, C_w) を、円筒型と同じ意味で持つ。

    G[1]      … secondary wick ↔ 第1層の中心 (直列: secondary 側 T_swick/2, primary 側 dH/2)
    G[2..N]   … 層の中心どうし  k_eff A / dH
    G[N+1]    … 第N層の中心 ↔ 蒸発界面  k_eff A / (dH/2)
    A         … primary wick の水平断面積 W L (側面からの熱流は無視)
k_s, k_l は固体 (SUS304) と液の熱伝導率で、primary / secondary それぞれの空隙率で
マクスウェルの複合則 (k_eff_maxwell) から k_eff を求める。
"""
function flatwick_params(geo::FlatWickGeometry, cc::FlatCCGeometry, sgeo::SecondaryWickGeometry, k_s, k_l)
    N, dH = geo.N, geo.dH
    k_eff = k_eff_maxwell(k_s, k_l, geo.φ)
    k_eff_s = k_eff_maxwell(k_s, k_l, sgeo.φ_s)
    A = geo.W * geo.L
    V = fill(A * dH, N)
    G_ps_sec = k_eff_s * cc.A_pr_se / (cc.T_swick / 2)
    G_ps_pri = k_eff * cc.A_pr_se / (dH / 2)
    G = Vector{Float64}(undef, N + 1)
    G[1] = 1 / (1 / G_ps_sec + 1 / G_ps_pri)
    for k in 2:N
        G[k] = k_eff * A / dH
    end
    G[N+1] = k_eff * A / (dH / 2)
    return (N = N, L = geo.L, V = V, G = G, k_eff = k_eff, A = A,
            por = geo.φ, ρ_w = geo.ρ_wick, C_w = geo.Cp_wick)
end

"""
    swick_params(sgeo::SecondaryWickGeometry, cc::FlatCCGeometry, k_s, k_l) -> NamedTuple

secondary wick の幾何量・熱容量用の量と、CCケース ↔ secondary wick の熱伝導の量。

    V       … 体積 [m³] (外寸の半円筒 D_swick × (dl+2T_swick) − CC流体の半円筒)
    k_eff   … 固相の実効熱伝導率 (マクスウェル、φ_s)
    A_cs    … CCケースとの接触面積 S_swick [m²]
    t_c, t_s… CCケース壁の厚さ (D_o−D_swick)/2、secondary wick の半径方向厚さ (D_swick−D_i)/2
"""
function swick_params(sgeo::SecondaryWickGeometry, cc::FlatCCGeometry, k_s, k_l)
    k_eff = k_eff_maxwell(k_s, k_l, sgeo.φ_s)
    L_sw = cc.dl + 2 * cc.T_swick                                   # secondary wick の全長 (両端面を含む)
    # 外寸: 直径 D_swick、全長 dl + 2 T_swick の半円筒 (平面側・両端面は厚さ T_swick の層)。
    # secondary wick の体積 = 外寸の半円筒 − CC流体の半円筒 (D_i, 長さ dl)
    L_sw = cc.dl + 2 * cc.T_swick
    V_env = 0.5 * (π / 4) * cc.D_swick^2 * L_sw
    V_fluid = cc.A_cross * cc.dl
    return (V = V_env - V_fluid, k_eff = k_eff,
            por = sgeo.φ_s, ρ_w = sgeo.ρ_swick, C_w = sgeo.Cp_swick,
            A_cs = cc.S_swick, t_c = (cc.D_o - cc.D_swick) / 2, t_s = (cc.D_swick - cc.D_i) / 2)
end

"""
    Q_case_swick(T_cw, T_s, sp) -> Q [W]

CCケース (温度 T_cw) から secondary wick (温度 T_s) への熱伝導 (CCケースが熱ければ正)。
ケース壁の半厚さと secondary wick の半厚さの直列 (k_sus304(T_cw) は温度依存)。
CCケースの残差と secondary wick の残差が同じ値を使うので、エネルギーが厳密に一致する。
"""
function Q_case_swick(T_cw, T_s, sp)
    G_c = k_sus304(T_cw) * sp.A_cs / (sp.t_c / 2)
    G_s = sp.k_eff * sp.A_cs / (sp.t_s / 2)
    return (T_cw - T_s) / (1 / G_c + 1 / G_s)
end

"""
    residual_swick!(resid, du, s, sp, ṁ, T_cc, T_w1, P_w, Q_cs, Q_nc, Q_leak, tables; ρ_l_ref, ε_m=1e-8)

secondary wick (1ノード) の残差。resid, du, s は (E, T) を持つ ComponentVector (各1要素)。

  ṁ      : wick_interface の ṁ (CC → secondary wick → primary wick を通る蒸発量)
  T_cc   : CC内の液温度 (ṁ>0 のとき流入する液の温度)
  T_w1   : primary wick 第1層の温度
  Q_cs   : Q_case_swick (CCケース → secondary wick)
  Q_nc   : flatcc_natural_convection の Q_fluid (secondary wick → CC流体)
  Q_leak : residual_wick! の戻り値 (primary wick → secondary wick)
"""
function residual_swick!(resid, du, s, sp, ṁ, T_cc, T_w1, P_w, Q_cs, Q_nc, Q_leak,
    tables::FluidTables; ρ_l_ref, ε_m=1e-8)

    T = s.T[1]
    wm = smooth_step(ṁ, ε_m)
    h_cc = liquid_props(T_cc, P_w, tables).h
    lp = liquid_props(T, P_w, tables)
    h_w1 = liquid_props(T_w1, P_w, tables).h
    hf_in = wm * h_cc + (1 - wm) * lp.h
    hf_out = wm * lp.h + (1 - wm) * h_w1

    resid.E[1] = du.E[1] - (ṁ * (hf_in - hf_out) + Q_cs - Q_nc + Q_leak)
    resid.T[1] = s.E[1] - sp.V * ((1 - sp.por) * sp.ρ_w * sp.C_w * T +
                                  sp.por * ρ_l_ref * lp.e)
    return nothing
end