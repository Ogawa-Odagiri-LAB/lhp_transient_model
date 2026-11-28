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
#   V_s = 半円環 + 平面側の層 (D_swick × T_swick × 全長) + 両端面の2枚。
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

    V       … 体積 [m³] (半円環 + 平面側の層 + 両端面)
    k_eff   … 固相の実効熱伝導率 (マクスウェル、φ_s)
    A_cs    … CCケースとの接触面積 S_swick [m²]
    t_c, t_s… CCケース壁の厚さ (D_o−D_swick)/2、secondary wick の半径方向厚さ (D_swick−D_i)/2
"""
function swick_params(sgeo::SecondaryWickGeometry, cc::FlatCCGeometry, k_s, k_l)
    k_eff = k_eff_maxwell(k_s, k_l, sgeo.φ_s)
    L_sw = cc.dl + 2 * cc.T_swick                                   # secondary wick の全長 (両端面を含む)
    # 形状: CC流体 (半円 D_i) を囲む半円環 (半径方向 (D_swick−D_i)/2) に、
    # 厚さ T_swick の平面側の層 (primary wick との間) と両端面を足したもの
    L_sw = cc.dl + 2 * cc.T_swick                                   # secondary wick の全長 (両端面を含む)
    V_ring = 0.5 * (π / 4) * (cc.D_swick^2 - cc.D_i^2) * cc.dl      # 半円環
    V_flat = cc.D_swick * cc.T_swick * L_sw                         # 平面側の層
    V_ends = 2 * 0.5 * (π / 4) * cc.D_swick^2 * cc.T_swick          # 両端面 (半円板 × 厚さ) 2枚
    # --- 2 節点 (A: 平面側 + 両端面 = primary wick・液に接する側、B: 半円環 = 曲面 = 蒸気が凝縮する壁 (A_vap) の側) ---
    t_ring = (cc.D_swick - cc.D_i) / 2                                # 半円環の半径方向の厚さ
    r_m    = (cc.D_i + cc.D_swick) / 4                                # 半円環の平均半径
    # A–B の熱伝導 (swick の k_eff。幾何だけで決まり、フィット係数なし。断面と経路は図面からの見積り):
    #   (1) 平面側の層 ↔ 半円環: 半円環の端 2 つ (断面 2 t_ring dl)。経路 = 平面の半幅の中心まで D_swick/4 + 半円環の 1/4 周の中心まで (π/2) r_m / 2
    #   (2) 両端面 ↔ 半円環: 半円環の両端の縁 (断面 π r_m × t_ring × 2 枚)。経路 = 端面の厚さの半分 + 軸方向の 1/4 長さ (T_swick/2 + dl/4)
    A_edge = 2 * t_ring * cc.dl
    L_edge = cc.D_swick / 4 + (π / 2) * r_m / 2
    A_rim  = 2 * π * r_m * t_ring
    L_rim  = cc.T_swick / 2 + cc.dl / 4
    G_AB   = k_eff * (A_edge / L_edge + A_rim / L_rim)
    return (V = V_ring + V_flat + V_ends, k_eff = k_eff,
            por = sgeo.φ_s, ρ_w = sgeo.ρ_swick, C_w = sgeo.Cp_swick,
            A_cs = cc.S_swick, t_c = (cc.D_o - cc.D_swick) / 2, t_s = (cc.D_swick - cc.D_i) / 2,
            V_A = V_flat + V_ends, V_B = V_ring, G_AB = G_AB,
            A_cs_A = (π / 4) * cc.D_swick^2,  t_s_A = cc.T_swick,        # 両端面 (CCケースの端の壁に接する)
            A_cs_B = π * cc.D_swick * 0.5 * cc.dl, t_s_B = t_ring)        # 曲面 (CCケースの曲面の壁に接する)
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
    Q_case_swick2(T_cw, T_A, T_B, sp) -> (Q_A, Q_B, Q)

CCケース (温度 T_cw) から secondary wick の 2 節点 (A: 両端面、B: 曲面) への熱伝導 (CCケースが熱ければ正)。
Q_case_swick と同じ直列の式を、節点ごとの面積・厚さ (A_cs_A, t_s_A / A_cs_B, t_s_B) で使う。Q = Q_A + Q_B (CCケースが失う熱)。
"""
function Q_case_swick2(T_cw, T_A, T_B, sp)
    G_c = k_sus304(T_cw) / (sp.t_c / 2)
    QA = (T_cw - T_A) * sp.A_cs_A / (1 / G_c + (sp.t_s_A / 2) / sp.k_eff)
    QB = (T_cw - T_B) * sp.A_cs_B / (1 / G_c + (sp.t_s_B / 2) / sp.k_eff)
    return (Q_A = QA, Q_B = QB, Q = QA + QB)
end

"""
    residual_swick2!(resid, du, s, sp, ṁ, T_cc, T_w1, P_w, Q_cs_A, Q_cs_B, Q_nc_A, Q_nc_B, Q_wall, Q_leak, tables; ρ_l_ref, ε_m=1e-8)

secondary wick の 2 節点の残差 (s.E, s.T は各 2 要素: [1] = A (平面側 + 両端面)、[2] = B (半円環))。
  A: dE_A/dt = ṁ(⟨h⟩_in − ⟨h⟩_out) + Q_cs_A − Q_nc_A + Q_leak + Q_AB      (液の流れ ṁ・primary wick からの熱は A を通る)
  B: dE_B/dt = Q_cs_B − Q_nc_B + Q_wall − Q_AB                           (壁面凝縮の潜熱 Q_wall は B へ)
  Q_AB = G_AB (T_B − T_A)  (B → A)

via_B = true (上面加熱 :heat_top): CC の液は底の曲面 (B) に溜まるので、液の流れは CC → B → A → primary wick の順に通る。
  B: + ṁ(⟨h⟩_in,B − ⟨h⟩_BA),   A: + ṁ(⟨h⟩_BA − ⟨h⟩_out)
  ⟨h⟩_in,B = wm h(T_CC) + (1−wm) h(T_B),  ⟨h⟩_BA = wm h(T_B) + (1−wm) h(T_A),  ⟨h⟩_out は従来と同じ
  (B が失う分と A が得る分が同じ式なので、エネルギーは厳密に保存する)。via_B = false (既定) は従来どおり。

Q_wall_A (上面加熱): 平面側 (節点 A) での壁面凝縮の潜熱。A に入る (既定 0 = 従来)。Q_wall は曲面 (B) の分。
"""
function residual_swick2!(resid, du, s, sp, ṁ, T_cc, T_w1, P_w, Q_cs_A, Q_cs_B, Q_nc_A, Q_nc_B, Q_wall, Q_leak,
    tables::FluidTables; ρ_l_ref, ε_m=1e-8, via_B::Bool=false, Q_wall_A=0.0)

    T_A = s.T[1]
    T_B = s.T[2]
    wm = smooth_step(ṁ, ε_m)
    h_cc = liquid_props(T_cc, P_w, tables).h
    lpA = liquid_props(T_A, P_w, tables)
    lpB = liquid_props(T_B, P_w, tables)
    h_w1 = liquid_props(T_w1, P_w, tables).h
    hf_out = wm * lpA.h + (1 - wm) * h_w1
    Q_AB = sp.G_AB * (T_B - T_A)

    if via_B
        hf_inB = wm * h_cc + (1 - wm) * lpB.h
        hf_BA  = wm * lpB.h + (1 - wm) * lpA.h
        resid.E[1] = du.E[1] - (ṁ * (hf_BA - hf_out) + Q_cs_A - Q_nc_A + Q_leak + Q_AB + Q_wall_A)
        resid.E[2] = du.E[2] - (ṁ * (hf_inB - hf_BA) + Q_cs_B - Q_nc_B + Q_wall - Q_AB)
    else
        hf_in = wm * h_cc + (1 - wm) * lpA.h
        resid.E[1] = du.E[1] - (ṁ * (hf_in - hf_out) + Q_cs_A - Q_nc_A + Q_leak + Q_AB + Q_wall_A)
        resid.E[2] = du.E[2] - (Q_cs_B - Q_nc_B + Q_wall - Q_AB)
    end
    resid.T[1] = s.E[1] - sp.V_A * ((1 - sp.por) * sp.ρ_w * sp.C_w * T_A + sp.por * ρ_l_ref * lpA.e)
    resid.T[2] = s.E[2] - sp.V_B * ((1 - sp.por) * sp.ρ_w * sp.C_w * T_B + sp.por * ρ_l_ref * lpB.e)
    return nothing
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

"""
    pore_liquid_mass(x_wick, x_swick, wp, sp, tables; ρ_l_ref) -> m [kg]

primary wick (各層) と secondary wick (A, B) の孔の中の液の質量 Σ V φ ρ_l,sat(T)。孔の体積は一定なので、温度が上がれば密度が下がり、質量が減る
(押し出された液は CC の液へ戻る。pore_liquid_mass_rate)。診断用 (summary の m_pore)。ρ_l_ref は使わない (引数は互換のため)。
"""
function pore_liquid_mass(T_wick, T_swick, wp, sp, tables::FluidTables)
    m = 0.0
    for i in eachindex(T_wick)
        m += wp.V[i] * wp.por * tables.ρ_l_sat(T_wick[i])
    end
    m += sp.V_A * sp.por * tables.ρ_l_sat(T_swick[1]) + sp.V_B * sp.por * tables.ρ_l_sat(T_swick[2])
    return m
end

"""
    pore_liquid_mass_rate(d_wick, x_wick, d_swick, x_swick, wp, sp, tables; ρ_l_ref, δT=0.05) -> dm/dt [kg/s]

孔内の液の質量の変化率 Σ V φ (dρ_l,sat/dT) dT/dt (昇温で負)。dT/dt は、エネルギー E の時間変化 dE/dt を熱容量 ∂E/∂T で割って得る
(E = V[(1−φ) ρ_w C_w T + φ ρ_l_ref e_l(T)] の ∂E/∂T = V[(1−φ) ρ_w C_w + φ ρ_l_ref c_l])。
CC の液の質量の式には、CC から wick へ出る流量を ṁ_evap + dm/dt とする (孔内の液が増えれば CC から出る、減れば CC へ戻る)。
近似: wick のエネルギーの式の ρ_l_ref は一定のまま (押し出される液が運ぶエネルギーは、ΔT ≲ 0.1 K では µW 以下)。
"""
function pore_liquid_mass_rate(d_wick, x_wick, d_swick, x_swick, wp, sp, tables::FluidTables; ρ_l_ref, δT = 0.05)
    dρ(T) = (tables.ρ_l_sat(T + δT) - tables.ρ_l_sat(T - δT)) / (2δT)
    dm = 0.0
    for i in eachindex(x_wick.T)
        T = x_wick.T[i]
        C = wp.V[i] * ((1 - wp.por) * wp.ρ_w * wp.C_w + wp.por * ρ_l_ref * tables.cp_l_sat(T))
        dm += wp.V[i] * wp.por * dρ(T) * d_wick.E[i] / C
    end
    for (j, V) in ((1, sp.V_A), (2, sp.V_B))
        T = x_swick.T[j]
        C = V * ((1 - sp.por) * sp.ρ_w * sp.C_w + sp.por * ρ_l_ref * tables.cp_l_sat(T))
        dm += V * sp.por * dρ(T) * d_swick.E[j] / C
    end
    return dm
end