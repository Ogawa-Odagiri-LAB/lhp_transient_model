# =============================================================================
# residuals_housing.jl (平板型: 周囲との熱交換の遅い経路 = 断熱材 (エアロフレックス) と PEEK ブロック)
#
# 【背景】 実験は +300 s 以降も、ループ全体の温度が 0.9 mK/s で上がり続ける (HB, CC, 蒸気管, 液管とも)。
# CC は還流の過冷度で 0.02 W/K しか熱を吸収できない (感度 50 K/W)。0.013 W (入熱の 0.3 %) の入熱の増加で
# T_cc が +0.67 K 動く。モデルは HB・蒸発器ケース・CC が周囲と熱交換しない (熱橋 Q_br だけ) ので、これが無い。
#
# 【モデル】
#   (1) 断熱材 (エアロフレックス、厚さ hi_t、HB・蒸発器ケース・CC の上側を 1 枚で包む) を 3 節点 (HB 側・蒸発器ケース側・CC 側) の
#       集中容量にする。各節点: 壁 → 節点 (厚さの半分の伝導 G_in = ins_k A/(t/2))、節点 → 周囲 (厚さの半分の伝導 + 外表面 ins_h_o の直列 G_out)。
#       熱容量 C = hi_rcp hi_t A。節点どうしは、断熱材の中の横方向の熱伝導 (hi_G_hb_ev, hi_G_ev_cc) で結ぶ。
#       時定数 ≈ C/(G_in+G_out) ≈ 数百 s … 壁が暖まると、断熱材が暖まるまで熱を取られ、その後、取る熱が減る。
#   (2) PEEK ブロック (HB の下、HB と同じ面積、高さ peek_L、下端は下板 = T_amb) を 1 次元の熱伝導 (N_PEEK 節点、表面ほど細かい等比格子) にする。
#       PEEK の拡散時間 L²/α ≈ 2.6e4 s と長いので、半無限体のように、取る熱が時間とともにゆっくり減る。
# 断熱材・PEEK の状態は T (K) で、スケールなし (hb と同じ)。hi_on = peek_on = false のときは dT/dt = 0 に固定して熱を出さない。
# =============================================================================

hi_areas(c) = (c.hi_A_hb, c.hi_A_ev, c.hi_A_cc)
hi_G_in(c, A)  = c.ins_k * A / (c.hi_t / 2)
hi_G_out(c, A) = A / ((c.hi_t / 2) / c.ins_k + 1 / c.ins_h_o)
hi_C(c, A)     = c.hi_rcp * c.hi_t * A

"PEEK の節点の厚さ [m] (HB 側から。等比 PEEK_RATIO、合計 peek_L)"
function peek_dz(c)
    dz1 = c.peek_L * (PEEK_RATIO - 1) / (PEEK_RATIO^N_PEEK - 1)
    return [dz1 * PEEK_RATIO^(k - 1) for k in 1:N_PEEK]
end

"""
    residual_housing!(r_hi, d_hi, x_hi, r_pk, d_pk, x_pk, T_hb, T_ec, T_ccw, c) -> (Q_hb, Q_ev, Q_cc)

断熱材 3 節点と PEEK の残差を書き込み、HB・蒸発器ケース・CC ケースが失う熱 [W] (失う向きが正) を返す。
"""
function residual_housing!(r_hi, d_hi, x_hi, r_pk, d_pk, x_pk, T_hb, T_ec, T_ccw, c; Q_ll_in = 0.0)
    Q_hb = 0.0; Q_ev = 0.0; Q_cc = 0.0
    if c.hi_on
        Tw = (T_hb, T_ec, T_ccw)
        for (i, A) in enumerate(hi_areas(c))
            Gi, Go = hi_G_in(c, A), hi_G_out(c, A)
            (i == 3 && c.cc_adiabatic) && (Gi = 0.0)    # 【試験】 CC ケース → 節点 3 の伝導を切る (cc_adiabatic)。節点 3 の横方向の伝導・液管の熱・外表面の損失は残る
            Qw = Gi * (Tw[i] - x_hi[i])
            # 横方向 (断熱材の中) の熱伝導: HB 側 ↔ 蒸発器ケース側 ↔ CC 側
            Ql = i == 1 ? -c.hi_G_hb_ev * (x_hi[1] - x_hi[2]) :
                 i == 2 ? c.hi_G_hb_ev * (x_hi[1] - x_hi[2]) - c.hi_G_ev_cc * (x_hi[2] - x_hi[3]) :
                          c.hi_G_ev_cc * (x_hi[2] - x_hi[3])
            r_hi[i] = hi_C(c, A) * d_hi[i] - (Qw - Go * (x_hi[i] - c.T_amb) + Ql + (i == 3 ? Q_ll_in : 0.0))
            i == 1 && (Q_hb += Qw)
            i == 2 && (Q_ev += Qw)
            i == 3 && (Q_cc += Qw)
        end
    else
        for i in eachindex(r_hi)
            r_hi[i] = d_hi[i]
        end
    end
    if c.peek_on
        dz = peek_dz(c)
        kA = c.peek_k * c.peek_A
        rA = c.peek_rcp * c.peek_A
        Qin = kA / (dz[1] / 2) * (T_hb - x_pk[1])          # HB → PEEK 第 1 節点
        Q_hb += Qin
        for k in 1:N_PEEK
            Qup = k == 1 ? Qin : kA / ((dz[k - 1] + dz[k]) / 2) * (x_pk[k - 1] - x_pk[k])
            Qdn = k == N_PEEK ? kA / (dz[k] / 2) * (x_pk[k] - c.T_amb) : kA / ((dz[k] + dz[k + 1]) / 2) * (x_pk[k] - x_pk[k + 1])
            r_pk[k] = rA * dz[k] * d_pk[k] - (Qup - Qdn)
        end
    else
        for k in eachindex(r_pk)
            r_pk[k] = d_pk[k]
        end
    end
    return (Q_hb, Q_ev, Q_cc)
end

"初期の断熱材・PEEK の温度: 壁の初期温度 (HB・蒸発器ケース・CC ケース) と周囲の間の定常の値"
function init_housing!(u0, c)
    Tw = (u0.T_wall.hb[1], u0.T_wall.evap[1], u0.T_wall.cc[1])
    # 断熱材 3 節点の定常 (横方向の熱伝導を含む 3×3 の連立)
    Gi = [hi_G_in(c, A) for A in hi_areas(c)]
    c.cc_adiabatic && (Gi[3] = 0.0)   # 【試験】 節点 3 は CC ケースとつながない (residual_housing! と同じ)
    Go = [hi_G_out(c, A) for A in hi_areas(c)]
    a, b = c.hi_G_hb_ev, c.hi_G_ev_cc
    M = [Gi[1]+Go[1]+a  -a  0.0; -a  Gi[2]+Go[2]+a+b  -b; 0.0  -b  Gi[3]+Go[3]+b]
    rhs = [Gi[i] * Tw[i] + Go[i] * c.T_amb for i in 1:3]
    u0.hi.T .= M \ rhs
    dz = peek_dz(c)
    z = 0.0
    for k in 1:N_PEEK
        zc = z + dz[k] / 2
        u0.peek.T[k] = Tw[1] - (Tw[1] - c.T_amb) * zc / c.peek_L
        z += dz[k]
    end
    return nothing
end

"断熱材・PEEK に蓄えられた熱 [J] (0 °C 基準) と、周囲へ出る熱 [W]"
function housing_energy(c_state, cnd)
    (:hi in keys(c_state)) || return 0.0
    E = 0.0
    if cnd.hi_on
        for (i, A) in enumerate(hi_areas(cnd))
            E += hi_C(cnd, A) * (c_state.hi.T[i] - 273.15)
        end
    end
    if cnd.peek_on
        dz = peek_dz(cnd)
        for k in 1:N_PEEK
            E += cnd.peek_rcp * cnd.peek_A * dz[k] * (c_state.peek.T[k] - 273.15)
        end
    end
    return E
end
function housing_heat_out(c_state, cnd)
    (:hi in keys(c_state)) || return 0.0
    Q = 0.0
    if cnd.hi_on
        for (i, A) in enumerate(hi_areas(cnd))
            Q += hi_G_out(cnd, A) * (c_state.hi.T[i] - cnd.T_amb)
        end
    end
    if cnd.peek_on
        dz = peek_dz(cnd)
        Q += cnd.peek_k * cnd.peek_A / (dz[end] / 2) * (c_state.peek.T[end] - cnd.T_amb)
    end
    return Q
end