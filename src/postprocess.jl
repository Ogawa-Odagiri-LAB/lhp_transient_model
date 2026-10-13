# =============================================================================
# postprocess.jl (後処理: エネルギー収支)
#
# solve_transient_lhp! の結果 (sol, p) から、系全体のエネルギー収支を確認する。
#
#   ΔE = E(t_end) − E(t_init)           … 蓄積したエネルギーの変化
#   Q_in  = ∫ Q_apply dt                 … 印加熱
#   Q_out = ∫ (周囲への放熱 + 凝縮器プレートへの放熱) dt
#   残差 = ΔE − (Q_in − Q_out)
#
# 流体のエネルギーは保存形 (式3) の離散化で、面のフラックスを隣り合うセルが共有する
# ので、厳密に保存する。残差は、保存間隔 (saveat) の粗さによる台形則の積分誤差と、
# IDA の許容誤差による。残差が入力に対して数 % を超えるなら、どこかでエネルギーが
# 生成・消滅している。
#
# 【蓄積エネルギー E の定義】
#   流体: Σ ρm Em V      (Em = e + u²/2 + g z。保存形の左辺そのもの)
#   ウィック: Σ E_i       (式49)
#   壁: m ∫_{T_ref}^{T} c(T') dT'   (c は温度の関数なので、C(T)·T とは違う。
#       C(T)·T だと、温度が変わると T0 ΔC の誤差が出る: 凝縮器の壁が約9 K 冷えるとき約35 J)
# 基準温度 T_ref は、差だけが意味を持つので、何でもよい。
# =============================================================================

"""
    wall_internal_energy(T, V_w; T_ref=273.15, n=20) -> E [J]

SUS304 の壁 (体積 V_w) の内部エネルギー m ∫ c dT (台形則、n 分割)。
"""
function wall_internal_energy(T, V_w; T_ref=273.15, n::Int=20)
    Ts = range(T_ref, T, length=n + 1)
    f = [ρ_sus304(t) * c_sus304(t) for t in Ts]
    Iint = sum(0.5 * (f[k] + f[k+1]) * (Ts[k+1] - Ts[k]) for k in 1:n)
    return V_w * Iint
end

"""
    stored_energy(c, p) -> E [J]

状態 c (ComponentArray) に蓄えられている、系全体のエネルギー。
"""
function stored_energy(c, p::AbstractLHPParams)
    E = 0.0
    for (comp, pc) in ((:grv, p.p_grv), (:vl, p.p_vl), (:cond, p.p_cond), (:ll, p.p_ll), (:cc, p.p_cc))
        cm = getproperty(c, comp)
        E += sum(cm.ρm .* cm.Em) * pc.V
    end
    E += sum(c.wick.E)
    :swick in keys(c) && (E += sum(c.swick.E))   # 平板型: secondary wick
    for (wall, wp) in ((:vl, p.wp_vl), (:cond, p.wp_cond), (:ll, p.wp_ll))
        E += sum(wall_internal_energy(T, wp.V_w) for T in getproperty(c.T_wall, wall))
    end
    E += wall_internal_energy(c.T_wall.cc[1], p.cw.V_w)
    E += wall_internal_energy(c.T_wall.evap[1], p.ec.V_w)
    return E
end

"""
    heat_out(c, p) -> Q [W]

系から出ていく熱。蒸気管・液管の壁→周囲 (自然対流)、凝縮器の壁→凝縮器プレート。
CCケース、蒸発器ケースは断熱。Q_amb_* は「壁が受け取る向き」が正なので、符号を反転する。
残差 (residuals_wall.jl) と同じ関数を使うので、残差が仮定する放熱と一致する。
"""
function heat_out(c, p::AbstractLHPParams)
    cnd = p.cond
    Q = 0.0
    for (wall, wp) in ((:vl, p.wp_vl), (:ll, p.wp_ll))
        for T in getproperty(c.T_wall, wall)
            Q -= Q_amb_pipe(T, wp, cnd.T_amb)
        end
    end
    for T in c.T_wall.cond
        Q -= Q_amb_condenser(T, p.wp_cond, cnd.T_cond; h_con = cnd.h_con)
    end
    return Q
end

"""
    energy_balance(sol, p; verbose=true) -> NamedTuple

系全体のエネルギー収支。保存点 (condition.saveat) の値を、台形則で積分する。
戻り値: (ΔE, Q_in, Q_out, residual, E, Q_out_t)。E, Q_out_t は、保存点ごとの時系列。
"""
function energy_balance(sol, p::AbstractLHPParams; verbose::Bool=true)
    n = length(sol.t)
    E = [stored_energy(components_at(sol, p, k), p) for k in 1:n]
    Q_in_t = [applied_heat(t, p.cond) for t in sol.t]
    Q_out_t = [heat_out(components_at(sol, p, k), p) for k in 1:n]
    trapz(f) = sum(0.5 * (f[k] + f[k+1]) * (sol.t[k+1] - sol.t[k]) for k in 1:n-1)
    ΔE = E[end] - E[1]
    Q_in = trapz(Q_in_t)
    Q_out = trapz(Q_out_t)
    residual = ΔE - (Q_in - Q_out)
    if verbose
        println("エネルギー収支 (t = ", sol.t[1], " 〜 ", sol.t[end], " s):")
        println("  蓄積したエネルギーの変化 ΔE            = ", round(ΔE, digits=2), " J")
        println("  印加熱の積分 ∫Q_apply dt               = ", round(Q_in, digits=2), " J")
        println("  放熱の積分 ∫Q_out dt (周囲+凝縮器プレート) = ", round(Q_out, digits=2), " J")
        println("  残差 ΔE − (∫Q_in − ∫Q_out)             = ", round(residual, digits=2), " J",
                "  (入力に対して ", round(100 * residual / max(Q_in, eps()), digits=2), " %)")
    end
    return (ΔE = ΔE, Q_in = Q_in, Q_out = Q_out, residual = residual, E = E, Q_out_t = Q_out_t)
end