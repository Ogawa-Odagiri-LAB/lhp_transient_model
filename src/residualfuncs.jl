function f_mass(V, T, T_old, α, α_old, m, m_2, FLUID, dt) # ノードiにおける質量保存則，ジャンクションjからj+1までを考慮する
    ρ = calc_ρ_2p(T, α, FLUID)
    ρ_old = calc_ρ_2p(T_old, α_old, FLUID)
    return V * (ρ - ρ_old) - dt * (m - m_2)
end

function f_mom(dl, A_cross, D_i, m, m_old, T_1, P_1, α_1, P, FLUID, dt) # ジャンクションjにおける運動量保存則，ノードi-1からiまでを考慮する
    ΔP = 0.0
    if α_1 == 1.0 || α_1 == 0.0
        ΔP = calc_ΔP_1p(T_1, P_1, m, D_i, A_cross, FLUID)
    else
        ΔP = calc_ΔP_2p(T_1, α_1, m, D_i, A_cross, FLUID)
    end
    return dl * (m - m_old) - dt * (A_cross * (P_1 - P) - ΔP * dl)
end

function f_ene(V, A_cross, T, T_old, α, α_old, m, m_old, T_1, α_1, m_1, T_2, α_2, m_2, Q_wall, FLUID, dt) # ノードiにおけるエネルギー保存則，ノードi-1からノードi+1までを考慮する
    ρ = calc_ρ_2p(T, α, FLUID)
    e = calc_e_2p(T, α, FLUID)
    E = V * ρ * (e + 0.5 * (m/(ρ*A_cross))^2)
    ρ_old = calc_ρ_2p(T_old, α_old, FLUID)
    e_old = calc_e_2p(T_old, α_old, FLUID)
    E_old = V * ρ_old * (e_old + 0.5 * (m_old/(ρ_old*A_cross))^2)

    # ノード1
    h_1 = calc_h_2p(T_1, α_1, FLUID)
    ρ_1 = calc_ρ_2p(T_1, α_1, FLUID)
    E_1 = m_1 * (h_1 + 0.5 * (m_1/(ρ_1*A_cross))^2)

    # ノード2
    h_2 = calc_h_2p(T_2, α_2, FLUID)
    ρ_2 = calc_ρ_2p(T_2, α_2, FLUID)
    E_2 = m_2 * (h_2 + 0.5 * (m_2/(ρ_2*A_cross))^2)

    # Q_i
    Q_i = Q_wall

    return (E - E_old) - dt * (E_1 - E_2 + Q_i)
end

function f_state(T, P, α, FLUID)
    if α == 1.0
        f = α - 1.0
    elseif α == 0.0
        f = α - 0.0
    else
        f = T - calc_Tsat(P, FLUID)
    end
    return f
end
