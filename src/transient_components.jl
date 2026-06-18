function calc_transient_wick!(f, var, var_old, geo, condition)
    FLUID = condition.fluid
    h_ev = condition.h_ev
    dt = condition.dt
    Φ = 0.56
    
    T_bay_n = var.bay.T[geo.bay.N]
    T_in = var.wick.T_in # 独立変数
    T_in_old = var_old.wick.T_in
    T_out = calc_Tsat(var.grv.P[1], FLUID)
    T_out_old = calc_Tsat(var_old.grv.P[1], FLUID)
    T_ec = var.grv.T_wall[1]
    m_wick = var.wick.m_wick # 独立変数

    A_ec_wick = geo.grv.A_ec_wick
    A_mid = geo.wick.A_mid
    A_in = geo.wick.A_in
    V_in = geo.wick.V_in
    V_out = geo.wick.V_out
    dr = (geo.wick.D_out - geo.wick.D_in) / 2
    G_wick = calc_k_eff(T_out, FLUID) * A_mid / dr
    G_core_wick = k_sus304(T_bay_n) * Φ * A_in / dr
    
    h_bay_n = calc_h_1p(T_bay_n, var.bay.P[geo.bay.N], FLUID)
    h_wick_in = calc_h_2p(T_in, 0.0, FLUID)
    h_wick_out = calc_h_2p(T_out, 0.0, FLUID)
    E_in = calc_wick_energy(T_in, V_in, FLUID)
    E_in_old = calc_wick_energy(T_in_old, V_in, FLUID)
    E_out = calc_wick_energy(T_out, V_out, FLUID)
    E_out_old = calc_wick_energy(T_out_old, V_out, FLUID)
    λ = calc_latent_heat(T_out, FLUID)

    f[1] = E_out - E_out_old - dt * (h_ev * A_ec_wick * (T_ec - T_out) - G_wick * (T_out - T_in) - m_wick * (λ + h_wick_out - h_wick_in)) # m_wickのための残差方程式
    f[2] = E_in - E_in_old - dt * (m_wick * h_bay_n - m_wick * h_wick_out + G_wick * (T_out - T_in) + G_core_wick * (T_bay_n - T_in)) # T_inのための残差方程式

end

function calc_transient_cc!(f, var, var_old, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt

    if var.cc.α_cc > 1.0
        var.cc.α_cc = 1.0
    elseif var.cc.α_cc < 0.0
        var.cc.α_cc = 0.0
    end

    T_cc = var.cc.T_cc
    T_bay = var.bay.T[geo.bay.N]
    T_cc_old = var_old.cc.T_cc
    T_cc_case = var.cc.T_wall
    P_cc = calc_Psat(T_cc, FLUID)
    α_cc = var.cc.α_cc
    α_cc_old = var_old.cc.α_cc
    m_bay = var.bay.m[geo.bay.N]
    m_wick = var.wick.m_wick

    V = geo.cc.V
    S_inner = geo.cc.S_inner

    ρ_cc = calc_ρ_2p(T_cc, α_cc, FLUID)
    ρ_cc_old = calc_ρ_2p(T_cc_old, α_cc_old, FLUID)
    c_p = calc_cp_2p(T_cc, α_cc, FLUID)
    c_p_old = calc_cp_2p(T_cc_old, α_cc_old, FLUID)

    E = V * ρ_cc * c_p * T_cc
    E_old = V * ρ_cc_old * c_p_old * T_cc_old
    h_bay = calc_h_2p(var.bay.T[geo.bay.N], var.bay.α[geo.bay.N], FLUID)
    h_wick = calc_h_2p(var.wick.T_in, 0.0, FLUID)

    G_cc = calc_k_2p(T_cc, α_cc, FLUID) * S_inner / (geo.cc.D_o - geo.cc.D_i)
    G_cc_bay = calc_k_2p(T_cc, α_cc, FLUID) * (π * geo.bay.D_o * geo.bay.L) / geo.bay.D_o 
    
    f[1] = (E - E_old) - dt * (m_bay * h_bay - m_wick * h_wick + G_cc * (T_cc_case - T_cc) + G_cc_bay * (T_bay - T_cc)) # T_ccのための残差方程式
    f[2] = V * (ρ_cc - ρ_cc_old) - dt * (m_bay - m_wick) # α_ccのための残差方程式

end

function calc_transient_bay!(f_mass, f_mom, f_ene, f_state, var, var_old, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt

    # 二相流チェック
    for i in 1:geo.bay.N
        check = TPcheck(var.bay.T[i], var.bay.P[i], FLUID)
        if check == :tp
            var.bay.T[i] = calc_Tsat(var.bay.P[i], FLUID)
            if var.bay.α[i] > 1.0
                var.bay.α[i] = 1.0
            elseif var.bay.α[i] < 0.0
                var.bay.α[i] = 0.0
            end
        elseif check == :vapor
            var.bay.α[i] = 1.0
        elseif check == :liquid
            var.bay.α[i] = 0.0
        else
            error("invalid phase")
        end
    end

    T_ll = var.ll.T[geo.ll.N]
    T_cc = var.cc.T_cc
    m_ll = var.ll.m[geo.ll.N]
    m_bay = var.bay.m[geo.bay.N]
    P_ll = var.ll.P[geo.ll.N]
    P_cc = calc_Psat(T_cc, FLUID)
    ρ_1p_ll = calc_ρ_1p(T_ll, P_ll, FLUID)
    ρ_2p_ll = calc_ρ_2p(T_ll, var.ll.α[geo.ll.N], FLUID)
    ρ_1p_cc = calc_ρ_1p(T_cc, P_cc, FLUID)
    ρ_2p_cc = calc_ρ_2p(T_cc, var.cc.α_cc, FLUID)
    h_1p_cc = calc_h_1p(T_cc, P_cc, FLUID)
    h_2p_cc = calc_h_2p(T_cc, var.cc.α_cc, FLUID)    

    dl = geo.bay.dl
    D_i = geo.bay.D_i
    D_o = geo.bay.D_o
    A_cross = geo.bay.A_cross
    S_inner = geo.bay.S_inner
    V = geo.bay.V

    check::Symbol = :init

    ρ_1p(idx) = calc_ρ_1p(var.bay.T[idx], var.bay.P[idx], FLUID)
    ρ_2p(idx) = calc_ρ_2p(var.bay.T[idx], var.bay.α[idx], FLUID)
    ρ_1p_old(idx) = calc_ρ_1p(var_old.bay.T[idx], var_old.bay.P[idx], FLUID)
    ρ_2p_old(idx) = calc_ρ_2p(var_old.bay.T[idx], var_old.bay.α[idx], FLUID)
    h_1p(idx) = calc_h_1p(var.bay.T[idx], var.bay.P[idx], FLUID)
    h_2p(idx) = calc_h_2p(var.bay.T[idx], var.bay.α[idx], FLUID)
    e_1p(idx) = calc_e_1p(var.bay.T[idx], var.bay.P[idx], FLUID)
    e_2p(idx) = calc_e_2p(var.bay.T[idx], var.bay.α[idx], FLUID)
    e_1p_old(idx) = calc_e_1p(var_old.bay.T[idx], var_old.bay.P[idx], FLUID)
    e_2p_old(idx) = calc_e_2p(var_old.bay.T[idx], var_old.bay.α[idx], FLUID)
    h_wall_f_1p(idx) = calc_h_wall_f_1p(var.bay.T[idx], var.bay.P[idx], var.bay.m[idx], D_i, A_cross, FLUID)
    h_wall_f_2p(idx) = calc_h_wall_f_2p(var.bay.T[idx], var.bay.α[idx], var.bay.m[idx], D_i, A_cross, FLUID)
    τ_1p(idx) = calc_τ_wall(var.bay.T[idx], var.bay.P[idx], var.bay.m[idx], D_i, A_cross, FLUID)
    ΔP_1p(idx) = τ_1p(idx) * S_inner
    ΔP_2p(idx) = calc_ΔP_2p(var.bay.T[idx], var.bay.α[idx], var.bay.m[idx], D_i, A_cross, FLUID)
    E_1p(idx) = V * ρ_1p(idx) * (e_1p(idx) + 0.5 * (var.bay.m[idx] / (ρ_1p(idx) * A_cross))^2)
    E_2p(idx) = V * ρ_2p(idx) * (e_2p(idx) + 0.5 * (var.bay.m[idx] / (ρ_2p(idx) * A_cross))^2)
    E_1p_old(idx) = V * ρ_1p_old(idx) * (e_1p_old(idx) + 0.5 * (var.bay.m[idx] / (ρ_1p_old(idx) * A_cross))^2)
    E_2p_old(idx) = V * ρ_2p_old(idx) * (e_2p_old(idx) + 0.5 * (var.bay.m[idx] / (ρ_2p_old(idx) * A_cross))^2)
    Q_1p(idx) = h_wall_f_1p(idx) * S_inner * (T_cc - var.bay.T[idx]) # バイヨネット管の管壁は十分に薄く，Tccと同じ温度になっているとした．
    Q_2p(idx) = h_wall_f_2p(idx) * S_inner * (T_cc - var.bay.T[idx])

    # 残差計算
    if var_old.bay.α[1] < 1.0 && var_old.bay.α[1] > 0.0
        f_mass[1]   = V * (ρ_2p(1) - ρ_2p_old(1)) - dt * (var.bay.m[1] - var.bay.m[2])
        f_mom[1]    = dl * (var.bay.m[1] - var_old.bay.m[1]) - dt * ( (m_ll^2)/(ρ_2p_ll*A_cross) - ((var.bay.m[1])^2)/(ρ_2p(1)*A_cross) + A_cross * (P_ll - var.bay.P[1]) - dl * ΔP_2p(1))
        f_ene[1]    = E_2p(1) - E_2p_old(1) - dt * (var.bay.m[1] * (h_2p(1) + 0.5 * (var.bay.m[1]/(ρ_2p(1)*A_cross))^2) - var.bay.m[2] * (h_2p(2) + 0.5 * (var.bay.m[2]/(ρ_2p(2)*A_cross))^2) + Q_2p(1))
        f_state[1]  = var.bay.T[1] - calc_Tsat(var.bay.P[1], FLUID)
    else
        f_mass[1]   = V * (ρ_1p(1) - ρ_1p_old(1)) - dt * (var.bay.m[1] - var.bay.m[2])
        f_mom[1]    = dl * (var.bay.m[1] - var_old.bay.m[1]) - dt * ( (m_ll^2)/(ρ_1p_ll*A_cross) - ((var.bay.m[1])^2)/(ρ_1p(1)*A_cross) + A_cross * (P_ll - var.bay.P[1]) - dl * ΔP_1p(1))
        f_ene[1]    = E_1p(1) - E_1p_old(1) - dt * (var.bay.m[1] * (h_1p(1) + 0.5 * (var.bay.m[1]/(ρ_1p(1)*A_cross))^2) - var.bay.m[2] * (h_1p(2) + 0.5 * (var.bay.m[2]/(ρ_1p(2)*A_cross))^2) + Q_1p(1))
        if var_old.bay.α[1] == 1.0
            f_state[1] = var.bay.α[1] - 1.0
        elseif var_old.bay.α[1] == 0.0
            f_state[1] = var.bay.α[1] - 0.0
        end
    end

    for i in 2:geo.bay.N-1
        if var_old.bay.α[i] < 1.0 && var_old.bay.α[i] > 0.0 # 二相流
            f_mass[i]   = V * (ρ_2p(i) - ρ_2p_old(i)) - dt * (var.bay.m[i] - var.bay.m[i+1])
            f_mom[i]    = dl * (var.bay.m[i] - var_old.bay.m[i]) - dt * ( ((var.bay.m[i-1])^2)/(ρ_2p(i-1)*A_cross) - ((var.bay.m[i])^2)/(ρ_2p(i)*A_cross) + A_cross * (var.bay.P[i-1] - var.bay.P[i]) - dl * ΔP_2p(i))
            f_ene[i]    = E_2p(i) - E_2p_old(i) - dt * (var.bay.m[i] * (h_2p(i) + 0.5 * (var.bay.m[i]/(ρ_2p(i)*A_cross))^2) - var.bay.m[i+1] * (h_2p(i+1) + 0.5 * (var.bay.m[i+1]/(ρ_2p(i+1)*A_cross))^2) + Q_2p(i))
            f_state[i]  = var.bay.T[i] - calc_Tsat(var.bay.P[i], FLUID)
        else # 単相流
            f_mass[i]   = V * (ρ_1p(i) - ρ_1p_old(i)) - dt * (var.bay.m[i] - var.bay.m[i+1])
            f_mom[i]    = dl * (var.bay.m[i] - var_old.bay.m[i]) - dt * ( ((var.bay.m[i-1])^2)/(ρ_1p(i-1)*A_cross) - ((var.bay.m[i])^2)/(ρ_1p(i)*A_cross) + A_cross * (var.bay.P[i-1] - var.bay.P[i]) - dl * ΔP_1p(i))
            f_ene[i]    = E_1p(i) - E_1p_old(i) - dt * (var.bay.m[i] * (h_1p(i) + 0.5 * (var.bay.m[i]/(ρ_1p(i)*A_cross))^2) - var.bay.m[i+1] * (h_1p(i+1) + 0.5 * (var.bay.m[i+1]/(ρ_1p(i+1)*A_cross))^2) + Q_1p(i))
            if var_old.bay.α[i] == 1.0
                f_state[i] = var.bay.α[i] - 1.0
            elseif var_old.bay.α[i] == 0.0
                f_state[i] = var.bay.α[i] - 0.0
            end
        end
    end

    if var_old.bay.α[geo.bay.N] < 1.0 && var_old.bay.α[geo.bay.N] > 0.0
        f_mass[geo.bay.N] = V * (ρ_2p(geo.bay.N) - ρ_2p_old(geo.bay.N)) - dt * (var.bay.m[geo.bay.N] - m_bay)
        f_mom[geo.bay.N]  = dl * (var.bay.m[geo.bay.N] - var_old.bay.m[geo.bay.N]) - dt * ( (var.bay.m[geo.bay.N-1]^2)/(ρ_2p(geo.bay.N-1)*A_cross) - ((var.bay.m[geo.bay.N])^2)/(ρ_2p(geo.bay.N)*A_cross) + A_cross * (var.bay.P[geo.bay.N-1] - var.bay.P[geo.bay.N]) - dl * ΔP_2p(geo.bay.N))
        f_ene[geo.bay.N]  = E_2p(geo.bay.N) - E_2p_old(geo.bay.N) - dt * (var.bay.m[geo.bay.N] * (h_2p(geo.bay.N) + 0.5 * (var.bay.m[geo.bay.N]/(ρ_2p(geo.bay.N)*A_cross))^2) - m_bay * (h_2p_cc + 0.5 * (m_bay/(ρ_2p_cc*A_cross))^2) + Q_2p(geo.bay.N))
        f_state[geo.bay.N]  = var.bay.T[geo.bay.N] - calc_Tsat(var.bay.P[geo.bay.N], FLUID)
    else
        f_mass[geo.bay.N] = V * (ρ_1p(geo.bay.N) - ρ_1p_old(geo.bay.N)) - dt * (var.bay.m[geo.bay.N] - m_bay)
        f_mom[geo.bay.N]  = dl * (var.bay.m[geo.bay.N] - var_old.bay.m[geo.bay.N]) - dt * ( (var.bay.m[geo.bay.N-1]^2)/(ρ_1p(geo.bay.N-1)*A_cross) - ((var.bay.m[geo.bay.N])^2)/(ρ_1p(geo.bay.N)*A_cross) + A_cross * (var.bay.P[geo.bay.N-1] - var.bay.P[geo.bay.N]) - dl * ΔP_1p(geo.bay.N))
        f_ene[geo.bay.N]  = E_1p(geo.bay.N) - E_1p_old(geo.bay.N) - dt * (var.bay.m[geo.bay.N] * (h_1p(geo.bay.N) + 0.5 * (var.bay.m[geo.bay.N]/(ρ_1p(geo.bay.N)*A_cross))^2) - m_bay * (h_1p_cc + 0.5 * (m_bay/(ρ_1p_cc*A_cross))^2) + Q_1p(geo.bay.N))
        if var_old.bay.α[geo.bay.N] == 1.0
            f_state[geo.bay.N] = var.bay.α[geo.bay.N] - 1.0
        elseif var_old.bay.α[geo.bay.N] == 0.0
            f_state[geo.bay.N] = var.bay.α[geo.bay.N] - 0.0
        end
    end

end

function calc_transient_grv!(f_mass, f_mom, f_ene, f_state, var, var_old, geo, condition) # グルーブでは1つの溝あたりで計算していることに注意
    FLUID = condition.fluid
    dt = condition.dt

    # 完全に気相であることを仮定．
    for i in 1:geo.grv.N
        var.grv.α[i] = 1.0
    end
    
    n_grv = geo.grv.n_grv
    dl = geo.grv.dl
    width = geo.grv.width
    height = geo.grv.height
    A_grv_cross = geo.grv.A_grv_cross
    v_grv = width * height * dl

    T_wick = calc_Tsat(var.grv.P[1], FLUID)
    T_vl   = var.vl.T[1]
    m_wick = var.wick.m_wick / n_grv
    m_vl   = var.vl.m[1] / n_grv
    P_wick = var.grv.P[1]
    P_vl   = var.vl.P[1]
    ρ_1p_wick = calc_ρ_1p(T_wick, P_wick, FLUID)
    ρ_1p_vl   = calc_ρ_1p(T_vl, P_vl, FLUID)
    h_1p_vl   = calc_h_1p(T_vl, P_vl, FLUID)

    # check::Symbol = :init

    ρ_1p(idx) = calc_ρ_1p(var.grv.T[idx], var.grv.P[idx], FLUID)
    ρ_1p_old(idx) = calc_ρ_1p(var_old.grv.T[idx], var_old.grv.P[idx], FLUID)
    h_1p(idx) = calc_h_1p(var.grv.T[idx], var.grv.P[idx], FLUID)
    e_1p(idx) = calc_e_1p(var.grv.T[idx], var.grv.P[idx], FLUID)
    e_1p_old(idx) = calc_e_1p(var_old.grv.T[idx], var_old.grv.P[idx], FLUID)
    h_ec_gr(idx) = calc_h_ec_gr(var.grv.T[idx], var.grv.P[idx], var.grv.m[idx], height, A_grv_cross, geo.grv.L, FLUID)
    τ_1p(idx) = calc_τ_wall(var.grv.T[idx], var.grv.P[idx], var.grv.m[idx], height, A_grv_cross, FLUID)
    ΔP_1p(idx) = τ_1p(idx) * (width + height) * 2
    E_1p(idx) = v_grv * ρ_1p(idx) * (e_1p(idx) + 0.5 * (var.grv.m[idx] / (ρ_1p(idx) * A_grv_cross))^2)
    E_1p_old(idx) = v_grv * ρ_1p_old(idx) * (e_1p_old(idx) + 0.5 * (var.grv.m[idx] / (ρ_1p_old(idx) * A_grv_cross))^2)
    Q_1p(idx) = h_ec_gr(idx) * (width * dl) * (var.grv.T_wall[idx] - var.grv.T[idx])

    

    # 残差計算
    # 境界条件（前）
    f_mass[1]   = v_grv * (ρ_1p(1) - ρ_1p_old(1)) - dt * (var.grv.m[1] - var.grv.m[2])
    f_mom[1]    = dl * (var.grv.m[1] - var_old.grv.m[1]) - dt * ( (m_wick^2)/(ρ_1p_wick*A_grv_cross) - ((var.grv.m[1])^2)/(ρ_1p(1)*A_grv_cross) + A_grv_cross * (P_wick - var.grv.P[1]) - dl * ΔP_1p(1))
    f_ene[1]    = E_1p(1) - E_1p_old(1) - dt * (var.grv.m[1] * (h_1p(1) + 0.5 * (var.grv.m[1]/(ρ_1p(1)*A_grv_cross))^2) - var.grv.m[2] * (h_1p(2) + 0.5 * (var.grv.m[2]/(ρ_1p(2)*A_grv_cross))^2) + Q_1p(1))
    f_state[1] = var.grv.α[1] - 1.0


    for i in 2:geo.grv.N-1
        f_mass[i]   = v_grv * (ρ_1p(i) - ρ_1p_old(i)) - dt * (var.grv.m[i] - var.grv.m[i+1])
        f_mom[i]    = dl * (var.grv.m[i] - var_old.grv.m[i]) - dt * ( ((var.grv.m[i-1])^2)/(ρ_1p(i-1)*A_grv_cross) - ((var.grv.m[i])^2)/(ρ_1p(i)*A_grv_cross) + A_grv_cross * (var.grv.P[i-1] - var.grv.P[i]) - dl * ΔP_1p(i))
        f_ene[i]    = E_1p(i) - E_1p_old(i) - dt * (var.grv.m[i] * (h_1p(i) + 0.5 * (var.grv.m[i]/(ρ_1p(i)*A_grv_cross))^2) - var.grv.m[i+1] * (h_1p(i+1) + 0.5 * (var.grv.m[i+1]/(ρ_1p(i+1)*A_grv_cross))^2) + Q_1p(i))
        f_state[i] = var.grv.α[i] - 1.0
    end

    # 境界条件（後）
    f_mass[geo.grv.N] = v_grv * (ρ_1p(geo.grv.N) - ρ_1p_old(geo.grv.N)) - dt * (var.grv.m[geo.grv.N] - m_vl)
    f_mom[geo.grv.N]  = dl * (var.grv.m[geo.grv.N] - var_old.grv.m[geo.grv.N]) - dt * ( (var.grv.m[geo.grv.N-1]^2)/(ρ_1p(geo.grv.N-1)*A_grv_cross) - ((var.grv.m[geo.grv.N])^2)/(ρ_1p(geo.grv.N)*A_grv_cross) + A_grv_cross * (var.grv.P[geo.grv.N-1] - var.grv.P[geo.grv.N]) - dl * ΔP_1p(geo.grv.N))
    f_ene[geo.grv.N]  = E_1p(geo.grv.N) - E_1p_old(geo.grv.N) - dt * (var.grv.m[geo.grv.N] * (h_1p(geo.grv.N) + 0.5 * (var.grv.m[geo.grv.N]/(ρ_1p(geo.grv.N)*A_grv_cross))^2) - m_vl * (h_1p_vl + 0.5 * (m_vl/(ρ_1p_vl*A_grv_cross))^2) + Q_1p(geo.grv.N))
    f_state[geo.grv.N] = var.grv.α[geo.grv.N] - 1.0
end

function calc_transient_vl!(f_mass, f_mom, f_ene, f_state, var, var_old, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt

    # 二相流チェック
    for i in 1:geo.vl.N
        check = TPcheck(var.vl.T[i], var.vl.P[i], FLUID)
        if check == :tp
            var.vl.T[i] = calc_Tsat(var.vl.P[i], FLUID)
            if var.vl.α[i] > 1.0
                var.vl.α[i] = 1.0
            elseif var.vl.α[i] < 0.0
                var.vl.α[i] = 0.0
            end
        elseif check == :vapor
            var.vl.α[i] = 1.0
        elseif check == :liquid
            var.vl.α[i] = 0.0
        else
            error("invalid phase")
        end
    end

    T_grv  = var.grv.T[geo.grv.N]
    T_cond = var.cond.T[1]
    m_grv  = var.grv.m[geo.grv.N] * geo.grv.n_grv # グルーブでは溝あたりで計算しているため，溝の数だけ倍する必要がある．
    m_cond = var.cond.m[1]
    P_grv  = var.grv.P[geo.grv.N]
    P_cond = var.cond.P[1]
    ρ_1p_grv  = calc_ρ_1p(T_grv, P_grv, FLUID)
    ρ_2p_grv  = calc_ρ_2p(T_grv, 1.0, FLUID)
    ρ_1p_cond = calc_ρ_1p(T_cond, P_cond, FLUID)
    ρ_2p_cond = calc_ρ_2p(T_cond, var.cond.α[1], FLUID)
    h_1p_cond = calc_h_1p(T_cond, P_cond, FLUID)
    h_2p_cond = calc_h_2p(T_cond, var.cond.α[1], FLUID)    

    dl = geo.vl.dl
    D_i = geo.vl.D_i
    D_o = geo.vl.D_o
    A_cross = geo.vl.A_cross
    S_inner = geo.vl.S_inner
    V = geo.vl.V

    check::Symbol = :init

    ρ_1p(idx) = calc_ρ_1p(var.vl.T[idx], var.vl.P[idx], FLUID)
    ρ_2p(idx) = calc_ρ_2p(var.vl.T[idx], var.vl.α[idx], FLUID)
    ρ_1p_old(idx) = calc_ρ_1p(var_old.vl.T[idx], var_old.vl.P[idx], FLUID)
    ρ_2p_old(idx) = calc_ρ_2p(var_old.vl.T[idx], var_old.vl.α[idx], FLUID)
    h_1p(idx) = calc_h_1p(var.vl.T[idx], var.vl.P[idx], FLUID)
    h_2p(idx) = calc_h_2p(var.vl.T[idx], var.vl.α[idx], FLUID)
    e_1p(idx) = calc_e_1p(var.vl.T[idx], var.vl.P[idx], FLUID)
    e_2p(idx) = calc_e_2p(var.vl.T[idx], var.vl.α[idx], FLUID)
    e_1p_old(idx) = calc_e_1p(var_old.vl.T[idx], var_old.vl.P[idx], FLUID)
    e_2p_old(idx) = calc_e_2p(var_old.vl.T[idx], var_old.vl.α[idx], FLUID)
    h_wall_f_1p(idx) = calc_h_wall_f_1p(var.vl.T[idx], var.vl.P[idx], var.vl.m[idx], D_i, A_cross, FLUID)
    h_wall_f_2p(idx) = calc_h_wall_f_2p(var.vl.T[idx], var.vl.α[idx], var.vl.m[idx], D_i, A_cross, FLUID)
    τ_1p(idx) = calc_τ_wall(var.vl.T[idx], var.vl.P[idx], var.vl.m[idx], D_i, A_cross, FLUID)
    ΔP_1p(idx) = τ_1p(idx) * S_inner
    ΔP_2p(idx) = calc_ΔP_2p(var.vl.T[idx], var.vl.α[idx], var.vl.m[idx], D_i, A_cross, FLUID)
    E_1p(idx) = V * ρ_1p(idx) * (e_1p(idx) + 0.5 * (var.vl.m[idx] / (ρ_1p(idx) * A_cross))^2)
    E_2p(idx) = V * ρ_2p(idx) * (e_2p(idx) + 0.5 * (var.vl.m[idx] / (ρ_2p(idx) * A_cross))^2)
    E_1p_old(idx) = V * ρ_1p_old(idx) * (e_1p_old(idx) + 0.5 * (var.vl.m[idx] / (ρ_1p_old(idx) * A_cross))^2)
    E_2p_old(idx) = V * ρ_2p_old(idx) * (e_2p_old(idx) + 0.5 * (var.vl.m[idx] / (ρ_2p_old(idx) * A_cross))^2)
    Q_1p(idx) = h_wall_f_1p(idx) * S_inner * (var.vl.T_wall[idx] - var.vl.T[idx])
    Q_2p(idx) = h_wall_f_2p(idx) * S_inner * (var.vl.T_wall[idx] - var.vl.T[idx])

    # 残差計算
    if var_old.vl.α[1] < 1.0 && var_old.vl.α[1] > 0.0
        f_mass[1]   = V * (ρ_2p(1) - ρ_2p_old(1)) - dt * (var.vl.m[1] - var.vl.m[2])
        f_mom[1]    = dl * (var.vl.m[1] - var_old.vl.m[1]) - dt * ( (m_grv^2)/(ρ_2p_grv*A_cross) - ((var.vl.m[1])^2)/(ρ_2p(1)*A_cross) + A_cross * (P_grv - var.vl.P[1]) - dl * ΔP_2p(1))
        f_ene[1]    = E_2p(1) - E_2p_old(1) - dt * (var.vl.m[1] * (h_2p(1) + 0.5 * (var.vl.m[1]/(ρ_2p(1)*A_cross))^2) - var.vl.m[2] * (h_2p(2) + 0.5 * (var.vl.m[2]/(ρ_2p(2)*A_cross))^2) + Q_2p(1))
        f_state[1]  = var.vl.T[1] - calc_Tsat(var.vl.P[1], FLUID)
    else
        f_mass[1]   = V * (ρ_1p(1) - ρ_1p_old(1)) - dt * (var.vl.m[1] - var.vl.m[2])
        f_mom[1]    = dl * (var.vl.m[1] - var_old.vl.m[1]) - dt * ( (m_grv^2)/(ρ_1p_grv*A_cross) - ((var.vl.m[1])^2)/(ρ_1p(1)*A_cross) + A_cross * (P_grv - var.vl.P[1]) - dl * ΔP_1p(1))
        f_ene[1]    = E_1p(1) - E_1p_old(1) - dt * (var.vl.m[1] * (h_1p(1) + 0.5 * (var.vl.m[1]/(ρ_1p(1)*A_cross))^2) - var.vl.m[2] * (h_1p(2) + 0.5 * (var.vl.m[2]/(ρ_1p(2)*A_cross))^2) + Q_1p(1))
        if var_old.vl.α[1] == 1.0
            f_state[1] = var.vl.α[1] - 1.0
        elseif var_old.vl.α[1] == 0.0
            f_state[1] = var.vl.α[1] - 0.0
        end
    end

    for i in 2:geo.vl.N-1
        if var_old.vl.α[i] < 1.0 && var_old.vl.α[i] > 0.0 # 二相流
            f_mass[i]   = V * (ρ_2p(i) - ρ_2p_old(i)) - dt * (var.vl.m[i] - var.vl.m[i+1])
            f_mom[i]    = dl * (var.vl.m[i] - var_old.vl.m[i]) - dt * ( ((var.vl.m[i-1])^2)/(ρ_2p(i-1)*A_cross) - ((var.vl.m[i])^2)/(ρ_2p(i)*A_cross) + A_cross * (var.vl.P[i-1] - var.vl.P[i]) - dl * ΔP_2p(i))
            f_ene[i]    = E_2p(i) - E_2p_old(i) - dt * (var.vl.m[i] * (h_2p(i) + 0.5 * (var.vl.m[i]/(ρ_2p(i)*A_cross))^2) - var.vl.m[i+1] * (h_2p(i+1) + 0.5 * (var.vl.m[i+1]/(ρ_2p(i+1)*A_cross))^2) + Q_2p(i))
            f_state[i]  = var.vl.T[i] - calc_Tsat(var.vl.P[i], FLUID)
        else # 単相流
            f_mass[i]   = V * (ρ_1p(i) - ρ_1p_old(i)) - dt * (var.vl.m[i] - var.vl.m[i+1])
            f_mom[i]    = dl * (var.vl.m[i] - var_old.vl.m[i]) - dt * ( ((var.vl.m[i-1])^2)/(ρ_1p(i-1)*A_cross) - ((var.vl.m[i])^2)/(ρ_1p(i)*A_cross) + A_cross * (var.vl.P[i-1] - var.vl.P[i]) - dl * ΔP_1p(i))
            f_ene[i]    = E_1p(i) - E_1p_old(i) - dt * (var.vl.m[i] * (h_1p(i) + 0.5 * (var.vl.m[i]/(ρ_1p(i)*A_cross))^2) - var.vl.m[i+1] * (h_1p(i+1) + 0.5 * (var.vl.m[i+1]/(ρ_1p(i+1)*A_cross))^2) + Q_1p(i))
            if var_old.vl.α[i] == 1.0
                f_state[i] = var.vl.α[i] - 1.0
            elseif var_old.vl.α[i] == 0.0
                f_state[i] = var.vl.α[i] - 0.0
            end
        end
    end

    if var_old.vl.α[geo.vl.N] < 1.0 && var_old.vl.α[geo.vl.N] > 0.0
        f_mass[geo.vl.N] = V * (ρ_2p(geo.vl.N) - ρ_2p_old(geo.vl.N)) - dt * (var.vl.m[geo.vl.N] - m_cond)
        f_mom[geo.vl.N]  = dl * (var.vl.m[geo.vl.N] - var_old.vl.m[geo.vl.N]) - dt * ( (var.vl.m[geo.vl.N-1]^2)/(ρ_2p(geo.vl.N-1)*A_cross) - ((var.vl.m[geo.vl.N])^2)/(ρ_2p(geo.vl.N)*A_cross) + A_cross * (var.vl.P[geo.vl.N-1] - var.vl.P[geo.vl.N]) - dl * ΔP_2p(geo.vl.N))
        f_ene[geo.vl.N]  = E_2p(geo.vl.N) - E_2p_old(geo.vl.N) - dt * (var.vl.m[geo.vl.N] * (h_2p(geo.vl.N) + 0.5 * (var.vl.m[geo.vl.N]/(ρ_2p(geo.vl.N)*A_cross))^2) - m_cond * (h_2p_cond + 0.5 * (m_cond/(ρ_2p_cond*A_cross))^2) + Q_2p(geo.vl.N))
        f_state[geo.vl.N]  = var.vl.T[geo.vl.N] - calc_Tsat(var.vl.P[geo.vl.N], FLUID)
    else
        f_mass[geo.vl.N] = V * (ρ_1p(geo.vl.N) - ρ_1p_old(geo.vl.N)) - dt * (var.vl.m[geo.vl.N] - m_cond)
        f_mom[geo.vl.N]  = dl * (var.vl.m[geo.vl.N] - var_old.vl.m[geo.vl.N]) - dt * ( (var.vl.m[geo.vl.N-1]^2)/(ρ_1p(geo.vl.N-1)*A_cross) - ((var.vl.m[geo.vl.N])^2)/(ρ_1p(geo.vl.N)*A_cross) + A_cross * (var.vl.P[geo.vl.N-1] - var.vl.P[geo.vl.N]) - dl * ΔP_1p(geo.vl.N))
        f_ene[geo.vl.N]  = E_1p(geo.vl.N) - E_1p_old(geo.vl.N) - dt * (var.vl.m[geo.vl.N] * (h_1p(geo.vl.N) + 0.5 * (var.vl.m[geo.vl.N]/(ρ_1p(geo.vl.N)*A_cross))^2) - m_cond * (h_1p_cond + 0.5 * (m_cond/(ρ_1p_cond*A_cross))^2) + Q_1p(geo.vl.N))
        if var_old.vl.α[geo.vl.N] == 1.0
            f_state[geo.vl.N] = var.vl.α[geo.vl.N] - 1.0
        elseif var_old.vl.α[geo.vl.N] == 0.0
            f_state[geo.vl.N] = var.vl.α[geo.vl.N] - 0.0
        end
    end

end

function calc_transient_cond!(f_mass, f_mom, f_ene, f_state, var, var_old, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt

    # 二相流チェック
    for i in 1:geo.cond.N
        check = TPcheck(var.cond.T[i], var.cond.P[i], FLUID)
        if check == :tp
            var.cond.T[i] = calc_Tsat(var.cond.P[i], FLUID)
            if var.cond.α[i] > 1.0
                var.cond.α[i] = 1.0
            elseif var.cond.α[i] < 0.0
                var.cond.α[i] = 0.0
            end
        elseif check == :vapor
            var.cond.α[i] = 1.0
        elseif check == :liquid
            var.cond.α[i] = 0.0   
        else
            error("invalid phase")
        end
    end

    T_vl   = var.vl.T[geo.vl.N]
    T_ll   = var.ll.T[1]
    m_vl   = var.vl.m[geo.vl.N]
    m_ll   = var.ll.m[1]
    P_vl   = var.vl.P[geo.vl.N]
    P_ll   = var.ll.P[1]
    ρ_1p_vl  = calc_ρ_1p(T_vl, P_vl, FLUID)
    ρ_2p_vl  = calc_ρ_2p(T_vl, var.vl.α[geo.vl.N], FLUID)
    ρ_1p_ll  = calc_ρ_1p(T_ll, P_ll, FLUID)
    ρ_2p_ll = calc_ρ_2p(T_ll, var.ll.α[1], FLUID)
    h_1p_ll = calc_h_1p(T_ll, P_ll, FLUID)
    h_2p_ll = calc_h_2p(T_ll, var.ll.α[1], FLUID)    

    dl = geo.cond.dl
    D_i = geo.cond.D_i
    D_o = geo.cond.D_o
    A_cross = geo.cond.A_cross
    S_inner = geo.cond.S_inner
    V = geo.cond.V

    check::Symbol = :init

    ρ_1p(idx) = calc_ρ_1p(var.cond.T[idx], var.cond.P[idx], FLUID)
    ρ_2p(idx) = calc_ρ_2p(var.cond.T[idx], var.cond.α[idx], FLUID)
    ρ_1p_old(idx) = calc_ρ_1p(var_old.cond.T[idx], var_old.cond.P[idx], FLUID)
    ρ_2p_old(idx) = calc_ρ_2p(var_old.cond.T[idx], var_old.cond.α[idx], FLUID)
    h_1p(idx) = calc_h_1p(var.cond.T[idx], var.cond.P[idx], FLUID)
    h_2p(idx) = calc_h_2p(var.cond.T[idx], var.cond.α[idx], FLUID)
    e_1p(idx) = calc_e_1p(var.cond.T[idx], var.cond.P[idx], FLUID)
    e_2p(idx) = calc_e_2p(var.cond.T[idx], var.cond.α[idx], FLUID)
    e_1p_old(idx) = calc_e_1p(var_old.cond.T[idx], var_old.cond.P[idx], FLUID)
    e_2p_old(idx) = calc_e_2p(var_old.cond.T[idx], var_old.cond.α[idx], FLUID)
    h_wall_f_1p(idx) = calc_h_wall_f_1p(var.cond.T[idx], var.cond.P[idx], var.cond.m[idx], D_i, A_cross, FLUID)
    h_wall_f_2p(idx) = calc_h_wall_f_2p(var.cond.T[idx], var.cond.α[idx], var.cond.m[idx], D_i, A_cross, FLUID)
    τ_1p(idx) = calc_τ_wall(var.cond.T[idx], var.cond.P[idx], var.cond.m[idx], D_i, A_cross, FLUID)
    ΔP_1p(idx) = τ_1p(idx) * S_inner
    ΔP_2p(idx) = calc_ΔP_2p(var.cond.T[idx], var.cond.α[idx], var.cond.m[idx], D_i, A_cross, FLUID)
    E_1p(idx) = V * ρ_1p(idx) * (e_1p(idx) + 0.5 * (var.cond.m[idx] / (ρ_1p(idx) * A_cross))^2)
    E_2p(idx) = V * ρ_2p(idx) * (e_2p(idx) + 0.5 * (var.cond.m[idx] / (ρ_2p(idx) * A_cross))^2)
    E_1p_old(idx) = V * ρ_1p_old(idx) * (e_1p_old(idx) + 0.5 * (var.cond.m[idx] / (ρ_1p_old(idx) * A_cross))^2)
    E_2p_old(idx) = V * ρ_2p_old(idx) * (e_2p_old(idx) + 0.5 * (var.cond.m[idx] / (ρ_2p_old(idx) * A_cross))^2)
    Q_1p(idx) = h_wall_f_1p(idx) * S_inner * (var.cond.T_wall[idx] - var.cond.T[idx])
    Q_2p(idx) = h_wall_f_2p(idx) * S_inner * (var.cond.T_wall[idx] - var.cond.T[idx])

    

    # 残差計算
    if var_old.cond.α[1] < 1.0 && var_old.cond.α[1] > 0.0
        f_mass[1]   = V * (ρ_2p(1) - ρ_2p_old(1)) - dt * (var.cond.m[1] - var.cond.m[2])
        f_mom[1]    = dl * (var.cond.m[1] - var_old.cond.m[1]) - dt * ( (m_vl^2)/(ρ_2p_vl*A_cross) - ((var.cond.m[1])^2)/(ρ_2p(1)*A_cross) + A_cross * (P_vl - var.cond.P[1]) - dl * ΔP_2p(1))
        f_ene[1]    = E_2p(1) - E_2p_old(1) - dt * (var.cond.m[1] * (h_2p(1) + 0.5 * (var.cond.m[1]/(ρ_2p(1)*A_cross))^2) - var.cond.m[2] * (h_2p(2) + 0.5 * (var.cond.m[2]/(ρ_2p(2)*A_cross))^2) + Q_2p(1))
        f_state[1]  = var.cond.T[1] - calc_Tsat(var.cond.P[1], FLUID)
    else
        f_mass[1]   = V * (ρ_1p(1) - ρ_1p_old(1)) - dt * (var.cond.m[1] - var.cond.m[2])
        f_mom[1]    = dl * (var.cond.m[1] - var_old.cond.m[1]) - dt * ( (m_vl^2)/(ρ_1p_vl*A_cross) - ((var.cond.m[1])^2)/(ρ_1p(1)*A_cross) + A_cross * (P_vl - var.cond.P[1]) - dl * ΔP_1p(1))
        f_ene[1]    = E_1p(1) - E_1p_old(1) - dt * (var.cond.m[1] * (h_1p(1) + 0.5 * (var.cond.m[1]/(ρ_1p(1)*A_cross))^2) - var.cond.m[2] * (h_1p(2) + 0.5 * (var.cond.m[2]/(ρ_1p(2)*A_cross))^2) + Q_1p(1))
        if var_old.cond.α[1] == 1.0
            f_state[1] = var.cond.α[1] - 1.0
        elseif var_old.cond.α[1] == 0.0
            f_state[1] = var.cond.α[1] - 0.0
        end
    end

    for i in 2:geo.cond.N-1
        if var_old.cond.α[i] < 1.0 && var_old.cond.α[i] > 0.0 # 二相流
            f_mass[i]   = V * (ρ_2p(i) - ρ_2p_old(i)) - dt * (var.cond.m[i] - var.cond.m[i+1])
            f_mom[i]    = dl * (var.cond.m[i] - var_old.cond.m[i]) - dt * ( ((var.cond.m[i-1])^2)/(ρ_2p(i-1)*A_cross) - ((var.cond.m[i])^2)/(ρ_2p(i)*A_cross) + A_cross * (var.cond.P[i-1] - var.cond.P[i]) - dl * ΔP_2p(i))
            f_ene[i]    = E_2p(i) - E_2p_old(i) - dt * (var.cond.m[i] * (h_2p(i) + 0.5 * (var.cond.m[i]/(ρ_2p(i)*A_cross))^2) - var.cond.m[i+1] * (h_2p(i+1) + 0.5 * (var.cond.m[i+1]/(ρ_2p(i+1)*A_cross))^2) + Q_2p(i))
            f_state[i]  = var.cond.T[i] - calc_Tsat(var.cond.P[i], FLUID)
        else # 単相流
            f_mass[i]   = V * (ρ_1p(i) - ρ_1p_old(i)) - dt * (var.cond.m[i] - var.cond.m[i+1])
            f_mom[i]    = dl * (var.cond.m[i] - var_old.cond.m[i]) - dt * ( ((var.cond.m[i-1])^2)/(ρ_1p(i-1)*A_cross) - ((var.cond.m[i])^2)/(ρ_1p(i)*A_cross) + A_cross * (var.cond.P[i-1] - var.cond.P[i]) - dl * ΔP_1p(i))
            f_ene[i]    = E_1p(i) - E_1p_old(i) - dt * (var.cond.m[i] * (h_1p(i) + 0.5 * (var.cond.m[i]/(ρ_1p(i)*A_cross))^2) - var.cond.m[i+1] * (h_1p(i+1) + 0.5 * (var.cond.m[i+1]/(ρ_1p(i+1)*A_cross))^2) + Q_1p(i))
            if var_old.cond.α[i] == 1.0
                f_state[i] = var.cond.α[i] - 1.0
            elseif var_old.cond.α[i] == 0.0
                f_state[i] = var.cond.α[i] - 0.0
            end
        end
    end

    if var_old.cond.α[geo.cond.N] < 1.0 && var_old.cond.α[geo.cond.N] > 0.0
        f_mass[geo.cond.N] = V * (ρ_2p(geo.cond.N) - ρ_2p_old(geo.cond.N)) - dt * (var.cond.m[geo.cond.N] - m_ll)
        f_mom[geo.cond.N]  = dl * (var.cond.m[geo.cond.N] - var_old.cond.m[geo.cond.N]) - dt * ( (var.cond.m[geo.cond.N-1]^2)/(ρ_2p(geo.cond.N-1)*A_cross) - ((var.cond.m[geo.cond.N])^2)/(ρ_2p(geo.cond.N)*A_cross) + A_cross * (var.cond.P[geo.cond.N-1] - var.cond.P[geo.cond.N]) - dl * ΔP_2p(geo.cond.N))
        f_ene[geo.cond.N]  = E_2p(geo.cond.N) - E_2p_old(geo.cond.N) - dt * (var.cond.m[geo.cond.N] * (h_2p(geo.cond.N) + 0.5 * (var.cond.m[geo.cond.N]/(ρ_2p(geo.cond.N)*A_cross))^2) - m_ll * (h_2p_ll + 0.5 * (m_ll/(ρ_2p_ll*A_cross))^2) + Q_2p(geo.cond.N))
        f_state[geo.cond.N]  = var.cond.T[geo.cond.N] - calc_Tsat(var.cond.P[geo.cond.N], FLUID)
    else
        f_mass[geo.cond.N] = V * (ρ_1p(geo.cond.N) - ρ_1p_old(geo.cond.N)) - dt * (var.cond.m[geo.cond.N] - m_ll)
        f_mom[geo.cond.N]  = dl * (var.cond.m[geo.cond.N] - var_old.cond.m[geo.cond.N]) - dt * ( (var.cond.m[geo.cond.N-1]^2)/(ρ_1p(geo.cond.N-1)*A_cross) - ((var.cond.m[geo.cond.N])^2)/(ρ_1p(geo.cond.N)*A_cross) + A_cross * (var.cond.P[geo.cond.N-1] - var.cond.P[geo.cond.N]) - dl * ΔP_1p(geo.cond.N))
        f_ene[geo.cond.N]  = E_1p(geo.cond.N) - E_1p_old(geo.cond.N) - dt * (var.cond.m[geo.cond.N] * (h_1p(geo.cond.N) + 0.5 * (var.cond.m[geo.cond.N]/(ρ_1p(geo.cond.N)*A_cross))^2) - m_ll * (h_1p_ll + 0.5 * (m_ll/(ρ_1p_ll*A_cross))^2) + Q_1p(geo.cond.N))
        if var_old.cond.α[geo.cond.N] == 1.0
            f_state[geo.cond.N] = var.cond.α[geo.cond.N] - 1.0
        elseif var_old.cond.α[geo.cond.N] == 0.0
            f_state[geo.cond.N] = var.cond.α[geo.cond.N] - 0.0
        end
    end

end

function calc_transient_ll!(f_mass, f_mom, f_ene, f_state, var, var_old, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt

    # 二相流チェック
    for i in 1:geo.ll.N
        check = TPcheck(var.ll.T[i], var.ll.P[i], FLUID)
        if check == :tp
            var.ll.T[i] = calc_Tsat(var.ll.P[i], FLUID)
            if var.ll.α[i] > 1.0
                var.ll.α[i] = 1.0
            elseif var.ll.α[i] < 0.0
                var.ll.α[i] = 0.0
            end
        elseif check == :vapor
            var.ll.α[i] = 1.0
        elseif check == :liquid
            var.ll.α[i] = 0.0   
        else
            error("invalid phase")
        end
    end

    T_ll   = var.ll.T[geo.ll.N]
    T_bay   = var.bay.T[1]
    m_ll   = var.ll.m[geo.ll.N]
    m_bay   = var.bay.m[1]
    P_ll   = var.ll.P[geo.ll.N]
    P_bay   = var.bay.P[1]
    ρ_1p_ll  = calc_ρ_1p(T_ll, P_ll, FLUID)
    ρ_2p_ll  = calc_ρ_2p(T_ll, var.ll.α[geo.ll.N], FLUID)
    ρ_1p_bay  = calc_ρ_1p(T_bay, P_bay, FLUID)
    ρ_2p_bay = calc_ρ_2p(T_bay, var.bay.α[1], FLUID)
    h_1p_bay = calc_h_1p(T_bay, P_bay, FLUID)
    h_2p_bay = calc_h_2p(T_bay, var.bay.α[1], FLUID)    

    dl = geo.ll.dl
    D_i = geo.ll.D_i
    D_o = geo.ll.D_o
    A_cross = geo.ll.A_cross
    S_inner = geo.ll.S_inner
    V = geo.ll.V

    check::Symbol = :init

    ρ_1p(idx) = calc_ρ_1p(var.ll.T[idx], var.ll.P[idx], FLUID)
    ρ_2p(idx) = calc_ρ_2p(var.ll.T[idx], var.ll.α[idx], FLUID)
    ρ_1p_old(idx) = calc_ρ_1p(var_old.ll.T[idx], var_old.ll.P[idx], FLUID)
    ρ_2p_old(idx) = calc_ρ_2p(var_old.ll.T[idx], var_old.ll.α[idx], FLUID)
    h_1p(idx) = calc_h_1p(var.ll.T[idx], var.ll.P[idx], FLUID)
    h_2p(idx) = calc_h_2p(var.ll.T[idx], var.ll.α[idx], FLUID)
    e_1p(idx) = calc_e_1p(var.ll.T[idx], var.ll.P[idx], FLUID)
    e_2p(idx) = calc_e_2p(var.ll.T[idx], var.ll.α[idx], FLUID)
    e_1p_old(idx) = calc_e_1p(var_old.ll.T[idx], var_old.ll.P[idx], FLUID)
    e_2p_old(idx) = calc_e_2p(var_old.ll.T[idx], var_old.ll.α[idx], FLUID)
    h_wall_f_1p(idx) = calc_h_wall_f_1p(var.ll.T[idx], var.ll.P[idx], var.ll.m[idx], D_i, A_cross, FLUID)
    h_wall_f_2p(idx) = calc_h_wall_f_2p(var.ll.T[idx], var.ll.α[idx], var.ll.m[idx], D_i, A_cross, FLUID)
    τ_1p(idx) = calc_τ_wall(var.ll.T[idx], var.ll.P[idx], var.ll.m[idx], D_i, A_cross, FLUID)
    ΔP_1p(idx) = τ_1p(idx) * S_inner
    ΔP_2p(idx) = calc_ΔP_2p(var.ll.T[idx], var.ll.α[idx], var.ll.m[idx], D_i, A_cross, FLUID)
    E_1p(idx) = V * ρ_1p(idx) * (e_1p(idx) + 0.5 * (var.ll.m[idx] / (ρ_1p(idx) * A_cross))^2)
    E_2p(idx) = V * ρ_2p(idx) * (e_2p(idx) + 0.5 * (var.ll.m[idx] / (ρ_2p(idx) * A_cross))^2)
    E_1p_old(idx) = V * ρ_1p_old(idx) * (e_1p_old(idx) + 0.5 * (var.ll.m[idx] / (ρ_1p_old(idx) * A_cross))^2)
    E_2p_old(idx) = V * ρ_2p_old(idx) * (e_2p_old(idx) + 0.5 * (var.ll.m[idx] / (ρ_2p_old(idx) * A_cross))^2)
    Q_1p(idx) = h_wall_f_1p(idx) * S_inner * (var.ll.T_wall[idx] - var.ll.T[idx])
    Q_2p(idx) = h_wall_f_2p(idx) * S_inner * (var.ll.T_wall[idx] - var.ll.T[idx])

    

    # 残差計算
    if var_old.ll.α[1] < 1.0 && var_old.ll.α[1] > 0.0
        f_mass[1]   = V * (ρ_2p(1) - ρ_2p_old(1)) - dt * (var.ll.m[1] - var.ll.m[2])
        f_mom[1]    = dl * (var.ll.m[1] - var_old.ll.m[1]) - dt * ( (m_ll^2)/(ρ_2p_ll*A_cross) - ((var.ll.m[1])^2)/(ρ_2p(1)*A_cross) + A_cross * (P_ll - var.ll.P[1]) - dl * ΔP_2p(1))
        f_ene[1]    = E_2p(1) - E_2p_old(1) - dt * (var.ll.m[1] * (h_2p(1) + 0.5 * (var.ll.m[1]/(ρ_2p(1)*A_cross))^2) - var.ll.m[2] * (h_2p(2) + 0.5 * (var.ll.m[2]/(ρ_2p(2)*A_cross))^2) + Q_2p(1))
        f_state[1]  = var.ll.T[1] - calc_Tsat(var.ll.P[1], FLUID)
    else
        f_mass[1]   = V * (ρ_1p(1) - ρ_1p_old(1)) - dt * (var.ll.m[1] - var.ll.m[2])
        f_mom[1]    = dl * (var.ll.m[1] - var_old.ll.m[1]) - dt * ( (m_ll^2)/(ρ_1p_ll*A_cross) - ((var.ll.m[1])^2)/(ρ_1p(1)*A_cross) + A_cross * (P_ll - var.ll.P[1]) - dl * ΔP_1p(1))
        f_ene[1]    = E_1p(1) - E_1p_old(1) - dt * (var.ll.m[1] * (h_1p(1) + 0.5 * (var.ll.m[1]/(ρ_1p(1)*A_cross))^2) - var.ll.m[2] * (h_1p(2) + 0.5 * (var.ll.m[2]/(ρ_1p(2)*A_cross))^2) + Q_1p(1))
        if var_old.ll.α[1] == 1.0
            f_state[1] = var.ll.α[1] - 1.0
        elseif var_old.ll.α[1] == 0.0
            f_state[1] = var.ll.α[1] - 0.0
        end
    end

    for i in 2:geo.ll.N-1
        if var_old.ll.α[i] < 1.0 && var_old.ll.α[i] > 0.0 # 二相流
            f_mass[i]   = V * (ρ_2p(i) - ρ_2p_old(i)) - dt * (var.ll.m[i] - var.ll.m[i+1])
            f_mom[i]    = dl * (var.ll.m[i] - var_old.ll.m[i]) - dt * ( ((var.ll.m[i-1])^2)/(ρ_2p(i-1)*A_cross) - ((var.ll.m[i])^2)/(ρ_2p(i)*A_cross) + A_cross * (var.ll.P[i-1] - var.ll.P[i]) - dl * ΔP_2p(i))
            f_ene[i]    = E_2p(i) - E_2p_old(i) - dt * (var.ll.m[i] * (h_2p(i) + 0.5 * (var.ll.m[i]/(ρ_2p(i)*A_cross))^2) - var.ll.m[i+1] * (h_2p(i+1) + 0.5 * (var.ll.m[i+1]/(ρ_2p(i+1)*A_cross))^2) + Q_2p(i))
            f_state[i]  = var.ll.T[i] - calc_Tsat(var.ll.P[i], FLUID)
        else # 単相流
            f_mass[i]   = V * (ρ_1p(i) - ρ_1p_old(i)) - dt * (var.ll.m[i] - var.ll.m[i+1])
            f_mom[i]    = dl * (var.ll.m[i] - var_old.ll.m[i]) - dt * ( ((var.ll.m[i-1])^2)/(ρ_1p(i-1)*A_cross) - ((var.ll.m[i])^2)/(ρ_1p(i)*A_cross) + A_cross * (var.ll.P[i-1] - var.ll.P[i]) - dl * ΔP_1p(i))
            f_ene[i]    = E_1p(i) - E_1p_old(i) - dt * (var.ll.m[i] * (h_1p(i) + 0.5 * (var.ll.m[i]/(ρ_1p(i)*A_cross))^2) - var.ll.m[i+1] * (h_1p(i+1) + 0.5 * (var.ll.m[i+1]/(ρ_1p(i+1)*A_cross))^2) + Q_1p(i))
            if var_old.ll.α[i] == 1.0
                f_state[i] = var.ll.α[i] - 1.0
            elseif var_old.ll.α[i] == 0.0
                f_state[i] = var.ll.α[i] - 0.0
            end
        end
    end

    if var_old.ll.α[geo.ll.N] < 1.0 && var_old.ll.α[geo.ll.N] > 0.0
        f_mass[geo.ll.N] = V * (ρ_2p(geo.ll.N) - ρ_2p_old(geo.ll.N)) - dt * (var.ll.m[geo.ll.N] - m_bay)
        f_mom[geo.ll.N]  = dl * (var.ll.m[geo.ll.N] - var_old.ll.m[geo.ll.N]) - dt * ( (var.ll.m[geo.ll.N-1]^2)/(ρ_2p(geo.ll.N-1)*A_cross) - ((var.ll.m[geo.ll.N])^2)/(ρ_2p(geo.ll.N)*A_cross) + A_cross * (var.ll.P[geo.ll.N-1] - var.ll.P[geo.ll.N]) - dl * ΔP_2p(geo.ll.N))
        f_ene[geo.ll.N]  = E_2p(geo.ll.N) - E_2p_old(geo.ll.N) - dt * (var.ll.m[geo.ll.N] * (h_2p(geo.ll.N) + 0.5 * (var.ll.m[geo.ll.N]/(ρ_2p(geo.ll.N)*A_cross))^2) - m_bay * (h_2p_bay + 0.5 * (m_bay/(ρ_2p_bay*A_cross))^2) + Q_2p(geo.ll.N))
        f_state[geo.ll.N]  = var.ll.T[geo.ll.N] - calc_Tsat(var.ll.P[geo.ll.N], FLUID)
    else
        f_mass[geo.ll.N] = V * (ρ_1p(geo.ll.N) - ρ_1p_old(geo.ll.N)) - dt * (var.ll.m[geo.ll.N] - m_bay)
        f_mom[geo.ll.N]  = dl * (var.ll.m[geo.ll.N] - var_old.ll.m[geo.ll.N]) - dt * ( (var.ll.m[geo.ll.N-1]^2)/(ρ_1p(geo.ll.N-1)*A_cross) - ((var.ll.m[geo.ll.N])^2)/(ρ_1p(geo.ll.N)*A_cross) + A_cross * (var.ll.P[geo.ll.N-1] - var.ll.P[geo.ll.N]) - dl * ΔP_1p(geo.ll.N))
        f_ene[geo.ll.N]  = E_1p(geo.ll.N) - E_1p_old(geo.ll.N) - dt * (var.ll.m[geo.ll.N] * (h_1p(geo.ll.N) + 0.5 * (var.ll.m[geo.ll.N]/(ρ_1p(geo.ll.N)*A_cross))^2) - m_bay * (h_1p_bay + 0.5 * (m_bay/(ρ_1p_bay*A_cross))^2) + Q_1p(geo.ll.N))
        if var_old.ll.α[geo.ll.N] == 1.0
            f_state[geo.ll.N] = var.ll.α[geo.ll.N] - 1.0
        elseif var_old.ll.α[geo.ll.N] == 0.0
            f_state[geo.ll.N] = var.ll.α[geo.ll.N] - 0.0
        end
    end

end

function pipe_eq_evap(var_old, geo, condition)
    FLUID = condition.fluid
    Q_apply = condition.Q_apply / geo.grv.N # 各ノードに均等に入熱を分割
    h_ev = condition.h_ev

    # Q_ev
    A_ec_wick = geo.grv.A_ec_wick
    T_ec = var_old.grv.T_wall[1]
    T_ev = var_old.wick.T_in
    Q_ev = h_ev * A_ec_wick * (T_ec - T_ev) / geo.grv.N # 各ノードに均等に入熱を分割

    # Q_gr
    n_grv = geo.grv.n_grv
    dl = geo.grv.dl
    width = geo.grv.width
    height = geo.grv.height
    A_grv_cross = geo.grv.A_grv_cross
    h_ec_gr(idx) = calc_h_ec_gr(var_old.grv.T[idx], var_old.grv.P[idx], var_old.grv.m[idx], height, A_grv_cross, geo.grv.L, FLUID)
    Q_gr(idx) = h_ec_gr(idx) * (width * dl) * (var_old.grv.T_wall[idx] - var_old.grv.T[idx]) * n_grv

    # Q_leak
    k_eff = calc_k_eff(T_ev, FLUID)
    L_wick = geo.grv.L
    A_ec_cc = (π / 4) * (geo.cc.D_o^2 - geo.cc.D_i^2)
    l_ec_cc = (geo.grv.L + geo.cc.L) / 2
    G_ec_cc = k_sus304(T_ec) * (A_ec_cc / l_ec_cc)
    T_cc = var_old.cc.T_cc
    
    A_ec_vl = (π / 4) * (geo.grv.D_o^2 - geo.grv.D_i^2)
    l_ec_vl = (geo.grv.L + geo.vl.dl) / 2
    G_ec_vl = k_sus304(T_ec) * (A_ec_vl / l_ec_vl)

    Q_cc_wick = ((2 * π * k_eff * L_wick) * (T_ev - T_cc) / (log(geo.grv.D_o / geo.grv.D_i)))/geo.grv.N # 各ノードに均等に入熱を分割
    Q_cc_pipe = G_ec_cc * (T_ec - T_cc)
    Q_vl = G_ec_vl * (T_ec - var_old.vl.T_wall[1])

    # C
    dT_wall = zeros(geo.grv.N)
    V_ec = π * (geo.grv.D_o^2 - geo.grv.D_i^2) * geo.grv.dl
    ρ_ec(idx) = ρ_sus304(var_old.grv.T_wall[idx])
    c_ec(idx) = c_sus304(var_old.grv.T_wall[idx])
    C(idx)    = V_ec * ρ_ec(idx) * c_ec(idx)

    dT_wall[1] = (Q_apply - Q_ev - Q_gr(1) - Q_cc_wick - Q_cc_pipe) / C(1)
    for i in 2:geo.grv.N-1
        dT_wall[i] = (Q_apply - Q_ev - Q_gr(i) - Q_cc_wick) / C(i)
    end
    dT_wall[geo.grv.N] = (Q_apply - Q_ev - Q_gr(geo.grv.N) - Q_cc_wick - Q_vl) / C(geo.grv.N)
    return dT_wall
end

function pipe_eq_vl(var_old, geo, condition)
    FLUID = condition.fluid

    D_i = geo.vl.D_i
    A_cross = geo.vl.A_cross
    S_inner = geo.vl.S_inner
    A = (π / 4) * (geo.vl.D_o^2 - geo.vl.D_i^2)
    l = geo.vl.dl
    G(idx) = k_sus304(var_old.vl.T_wall[idx]) * (A / l)

    h_wall_f_1p(idx) = calc_h_wall_f_1p(var_old.vl.T[idx], var_old.vl.P[idx], var_old.vl.m[idx], D_i, A_cross, FLUID)
    h_wall_f_2p(idx) = calc_h_wall_f_2p(var_old.vl.T[idx], var_old.vl.α[idx], var_old.vl.m[idx], D_i, A_cross, FLUID)

    dT_wall = zeros(geo.vl.N)
    V_ec = π * (geo.vl.D_o^2 - geo.vl.D_i^2) * geo.vl.dl
    ρ_ec(idx) = ρ_sus304(var_old.vl.T_wall[idx])
    c_ec(idx) = c_sus304(var_old.vl.T_wall[idx])
    C(idx)    = V_ec * ρ_ec(idx) * c_ec(idx)

    for i in 1:geo.vl.N
        Q_cond = 0
        if i == 1;          Q_cond += G(i) * (var_old.grv.T_wall[geo.grv.N] - var_old.vl.T_wall[i]); end
        if i > 1;           Q_cond += G(i) * (var_old.vl.T_wall[i-1] - var_old.vl.T_wall[i]); end
        if i < geo.vl.N;    Q_cond += G(i) * (var_old.vl.T_wall[i+1] - var_old.vl.T_wall[i]); end
        if i == geo.vl.N;   Q_cond += G(i) * (var_old.cond.T_wall[1] - var_old.vl.T_wall[i]); end

        if var_old.vl.α[i] < 1.0 && var_old.vl.α[i] > 0.0 # 二相流
            Q_conv = h_wall_f_2p(i) * S_inner * (var_old.vl.T[i] - var_old.vl.T_wall[i])
        else
            Q_conv = h_wall_f_1p(i) * S_inner * (var_old.vl.T[i] - var_old.vl.T_wall[i])
        end

        dT_wall[i] = (Q_cond + Q_conv) / C(i)
    end        
    return dT_wall
end

function pipe_eq_cond(var_old, geo, condition) # 使わないかもしれない
    FLUID = condition.fluid
    T_cond = condition.T_cond

    D_i = geo.cond.D_i
    A_cross = geo.cond.A_cross
    S_inner = geo.cond.S_inner
    A = (π / 4) * (geo.cond.D_o^2 - geo.cond.D_i^2)
    l = geo.cond.dl
    G(idx) = k_sus304(var_old.cond.T_wall[idx]) * (A / l)

    h_wall_f_1p(idx) = calc_h_wall_f_1p(var_old.cond.T[idx], var_old.cond.P[idx], var_old.cond.m[idx], D_i, A_cross, FLUID)
    h_wall_f_2p(idx) = calc_h_wall_f_2p(var_old.cond.T[idx], var_old.cond.α[idx], var_old.cond.m[idx], D_i, A_cross, FLUID)

    dT_wall = zeros(geo.cond.N)
    V_ec = π * (geo.cond.D_o^2 - geo.cond.D_i^2) * geo.cond.dl
    ρ_ec(idx) = ρ_sus304(var_old.cond.T_wall[idx])
    c_ec(idx) = c_sus304(var_old.cond.T_wall[idx])
    C(idx)    = V_ec * ρ_ec(idx) * c_ec(idx)

    for i in 1:geo.cond.N
        Q_cond = 0
        if i == 1;          Q_cond += G(i) * (var_old.vl.T_wall[geo.vl.N] - var_old.cond.T_wall[i]); end
        if i > 1;           Q_cond += G(i) * (var_old.cond.T_wall[i-1] - var_old.cond.T_wall[i]); end
        if i < geo.cond.N;    Q_cond += G(i) * (var_old.cond.T_wall[i+1] - var_old.cond.T_wall[i]); end
        if i == geo.cond.N;   Q_cond += G(i) * (var_old.ll.T_wall[1] - var_old.cond.T_wall[i]); end

        if var_old.cond.α[i] < 1.0 && var_old.cond.α[i] > 0.0 # 二相流
            Q_conv = h_wall_f_2p(i) * S_inner * (var_old.cond.T[i] - var_old.cond.T_wall[i])
        else
            Q_conv = h_wall_f_1p(i) * S_inner * (var_old.cond.T[i] - var_old.cond.T_wall[i])
        end

        dT_wall[i] = (Q_cond + Q_conv) / C(i)
    end
    return dT_wall        
end

function pipe_eq_ll(var_old, geo, condition)
    FLUID = condition.fluid

    D_i = geo.ll.D_i
    A_cross = geo.ll.A_cross
    S_inner = geo.ll.S_inner
    A = (π / 4) * (geo.ll.D_o^2 - geo.ll.D_i^2)
    l = geo.ll.dl
    G(idx) = k_sus304(var_old.ll.T_wall[idx]) * (A / l)

    h_wall_f_1p(idx) = calc_h_wall_f_1p(var_old.ll.T[idx], var_old.ll.P[idx], var_old.ll.m[idx], D_i, A_cross, FLUID)
    h_wall_f_2p(idx) = calc_h_wall_f_2p(var_old.ll.T[idx], var_old.ll.α[idx], var_old.ll.m[idx], D_i, A_cross, FLUID)

    dT_wall = zeros(geo.ll.N)
    V_ec = π * (geo.ll.D_o^2 - geo.ll.D_i^2) * geo.ll.dl
    ρ_ec(idx) = ρ_sus304(var_old.ll.T_wall[idx])
    c_ec(idx) = c_sus304(var_old.ll.T_wall[idx])
    C(idx)    = V_ec * ρ_ec(idx) * c_ec(idx)

    for i in 1:geo.ll.N
        Q_cond = 0
        if i == 1;          Q_cond += G(i) * (var_old.cond.T_wall[geo.cond.N] - var_old.ll.T_wall[i]); end
        if i > 1;           Q_cond += G(i) * (var_old.ll.T_wall[i-1] - var_old.ll.T_wall[i]); end
        if i < geo.ll.N;    Q_cond += G(i) * (var_old.ll.T_wall[i+1] - var_old.ll.T_wall[i]); end
        if i == geo.ll.N;   Q_cond += G(i) * (var_old.cc.T_wall[1] - var_old.ll.T_wall[i]); end

        if var_old.ll.α[i] < 1.0 && var_old.ll.α[i] > 0.0 # 二相流
            Q_conv = h_wall_f_2p(i) * S_inner * (var_old.ll.T[i] - var_old.ll.T_wall[i])
        else
            Q_conv = h_wall_f_1p(i) * S_inner * (var_old.ll.T[i] - var_old.ll.T_wall[i])
        end

        dT_wall[i] = (Q_cond + Q_conv) / C(i)
    end        
    return dT_wall
end

function pipe_eq_cc(var_old, geo, condition)
    FLUID = condition.fluid

    D_i = geo.cc.D_i
    A_cross = geo.cc.A_cross
    S_inner = geo.cc.S_inner
    A = (π / 4) * (geo.cc.D_o^2 - geo.cc.D_i^2)
    l = geo.cc.dl
    G = k_sus304(var_old.cc.T_wall) * (A / l)

    # Q_leak
    T_ec = var_old.grv.T_wall[1]
    T_ev = var_old.wick.T_in
    T_cc = var_old.cc.T_wall
    k_eff = calc_k_eff(T_ev, FLUID)
    L_wick = geo.grv.L
    A_ec_cc = (π / 4) * (geo.cc.D_o^2 - geo.cc.D_i^2)
    l_ec_cc = (geo.grv.L + geo.cc.L) / 2
    G_ec_cc = k_sus304(T_ec) * (A_ec_cc / l_ec_cc)
    Q_leak = G_ec_cc * (T_ec - T_cc)

    # Q_conv
    G_conv = calc_k_2p(var_old.cc.T_cc, var_old.cc.α_cc, FLUID) * S_inner / (geo.cc.D_o - geo.cc.D_i)

    dT_wall = zeros(geo.cc.N)
    V_ec = π * (geo.cc.D_o^2 - geo.cc.D_i^2) * geo.cc.L
    ρ_ec = ρ_sus304(var_old.cc.T_wall)
    c_ec = c_sus304(var_old.cc.T_wall)
    C    = V_ec * ρ_ec * c_ec


    Q_cond = 0
    Q_cond += G * (var_old.ll.T_wall[geo.ll.N] - var_old.cc.T_wall)
    Q_cond += G * (var_old.grv.T_wall[1] - var_old.cc.T_wall)

    Q_conv = G_conv * (var_old.cc.T_cc - var_old.cc.T_wall)

    dT_wall = (Q_cond + Q_conv) / C
    
    return dT_wall
end

function Euler_wall(vars_old, T_wall, T_wall_old, geo, condition)
    dt = condition.dt
    T_cond = condition.T_cond

    # 1. 既存の pipe_eq_~ が期待する【T_wall を含んだ古い状態構造】をその場で復元
    full_var_old = ComponentArray(
        wick = vars_old.wick,
        grv  = (T=vars_old.grv.T,   P=vars_old.grv.P,   α=vars_old.grv.α,   m=vars_old.grv.m,   T_wall=T_wall_old.grv),
        vl   = (T=vars_old.vl.T,    P=vars_old.vl.P,    α=vars_old.vl.α,    m=vars_old.vl.m,    T_wall=T_wall_old.vl),
        cond = (T=vars_old.cond.T,  P=vars_old.cond.P,  α=vars_old.cond.α,  m=vars_old.cond.m,  T_wall=T_wall_old.cond),
        ll   = (T=vars_old.ll.T,    P=vars_old.ll.P,    α=vars_old.ll.α,    m=vars_old.ll.m,    T_wall=T_wall_old.ll),
        bay  = (T=vars_old.bay.T,   P=vars_old.bay.P,   α=vars_old.bay.α,   m=vars_old.bay.m,   T_wall=T_wall_old.bay),
        cc   = (T_cc=vars_old.cc.T_cc[1], α_cc=vars_old.cc.α_cc[1], T_wall=T_wall_old.cc[1])
    )

    # 2. 既存の物理計算関数をそのまま呼び出す（中身の変更は一切不要！）
    k1_evap = pipe_eq_evap(full_var_old, geo, condition)
    k1_vl   = pipe_eq_vl(full_var_old, geo, condition)
    k1_cond = pipe_eq_cond(full_var_old, geo, condition) # 今回は固定値にするため使用しませんが、引数は適合させておきます
    k1_ll   = pipe_eq_ll(full_var_old, geo, condition)
    k1_cc   = pipe_eq_cc(full_var_old, geo, condition)

    # 3. 外部管理している T_wall の配列に対してオイラー法を適用し更新
    T_wall.grv  .= T_wall_old.grv  .+ dt .* k1_evap
    T_wall.vl   .= T_wall_old.vl   .+ dt .* k1_vl
    T_wall.cond .= T_cond                               # 凝縮器は指定の固定値
    T_wall.ll   .= T_wall_old.ll   .+ dt .* k1_ll
    T_wall.cc   = T_wall_old.cc   + dt * k1_cc
end

function calc_transient_residual!(f_res, var, var_old, geo, condition)
    calc_transient_wick!(f_res.wick, var, var_old, geo, condition)
    calc_transient_cc!(f_res.cc, var, var_old, geo, condition)
    calc_transient_bay!(f_res.bay.α, f_res.bay.m, f_res.bay.T, f_res.bay.P, var, var_old, geo, condition)
    calc_transient_grv!(f_res.grv.α, f_res.grv.m, f_res.grv.T, f_res.grv.P, var, var_old, geo, condition)
    calc_transient_vl!(f_res.vl.α, f_res.vl.m, f_res.vl.T, f_res.vl.P, var, var_old, geo, condition)
    calc_transient_cond!(f_res.cond.α, f_res.cond.T, f_res.cond.m, f_res.cond.P, var, var_old, geo, condition)
    calc_transient_ll!(f_res.ll.α, f_res.ll.T, f_res.ll.m, f_res.ll.P, var, var_old, geo, condition)
    @show (var.vl.T[1], var.vl.P[1], var.vl.α[1], var.vl.m[1], f_res.vl.T[1])
    @show (var.vl.T[end], var.vl.P[end], var.vl.α[end], var.vl.m[end], f_res.vl.T[end])    
end

function f_nonlinear_transient_solve!(f_res, vars, p)
    f_res .= 0
    # @show vars.grv.P
    # 既存の残差計算関数（calc_transient_residual!）が期待する【T_wall入りの構造】を復元
    full_vars = ComponentArray(
        wick = vars.wick,
        grv  = (T=vars.grv.T,   P=vars.grv.P,   α=vars.grv.α,   m=vars.grv.m,   T_wall=p.T_wall.grv),
        vl   = (T=vars.vl.T,    P=vars.vl.P,    α=vars.vl.α,    m=vars.vl.m,    T_wall=p.T_wall.vl),
        cond = (T=vars.cond.T,  P=vars.cond.P,  α=vars.cond.α,  m=vars.cond.m,  T_wall=p.T_wall.cond),
        ll   = (T=vars.ll.T,    P=vars.ll.P,    α=vars.ll.α,    m=vars.ll.m,    T_wall=p.T_wall.ll),
        bay  = (T=vars.bay.T,   P=vars.bay.P,   α=vars.bay.α,   m=vars.bay.m,   T_wall=p.T_wall.bay),
        cc   = (T_cc=vars.cc.T_cc[1], α_cc=vars.cc.α_cc[1], T_wall=p.T_wall.cc)
    )
    
    full_vars_old = ComponentArray(
        wick = p.vars_old.wick,
        grv  = (T=p.vars_old.grv.T,   P=p.vars_old.grv.P,   α=p.vars_old.grv.α,   m=p.vars_old.grv.m,   T_wall=p.T_wall_old.grv),
        vl   = (T=p.vars_old.vl.T,    P=p.vars_old.vl.P,    α=p.vars_old.vl.α,    m=p.vars_old.vl.m,    T_wall=p.T_wall_old.vl),
        cond = (T=p.vars_old.cond.T,  P=p.vars_old.cond.P,  α=p.vars_old.cond.α,  m=p.vars_old.cond.m,  T_wall=p.T_wall_old.cond),
        ll   = (T=p.vars_old.ll.T,    P=p.vars_old.ll.P,    α=p.vars_old.ll.α,    m=p.vars_old.ll.m,    T_wall=p.T_wall_old.ll),
        bay  = (T=p.vars_old.bay.T,   P=p.vars_old.bay.P,   α=p.vars_old.bay.α,   m=p.vars_old.bay.m,   T_wall=p.T_wall_old.bay),
        cc   = (T_cc=p.vars_old.cc.T_cc[1], α_cc=p.vars_old.cc.α_cc[1], T_wall=p.T_wall_old.cc)
    )

    calc_transient_residual!(f_res, full_vars, full_vars_old, p.geo, p.transient_condition)
end

function solve_transient_lhp!(geo, condition)
    t_init = condition.t_init
    t_end = condition.t_end
    dt = condition.dt

    full_vars_init = load_initial_vars_from_csv(geo, condition)
    vars = ComponentArray(
        wick = full_vars_init.wick, # (T_in, m_wick)
        grv  = (T=copy(full_vars_init.grv.T),   P=copy(full_vars_init.grv.P),   α=copy(full_vars_init.grv.α),   m=copy(full_vars_init.grv.m)),
        vl   = (T=copy(full_vars_init.vl.T),    P=copy(full_vars_init.vl.P),    α=copy(full_vars_init.vl.α),    m=copy(full_vars_init.vl.m)),
        cond = (T=copy(full_vars_init.cond.T),  P=copy(full_vars_init.cond.P),  α=copy(full_vars_init.cond.α),  m=copy(full_vars_init.cond.m)),
        ll   = (T=copy(full_vars_init.ll.T),    P=copy(full_vars_init.ll.P),    α=copy(full_vars_init.ll.α),    m=copy(full_vars_init.ll.m)),
        bay  = (T=copy(full_vars_init.bay.T),   P=copy(full_vars_init.bay.P),   α=copy(full_vars_init.bay.α),   m=copy(full_vars_init.bay.m)),
        cc   = (T_cc=copy(full_vars_init.cc.T_cc), α_cc=copy(full_vars_init.cc.α_cc))
    )
    
    vars_old = copy(vars)

    T_wall = ComponentArray(
        grv  = copy(full_vars_init.grv.T_wall),
        vl   = copy(full_vars_init.vl.T_wall),
        cond = copy(full_vars_init.cond.T_wall),
        ll   = copy(full_vars_init.ll.T_wall),
        bay  = copy(full_vars_init.bay.T_wall),
        cc   = copy(full_vars_init.cc.T_wall)
    )
    T_wall_old = copy(T_wall)

    lb = copy(vars)
    ub = copy(vars)
    for comp in (:grv, :vl, :cond, :ll, :bay)
        # lb.grv や lb.vl を直接参照として取得
        l_obj = getproperty(lb, comp)
        u_obj = getproperty(ub, comp)
        
        # 各変数の上限・下限を一括設定
        l_obj.T .= 200.0;  u_obj.T .= 400.0
        l_obj.P .= 100.0;  u_obj.P .= 1.0e10
        l_obj.α .= -1.0;    u_obj.α .= 2.0
        l_obj.m .= 0.0;  u_obj.m .= 10.0
    end

    lb.wick.T_in   = 200.0
    ub.wick.T_in   = 400.0
    lb.wick.m_wick = 0.0
    ub.wick.m_wick = 10.0
    lb.cc.T_cc     = 200.0
    ub.cc.T_cc     = 400.0
    lb.cc.α_cc     = -1.0
    ub.cc.α_cc     = 2.0


    params = (
        vars_old = vars_old, 
        T_wall = T_wall, 
        T_wall_old = T_wall_old, 
        geo = geo, 
        transient_condition = condition
    )

    for t in t_init:dt:t_end
        T_wall_old .= T_wall
        Euler_wall(params.vars_old, params.T_wall, params.T_wall_old, geo, condition)

        current_params = (
            vars_old = vars_old, # vars_old は vars とは独立したメモリ
            T_wall = T_wall, 
            T_wall_old = T_wall_old, 
            geo = geo, 
            transient_condition = condition
        )

        prob = NonlinearProblem(f_nonlinear_transient_solve!, vars, params; lb=lb, ub=ub)
        println("problem setting completed")
        sol = solve(prob, TrustRegion(autodiff=AutoFiniteDiff()); show_trace = Val(true), abstol=1.0e-5)

        vars .= sol.u
        vars_old .= vars
        save_transient_vars_to_csv(vars, params.T_wall, geo, condition, t)
    end
end