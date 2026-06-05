function calc_transient_wick!(f, var, var_old, geo, condition)
    FLUID = condition.fluid
    h_ev = condition.h_ev
    dt = condition.dt
    Φ = 0.56
    
    T_bay_n = var.bay.T[geo.bay.N]
    T_in = var.wick.T[1]
    T_in_old = var_old.wick.T[1]
    T_out = calc_Tsat(var.grv.P[1], FLUID)
    T_out_old = calc_Tsat(var_old.grv.P[1], FLUID)
    T_ec = var.grv.T_wall[1]

    A_ec_wick = geo.grv.A_ec_wick
    A_mid = geo.wick.A_mid
    A_in = geo.wick.A_in
    V_in = geo.wick.V_in
    V_out = geo.wick.V_out
    dr = (geo.wick.D_3 - geo.wick.D_1) / 2
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

    # ウィック内の質量流量を求める．
    m_wick = (h_ev * A_ec_wick * (T_ec - T_out) - G_wick * (T_out - T_in) - (E_out - E_out_old) / dt) / (λ + h_wick_out - h_wick_in)
    for i in 1:geo.wick.N
        var.wick.m[i] = m_wick
    end

    f = E_in - E_in_old - dt * (m_wick * h_bay_n - m_wick * h_wick_out + G_wick * (T_out - T_in) + G_core_wick * (T_bay_n - T_in))

end
