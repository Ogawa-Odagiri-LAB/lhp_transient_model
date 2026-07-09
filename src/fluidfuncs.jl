function fluid_wick(var, var_old, T_wall, geo, condition) # 返り値 m_evap
    FLUID = condition.fluid
    dt = condition.dt
    ΔT_startup = condition.ΔT_startup
    N = geo.wick.N
    D_in = geo.wick.D_i
    dD = geo.wick.dD

    # m_evap
    T_ev = var.wick[end]
    T_cc = var.cc.T[1]
    T_adj = var.wick[N-1]
    k_eff = calc_k_eff(T_ev, FLUID)
    l_w = geo.wick.L
    D_o = geo.wick.D_o
    D_i = geo.wick.D_o - geo.wick.dD
    Q_ev = calc_Q_ev(var, T_wall, geo, condition)
    Q_in = calc_Q_wick(k_eff, l_w, D_o, D_i, T_adj, T_ev)
    λ = calc_latent_heat(T_ev, FLUID)
    ΔT_ev_cc = abs(T_ev - T_cc)
    m_evap = (ΔT_ev_cc > ΔT_startup) ? (Q_ev + Q_in) / λ : 0.0 # startup

    # residual
    Q_1 = 0.0
    Q_2 = 0.0
    T_1 = 0.0
    T_2 = 0.0
    res = zeros(N)
    for i in 1:N
        T = var.wick[i]
        T_old = var_old.wick[i]
        k_eff = calc_k_eff(T, FLUID)
        D_o = D_in + i * dD
        D_i = D_in + (i - 1) * dD
        
        T_1 = (i == 1) ? var.cc.T[end] : var.wick[i-1]
        T_2 = (i == N) ? var.grv.T[1] : var.wick[i+1] 
        Q_1 = (i == 1) ? 0.0 : calc_Q_wick(k_eff, l_w, D_o, D_i, T_1, T) # 内側との熱伝導による熱授受（最内層では0）
        Q_2 = (i == N) ? Q_ev : calc_Q_wick(k_eff, l_w, D_o, D_i, T_2, T) # 外側との熱伝導による熱授受（最外層ではQ_ev）
        α_1 = (i == 1) ? var.cc.α[end] : 0.0
        α_2 = (i == N) ? var.grv.α[1] : 0.0
        h_1 = calc_h_2p(T_1, α_1, FLUID)
        h_2 = calc_h_2p(T_2, α_2, FLUID)
        if i == N
            V = (geo.wick.width * geo.wick.height + π * (D_o^2 - D_i^2) * 0.25) * geo.wick.L
        else
            V = π * (D_o^2 - D_i^2) * 0.25 * geo.wick.L
        end
        E = calc_wick_energy(T, V, FLUID)
        E_old = calc_wick_energy(T_old, V, FLUID)

        res[i] = (E - E_old) - dt * (m_evap * (h_1 - h_2) + Q_1 + Q_2)
    end
    return m_evap, res
end

function fluid_grv(res_mass, res_mom, res_ene, res_state, m_evap, var, var_old, T_wall, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt
    
    dl = geo.grv.dl
    V = geo.grv.A_ec_wick * geo.grv.height
    A_cross = geo.grv.A_grv_cross * geo.grv.n_grv
    D = geo.grv.height # 代表長さとして

    T = var.grv.T[1]
    T_old = var_old.grv.T[1]
    P = var.grv.P[1]
    α = α_cropper(var.grv.α[1])
    α_old = α_cropper(var_old.grv.α[1])
    m = var.grv.m[1]
    m_old = var_old.grv.m[1]

    T_1 = var.wick[end]
    P_1 = calc_Psat(T_1, FLUID)
    α_1 = 0.0
    m_1 = m_evap

    T_2 = var.vl.T[1]
    α_2 = α_cropper(var.vl.α[1])
    m_2 = var.vl.m[1]

    # Q_grv
    Q_grv = calc_Q_grv(m_evap, var, T_wall, geo, condition)

    # residuals
    res_mass  .= f_mass(V, T, T_old, α, α_old, m, m_2, FLUID, dt)
    res_mom   .= f_mom(dl, A_cross, D, m, m_old, T_1, P_1, α_1, P, FLUID, dt)
    res_ene   .= f_ene(V, A_cross, T, T_old, α, α_old, m, m_old, T_1, α_1, m_1, T_2, α_2, m_2, Q_grv, FLUID, dt)
    res_state .= f_state(T, P, α, FLUID)
end

function fluid_vl!(res_mass, res_mom, res_ene, res_state, var, var_old, T_wall, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt

    dl = geo.vl.dl
    V = geo.vl.V
    A_cross = geo.vl.A_cross
    S = geo.vl.S_inner
    D_i = geo.vl.D_i
    N = geo.vl.N

    for i in 1:N
        T = var.vl.T[i]
        T_old = var_old.vl.T[i]
        T_w = T_wall.vl[i]
        P = var.vl.P[i]
        α = α_cropper(var.vl.α[i])
        α_old = α_cropper(var_old.vl.α[i])
        m = var.vl.m[i]
        m_old = var_old.vl.m[i]

        T_1 = (i == 1) ? var.grv.T[end] : var.vl.T[i-1]
        P_1 = (i == 1) ? var.grv.P[end] : var.vl.P[i-1]
        α_1 = (i == 1) ? var.grv.α[end] : var.vl.α[i-1]
        m_1 = (i == 1) ? var.grv.m[end] : var.vl.m[i-1]

        T_2 = (i == N) ? var.cond.T[1] : var.vl.T[i+1]
        # P_2 = (i == N) ? var.cond.P[1] : var.vl.P[i+1]
        α_2 = (i == N) ? var.cond.α[1] : var.vl.α[i+1]
        m_2 = (i == N) ? var.cond.m[1] : var.vl.m[i+1]

        # Q_wall
        h = calc_h_wall_f_1p(T, P, m, D_i, A_cross, FLUID)
        Q_wall = calc_Q_conv(h, S, T, T_w)

        # residuals
        res_mass[i]  = f_mass(V, T, T_old, α, α_old, m, m_2, FLUID, dt)
        res_mom[i]   = f_mom(dl, A_cross, D_i, m, m_old, T_1, P_1, α_1, P, FLUID, dt)
        res_ene[i]   = f_ene(V, A_cross, T, T_old, α, α_old, m, m_old, T_1, α_1, m_1, T_2, α_2, m_2, Q_wall, FLUID, dt)
        res_state[i] = f_state(T, P, α, FLUID) 
    end

end

function fluid_cond!(res_mass, res_mom, res_ene, res_state, var, var_old, T_wall, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt

    dl = geo.cond.dl
    V = geo.cond.V
    A_cross = geo.cond.A_cross
    S = geo.cond.S_inner
    D_i = geo.cond.D_i
    N = geo.cond.N

    for i in 1:N
        T = var.cond.T[i]
        T_old = var_old.cond.T[i]
        T_w = T_wall.cond[i]
        P = var.cond.P[i]
        α = α_cropper(var.cond.α[i])
        α_old = α_cropper(var_old.cond.α[i])
        m = var.cond.m[i]
        m_old = var_old.cond.m[i]

        T_1 = (i == 1) ? var.vl.T[end] : var.cond.T[i-1]
        P_1 = (i == 1) ? var.vl.P[end] : var.cond.P[i-1]
        α_1 = (i == 1) ? var.vl.α[end] : var.cond.α[i-1]
        m_1 = (i == 1) ? var.vl.m[end] : var.cond.m[i-1]

        T_2 = (i == N) ? var.ll.T[1] : var.cond.T[i+1]
        # P_2 = (i == N) ? var.ll.P[1] : var.cond.P[i+1]
        α_2 = (i == N) ? var.ll.α[1] : var.cond.α[i+1]
        m_2 = (i == N) ? var.ll.m[1] : var.cond.m[i+1]

        # Q_wall
        if α == 1.0 || α == 0.0
            h = calc_h_wall_f_1p(T, P, m, D_i, A_cross, FLUID)
        else
            h = calc_h_condense(T, P, α, m, D_i, A_cross, FLUID)
        end
        Q_wall = calc_Q_conv(h, S, T, T_w)

        # residuals
        res_mass[i]  = f_mass(V, T, T_old, α, α_old, m, m_2, FLUID, dt)
        res_mom[i]   = f_mom(dl, A_cross, D_i, m, m_old, T_1, P_1, α_1, P, FLUID, dt)
        res_ene[i]   = f_ene(V, A_cross, T, T_old, α, α_old, m, m_old, T_1, α_1, m_1, T_2, α_2, m_2, Q_wall, FLUID, dt)
        res_state[i] = f_state(T, P, α, FLUID) 
    end

end

function fluid_ll!(res_mass, res_mom, res_ene, res_state, var, var_old, T_wall, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt

    dl = geo.ll.dl
    V = geo.ll.V
    A_cross = geo.ll.A_cross
    S = geo.ll.S_inner
    D_i = geo.ll.D_i
    N = geo.ll.N

    for i in 1:N
        T = var.ll.T[i]
        T_old = var_old.ll.T[i]
        T_w = T_wall.ll[i]
        P = var.ll.P[i]
        α = α_cropper(var.ll.α[i])
        α_old = α_cropper(var_old.ll.α[i])
        m = var.ll.m[i]
        m_old = var_old.ll.m[i]

        T_1 = (i == 1) ? var.cond.T[end] : var.ll.T[i-1]
        P_1 = (i == 1) ? var.cond.P[end] : var.ll.P[i-1]
        α_1 = (i == 1) ? var.cond.α[end] : var.ll.α[i-1]
        m_1 = (i == 1) ? var.cond.m[end] : var.ll.m[i-1]

        T_2 = (i == N) ? var.cc.T[1] : var.ll.T[i+1]
        # P_2 = (i == N) ? var.cc.P : var.ll.P[i+1]
        α_2 = (i == N) ? var.cc.α[1] : var.ll.α[i+1]
        m_2 = (i == N) ? var.cc.m[1] : var.ll.m[i+1]

        # Q_wall
        if α == 1.0 || α == 0.0
            h = calc_h_wall_f_1p(T, P, m, D_i, A_cross, FLUID)
        else
            h = calc_h_condense(T, P, α, m, D_i, A_cross, FLUID)
        end
        Q_wall = calc_Q_conv(h, S, T, T_w)

        # residuals
        res_mass[i]  = f_mass(V, T, T_old, α, α_old, m, m_2, FLUID, dt)
        res_mom[i]   = f_mom(dl, A_cross, D_i, m, m_old, T_1, P_1, α_1, P, FLUID, dt)
        res_ene[i]   = f_ene(V, A_cross, T, T_old, α, α_old, m, m_old, T_1, α_1, m_1, T_2, α_2, m_2, Q_wall, FLUID, dt)
        res_state[i] = f_state(T, P, α, FLUID) 
    end

end

function fluid_cc(res_mass, res_mom, res_ene, res_state, m_evap, var, var_old, T_wall, geo, condition)
    FLUID = condition.fluid
    dt = condition.dt
    
    dl = geo.cc.dl
    V = geo.cc.V
    A_cross = geo.cc.A_cross
    D_i = geo.cc.D_i # 代表長さとして
    D_o = geo.cc.D_o

    T = var.cc.T[1]
    T_old = var_old.cc.T[1]
    T_w = T_wall.cc
    P = var.cc.P[1]
    α = α_cropper(var.cc.α[1])
    α_old = α_cropper(var_old.cc.α[1])
    m = var.cc.m[1]
    m_old = var_old.cc.m[1]

    T_1 = var.ll.T[end]
    P_1 = var.ll.P[end]
    α_1 = var.ll.α[end]
    m_1 = var.ll.m[end]

    T_2 = var.wick[1]
    α_2 = 0.0
    m_2 = m_evap

    # Q_wall
    k = calc_k_2p(T, α, FLUID)
    Q_wall = calc_Q_fluid(k, D_o, D_i, T, T_w)

    # residuals
    res_mass  .= f_mass(V, T, T_old, α, α_old, m_1, m_2, FLUID, dt)
    res_mom   .= f_mom(dl, A_cross, D_i, m, m_old, T_1, P_1, α_1, P, FLUID, dt)
    res_ene   .= f_ene(V, A_cross, T, T_old, α, α_old, m, m_old, T_1, α_1, m_1, T_2, α_2, m_2, Q_wall, FLUID, dt)
    res_state .= f_state(T, P, α, FLUID)
end