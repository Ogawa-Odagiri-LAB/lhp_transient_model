function restore_var(scaled_var, geo, condition)
    T_scale = condition.T_scale
    P_scale = condition.P_scale
    α_scale = condition.α_scale
    m_scale = condition.m_scale

    restored_var = initialize_transient_var(geo)

    # T
    restored_var.wick   .= scaled_var.wick   .* T_scale
    restored_var.grv.T   = scaled_var.grv.T   * T_scale
    restored_var.vl.T   .= scaled_var.vl.T   .* T_scale
    restored_var.cond.T .= scaled_var.cond.T .* T_scale
    restored_var.ll.T   .= scaled_var.ll.T   .* T_scale
    restored_var.cc.T    = scaled_var.cc.T    * T_scale

    # P
    # wick
    restored_var.grv.P   = scaled_var.grv.P   * P_scale
    restored_var.vl.P   .= scaled_var.vl.P   .* P_scale
    restored_var.cond.P .= scaled_var.cond.P .* P_scale
    restored_var.ll.P   .= scaled_var.ll.P   .* P_scale
    restored_var.cc.P    = scaled_var.cc.P    * P_scale

    # α
    # wick
    restored_var.grv.α   = scaled_var.grv.α   * α_scale
    restored_var.vl.α   .= scaled_var.vl.α   .* α_scale
    restored_var.cond.α .= scaled_var.cond.α .* α_scale
    restored_var.ll.α   .= scaled_var.ll.α   .* α_scale
    restored_var.cc.α    = scaled_var.cc.α    * α_scale

    # m
    # wick
    restored_var.grv.m   = scaled_var.grv.m   * m_scale
    restored_var.vl.m   .= scaled_var.vl.m   .* m_scale
    restored_var.cond.m .= scaled_var.cond.m .* m_scale
    restored_var.ll.m   .= scaled_var.ll.m   .* m_scale
    restored_var.cc.m    = scaled_var.cc.m    * m_scale

    return restored_var
end

function scale_var(var, geo, condition)
    T_scale = condition.T_scale
    P_scale = condition.P_scale
    α_scale = condition.α_scale
    m_scale = condition.m_scale

    scaled_var = initialize_transient_var(geo)

    # T
    scaled_var.wick   .= var.wick   ./ T_scale
    scaled_var.grv.T   = var.grv.T   / T_scale
    scaled_var.vl.T   .= var.vl.T   ./ T_scale
    scaled_var.cond.T .= var.cond.T ./ T_scale
    scaled_var.ll.T   .= var.ll.T   ./ T_scale
    scaled_var.cc.T    = var.cc.T    / T_scale

    # P
    # wick
    scaled_var.grv.P   = var.grv.P   / P_scale
    scaled_var.vl.P   .= var.vl.P   ./ P_scale
    scaled_var.cond.P .= var.cond.P ./ P_scale
    scaled_var.ll.P   .= var.ll.P   ./ P_scale
    scaled_var.cc.P    = var.cc.P    / P_scale

    # α
    # wick
    scaled_var.grv.α   = var.grv.α   / α_scale
    scaled_var.vl.α   .= var.vl.α   ./ α_scale
    scaled_var.cond.α .= var.cond.α ./ α_scale
    scaled_var.ll.α   .= var.ll.α   ./ α_scale
    scaled_var.cc.α    = var.cc.α    / α_scale

    # m
    # wick
    scaled_var.grv.m   = var.grv.m   / m_scale
    scaled_var.vl.m   .= var.vl.m   ./ m_scale
    scaled_var.cond.m .= var.cond.m ./ m_scale
    scaled_var.ll.m   .= var.ll.m   ./ m_scale
    scaled_var.cc.m    = var.cc.m    / m_scale

    return scaled_var
end

function f_transient_solve!(f, scaled_var, p)
    geo = p.geo
    transient_condition = p.condition
    T_wall = p.T_wall
    var_old = p.var_old

    restored_var = restore_var(scaled_var, geo, transient_condition)

    m_evap, f.wick = fluid_wick(restored_var, var_old, T_wall, geo, transient_condition)
    f.grv.α, f.grv.m, f.grv.T, f.grv.P = fluid_grv(m_evap, restored_var, var_old, T_wall, geo, transient_condition)
    fluid_vl!(f.vl.α, f.vl.m, f.vl.T, f.vl.P, restored_var, var_old, T_wall, geo, transient_condition)
    fluid_cond!(f.cond.α, f.cond.m, f.cond.T, f.cond.P, restored_var, var_old, T_wall, geo, transient_condition)
    fluid_ll!(f.ll.α, f.ll.m, f.ll.T, f.ll.P, restored_var, var_old, T_wall, geo, transient_condition)
    f.cc.α, f.cc.m, f.cc.T, f.cc.P = fluid_cc(m_evap, restored_var, var_old, T_wall, geo, transient_condition)

    p.m_evap[1] = m_evap

    # save_f_to_csv(f, geo, "residuals")
    # save_var_to_csv(m_evap, restored_var, geo, "var") # debug
end

function solve_transient_lhp!(geo, transient_condition)
    # condition
    t_init = transient_condition.t_init
    t_end = transient_condition.t_end
    dt = transient_condition.dt

    # initialize
    var = set_values_on_var_from_steady(geo, 20)
    var_old = copy(var)
    lb = copy(var)
    ub = copy(var)
    set_values_on_lb!(lb)
    set_values_on_ub!(ub)

    T_wall = set_values_on_T_wall_from_steady(geo, transient_condition)
    T_wall_old = copy(T_wall)

    # scaling
    scaled_var = scale_var(var, geo, transient_condition)
    lb_scaled = scale_var(lb, geo, transient_condition)
    ub_scaled = scale_var(ub, geo, transient_condition)

    # Parameters
    params = (
        geo = geo,
        condition = transient_condition, # struct もそのまま格納できます
        m_evap = [0.0],                    # 外側で初期値（0.0など）を入れておく
        T_wall = T_wall,
        var_old = var_old
    )

    prob = NonlinearProblem(f_transient_solve!, scaled_var, params; lb=lb_scaled, ub=ub_scaled) # debug

    # main
    for t in t_init:dt:t_end
        @show t
        RK4_wall!(params.m_evap[1], params.var_old, params.T_wall, T_wall_old, params.geo, params.condition)

        # prob = NonlinearProblem(f_transient_solve!, scaled_var, params; lb=lb_scaled, ub=ub_scaled, jac_prototype=jac_prototype)
        prob = remake(prob; u0 = scaled_var, p = params) # debug
        algorithm = TrustRegion(
            autodiff     = AutoFiniteDiff(),
            # concrete_jac = Val(true),
            vjp_autodiff = nothing,        
            jvp_autodiff = nothing          
        )
        # algorithm = LevenbergMarquardt(
        #     autodiff = AutoFiniteDiff(),
        #     vjp_autodiff = nothing
        # )
        sol = solve(prob, algorithm; show_trace = Val(true), abstol=2.5e-3, maxiters=20)

        T_wall_old .= params.T_wall
        scaled_var .= sol.u
        params.var_old .= restore_var(scaled_var, params.geo, params.condition)

        save_T_wall_to_csv(params.T_wall, geo, "T_wall_$(t)s")
        save_var_to_csv(params.m_evap[1], params.var_old, geo, "var_$(t)s")
        save_var_to_csv(params.m_evap[1], scaled_var, geo, "scaled_var_$(t)s")
    end
end