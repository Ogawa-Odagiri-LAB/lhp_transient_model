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

const M_REF = 1.0e-5   # [kg/s] 流量の基準
const Q_REF = 1.0      # [W]    熱量の基準
const T_REF = 1.0      # [K]    温度の基準
 
function residual_scales(condition)
    dt = condition.dt
    return (
        mass  = dt * M_REF,   # [kg]   f_mass の残差の単位
        mom   = M_REF,        # [kg/s] f_mom の残差の単位（下の修正版 f_mom を使う前提）
        ene   = dt * Q_REF,   # [J]    f_ene・ウィックの残差の単位
        state = T_REF,        # [K]    f_state の残差の単位
    )
end
 
function f_transient_solve!(f, scaled_var, p)
    geo = p.geo
    transient_condition = p.condition
    T_wall = p.T_wall
    var_old = p.var_old
    s = residual_scales(transient_condition)
 
    restored_var = restore_var(scaled_var, geo, transient_condition)
 
    # ウィック（エネルギー式のみ）
    m_evap, res_wick = fluid_wick(restored_var, var_old, T_wall, geo, transient_condition)
    f.wick .= res_wick ./ s.ene
 
    # グルーブ
    r_mass, r_mom, r_ene, r_state = fluid_grv(m_evap, restored_var, var_old, T_wall, geo, transient_condition)
    f.grv.α = r_mass  / s.mass
    f.grv.m = r_mom   / s.mom
    f.grv.T = r_ene   / s.ene
    f.grv.P = r_state / s.state
 
    # 蒸気管・凝縮器・液管（各関数が f のビューに直接書き込んだ後でスケーリング）
    fluid_vl!(f.vl.α, f.vl.m, f.vl.T, f.vl.P, restored_var, var_old, T_wall, geo, transient_condition)
    fluid_cond!(f.cond.α, f.cond.m, f.cond.T, f.cond.P, restored_var, var_old, T_wall, geo, transient_condition)
    fluid_ll!(f.ll.α, f.ll.m, f.ll.T, f.ll.P, restored_var, var_old, T_wall, geo, transient_condition)
    for c in (:vl, :cond, :ll)
        comp = getproperty(f, c)
        comp.α ./= s.mass
        comp.m ./= s.mom
        comp.T ./= s.ene
        comp.P ./= s.state
    end
 
    # CC
    r_mass, r_mom, r_ene, r_state = fluid_cc(m_evap, restored_var, var_old, T_wall, geo, transient_condition)
    f.cc.α = r_mass  / s.mass
    f.cc.m = r_mom   / s.mom
    f.cc.T = r_ene   / s.ene
    f.cc.P = r_state / s.state
 
    p.m_evap[1] = m_evap
 
    # save_f_to_csv(f, geo, "residuals")
    # save_var_to_csv(m_evap, restored_var, geo, "var") # debug
    return nothing
end

function report_residuals(f)
    println("最大残差（無次元、abstol と比較）")
    println("  wick : ene = ", round(maximum(abs, f.wick), sigdigits=3))
    for c in (:grv, :vl, :cond, :ll, :cc)
        comp = getproperty(f, c)
        println("  ", rpad(string(c), 5),
                ": mass = ",  round(maximum(abs, comp.α), sigdigits=3),
                ", mom = ",   round(maximum(abs, comp.m), sigdigits=3),
                ", ene = ",   round(maximum(abs, comp.T), sigdigits=3),
                ", state = ", round(maximum(abs, comp.P), sigdigits=3))
    end
end

# 前回の計算結果（CSV）を削除する。figs/ などのサブフォルダと画像は残す
function clear_old_results(res_dir = "./res/vars_transient/")
    isdir(res_dir) || return nothing
    n_removed = 0
    for f in readdir(res_dir)
        path = joinpath(res_dir, f)
        if isfile(path) && endswith(f, ".csv")
            rm(path)
            n_removed += 1
        end
    end
    n_removed > 0 && @info "古い計算結果の CSV を $(n_removed) 個削除しました（$(res_dir)）"
    return nothing
end

function solve_transient_lhp!(geo, transient_condition)
    # condition
    t_init = transient_condition.t_init
    t_end = transient_condition.t_end
    dt = transient_condition.dt

    # 古い結果の削除（別条件の結果が混ざらないようにする）
    clear_old_results("./res/vars_transient/")

    # initialize
    # var = set_values_on_var_from_steady(geo, 20)

    # var = initialize_transient_var(geo)
    # set_values_on_var!(var, transient_condition)
    var = set_smooth_values_on_var(geo, transient_condition)
    var_old = copy(var)
    lb = copy(var)
    ub = copy(var)
    set_values_on_lb!(lb)
    set_values_on_ub!(ub)

    # T_wall = set_values_on_T_wall_from_steady(geo, transient_condition)
    T_wall = initialize_T_wall(geo)
    set_values_on_T_wall!(T_wall, geo, transient_condition)
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
        sol = solve(prob, algorithm; show_trace = Val(true), abstol=1.0e-3, maxiters=40)
        if !SciMLBase.successful_retcode(sol)
            @warn "t = $t s で未収束" sol.retcode
        end
        report_residuals(sol.resid)

        T_wall_old .= params.T_wall
        scaled_var .= sol.u
        params.var_old .= restore_var(scaled_var, params.geo, params.condition)
        m_evap, _ = fluid_wick(params.var_old, params.var_old, params.T_wall, params.geo, params.condition)
        params.m_evap[1] = m_evap

        save_T_wall_to_csv(params.T_wall, geo, "T_wall_$(t)s")
        save_var_to_csv(params.m_evap[1], params.var_old, geo, "var_$(t)s")
        save_var_to_csv(params.m_evap[1], scaled_var, geo, "scaled_var_$(t)s")
    end
end