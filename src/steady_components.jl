function calc_grv!(u, vars, geo, h_ev, FLUID)
    n_grv = geo.grv.n_grv
    dl = geo.grv.dl
    width = geo.grv.width
    height = geo.grv.height
    A_grv_cross = geo.grv.A_grv_cross
    A_ec_wick = geo.grv.A_ec_wick
    v_grv = width * height * dl

    m = calc_mass_flow_init(u.T_ev, u.T_ec, A_ec_wick, h_ev, FLUID)
    m_grv = m / n_grv

    vars.grv.T[1] = u.T_ev
    vars.grv.P[1] = calc_Psat(u.T_ev, FLUID)

    ρ(idx) = calc_ρ_1p(vars.grv.T[idx], vars.grv.P[idx], FLUID)
    cp(idx) = calc_cp_1p(vars.grv.T[idx], vars.grv.P[idx], FLUID)
    h_ec_gr(idx) = calc_h_ec_gr(vars.grv.T[idx], vars.grv.P[idx], m_grv, height, A_grv_cross, geo.grv.L, FLUID)
    τ(idx) = calc_τ_wall(vars.grv.T[idx], vars.grv.P[idx], m_grv, height, A_grv_cross, FLUID)
    Q_gr_t = zeros(geo.grv.N)

    # println("grv!")
    for i in 2:geo.grv.N
        Q_gr_t[i-1]   = h_ec_gr(i-1) * width * (u.T_ec - vars.grv.T[i-1])
        vars.grv.T[i] = vars.grv.T[i-1] + (Q_gr_t[i-1] + (width + height) * 2 * τ(i-1) * m_grv / (ρ(i-1) * A_grv_cross)) * (dl / (m_grv * cp(i-1)))
        vars.grv.P[i] = vars.grv.P[i-1] - ((width + height) * 2 * τ(i-1) / A_grv_cross) * dl
        vars.grv.mass[i] = ρ(i) * v_grv
        # println("T=$(vars.grv.T[i])")
    end
    Q_gr = sum(Q_gr_t) * dl
    return m, Q_gr
end

function calc_vl!(vars, m, geo, FLUID)
    dl = geo.vl.dl
    D_i = geo.vl.D_i
    A_cross = geo.vl.A_cross
    v_vl = A_cross * dl

    vars.vl.T[1] = vars.grv.T[geo.grv.N]
    vars.vl.P[1] = vars.grv.P[geo.grv.N]

    ρ(idx) = calc_ρ_1p(vars.vl.T[idx], vars.vl.P[idx], FLUID)
    cp(idx) = calc_cp_1p(vars.vl.T[idx], vars.vl.P[idx], FLUID)
    τ(idx) = calc_τ_wall(vars.vl.T[idx], vars.vl.P[idx], m, D_i, A_cross, FLUID)

    # println("vl!")
    for i in 2:geo.vl.N
        vars.vl.T[i] = vars.vl.T[i-1] + (π * D_i * τ(i-1) * m / (ρ(i-1) * A_cross)) *(dl / (m * cp(i-1)))
        vars.vl.P[i] = vars.vl.P[i-1] - (π * D_i * τ(i-1) / A_cross) * dl
        vars.vl.mass[i] = ρ(i) * v_vl
        # println("T=$(vars.vl.T[i])")
    end
end

function calc_cond!(vars, m, geo, T_cond, FLUID)
    dl = geo.cond.dl
    D_i = geo.cond.D_i
    A_cross = geo.cond.A_cross
    v_cond = A_cross * dl

    vars.cond.T[1] = vars.vl.T[geo.vl.N]
    vars.cond.P[1] = vars.vl.P[geo.vl.N]

    ρ_1p(idx) = calc_ρ_1p(vars.cond.T[idx], vars.cond.P[idx], FLUID)
    cp_1p(idx) = calc_cp_1p(vars.cond.T[idx], vars.cond.P[idx], FLUID)
    τ_1p(idx) = calc_τ_wall(vars.cond.T[idx], vars.cond.P[idx], m, D_i, A_cross, FLUID)
    h_wall_f_1p(idx) = calc_h_wall_f_1p(vars.cond.T[idx], vars.cond.P[idx], m, D_i, A_cross, FLUID)

    ρ_2p(idx) = calc_ρ_2p(vars.cond.T[idx], vars.cond.P[idx], FLUID)
    cp_2p(idx) = calc_cp_2p(vars.cond.T[idx], vars.cond.P[idx], FLUID)
    α(idx) = calc_void_rate(vars.cond.T[idx], vars.cond.x[idx], FLUID)
    λ(idx) = calc_latent_heat(vars.cond.T[idx], FLUID)
    ΔP_2p(idx) = calc_ΔP_2p(vars.cond.T[idx], α(idx), m, D_i, A_cross, FLUID)
    h_wall_f_2p(idx) = calc_h_wall_f_2p(vars.cond.T[idx], α(idx), m, D_i, A_cross, FLUID)
    tp_flag::Symbol = :v # 二相状態の計算に入るかどうか, 初期値は:v(vapor), 二相状態になったら:tp(TwoPhase), 液相になったら:l(liquid)

    # println("cond!")

    for i in 2:geo.cond.N
        # 通常計算
        if TPcheck(vars.cond.T[i-1], vars.cond.P[i-1], FLUID) == :tp && tp_flag == :v # 二相流の入り口を検知(最初は気相)
            tp_flag = :tp
            # println("vapor ends!")
        elseif tp_flag != :tp
            vars.cond.T[i] = vars.cond.T[i-1] + (h_wall_f_1p(i-1) * π * D_i * (T_cond - vars.cond.T[i-1]) + π * D_i * τ_1p(i-1) * m / (ρ_1p(i-1) * A_cross)) * (dl / (m * cp_1p(i-1)))
            # println("conv=$((h_wall_f_1p(i-1) * π * D_i * (T_cond - vars.cond.T[i-1]))*(dl / (m * cp_1p(i-1)))), m=$(m)")
            vars.cond.P[i] = vars.cond.P[i-1] - (π * D_i * τ_1p(i-1) / A_cross) * dl
            vars.cond.mass[i] = ρ_1p(i) * v_cond
        end

        # 二相流に関する計算
        if tp_flag == :tp # 二相流領域
            # println("tp!")
            vars.cond.x[i] = vars.cond.x[i-1] + (h_wall_f_2p(i-1) * π * D_i * (T_cond - vars.cond.T[i-1]) / (m * λ(i-1))) * dl
            vars.cond.P[i] = vars.cond.P[i-1] - ΔP_2p(i-1) * dl
            vars.cond.T[i] = calc_Tsat(vars.cond.P[i], FLUID)
            vars.cond.mass[i] = ρ_2p(i) * v_cond
            vars.cond.α[i] = calc_void_rate(vars.cond.T[i], vars.cond.x[i], FLUID)
            if vars.cond.x[i] < 0.0 # 直前の計算結果から二相流の出口を検知
                tp_flag = :l # 液相開始
                # println("liquid starts!")
                vars.cond.x[i] = 0.0
            end
        elseif tp_flag == :l # 液相領域
            vars.cond.x[i] = 0.0
        end
        # println("T=$(vars.cond.T[i]), P=$(vars.cond.P[i]), x=$(vars.cond.x[i])")
    end
end

function calc_ll!(vars, m, geo, FLUID)
    dl = geo.ll.dl
    D_i = geo.ll.D_i
    A_cross = geo.ll.A_cross
    v_ll = A_cross * dl

    vars.ll.T[1] = vars.cond.T[geo.cond.N]
    vars.ll.P[1] = vars.cond.P[geo.cond.N]

    ρ(idx) = calc_ρ_1p(vars.ll.T[idx], vars.ll.P[idx], FLUID)
    cp(idx) = calc_cp_1p(vars.ll.T[idx], vars.ll.P[idx], FLUID)
    τ(idx) = calc_τ_wall(vars.ll.T[idx], vars.ll.P[idx], m, D_i, A_cross, FLUID)

    # println("ll!")
    for i in 2:geo.ll.N
        vars.ll.T[i] = vars.ll.T[i-1] + (π * D_i * τ(i-1) * m / (ρ(i-1) * A_cross)) *(dl / (m * cp(i-1)))
        vars.ll.P[i] = vars.ll.P[i-1] - (π * D_i * τ(i-1) / A_cross) * dl
        vars.ll.mass[i] = ρ(i) * v_ll
        # println("T=$(vars.ll.T[i]), P=$(vars.ll.P[i])")
    end
end

function calc_bay!(u, vars, m, geo, FLUID)
    dl = geo.bay.dl
    D_i = geo.bay.D_i
    A_cross = geo.bay.A_cross

    vars.bay.T[1] = vars.ll.T[geo.ll.N]
    vars.bay.P[1] = vars.ll.P[geo.ll.N]

    ρ_1p(idx) = calc_ρ_1p(vars.bay.T[idx], vars.bay.P[idx], FLUID)
    cp_1p(idx) = calc_cp_1p(vars.bay.T[idx], vars.bay.P[idx], FLUID)
    τ_1p(idx) = calc_τ_wall(vars.bay.T[idx], vars.bay.P[idx], m, D_i, A_cross, FLUID)
    h_wall_f_1p(idx) = calc_h_wall_f_1p(vars.bay.T[idx], vars.bay.P[idx], m, D_i, A_cross, FLUID)

    ρ_2p(idx) = calc_ρ_2p(vars.bay.T[idx], vars.bay.P[idx], FLUID)
    cp_2p(idx) = calc_cp_2p(vars.bay.T[idx], vars.bay.P[idx], FLUID)
    α(idx) = calc_void_rate(vars.bay.T[idx], vars.bay.x[idx], FLUID)
    λ(idx) = calc_latent_heat(vars.bay.T[idx], FLUID)
    ΔP_2p(idx) = calc_ΔP_2p(vars.bay.T[idx], α(idx), m, D_i, A_cross, FLUID)
    h_wall_f_2p(idx) = calc_h_wall_f_2p(vars.bay.T[idx], α(idx), m, D_i, A_cross, FLUID)
    tp_flag::Symbol = :l # 二相状態の計算に入るかどうか, 初期値は:l(liquid), 二相状態になったら:tp(TwoPhase), 気相になったら:v(vapor)

    # println("bay!")

    for i in 2:geo.bay.N
        # 通常計算
        if TPcheck(vars.bay.T[i-1], vars.bay.P[i-1], FLUID) == :tp && tp_flag == :l # 二相流の入り口を検知(最初は液相)
            tp_flag = :tp
            # println("liquid ends!")
        elseif tp_flag != :tp
            vars.bay.T[i] = vars.bay.T[i-1] + (h_wall_f_1p(i-1) * π * D_i * (u.T_cc - vars.bay.T[i-1]) + π * D_i * τ_1p(i-1) * m / (ρ_1p(i-1) * A_cross)) * (dl / (m * cp_1p(i-1)))
            vars.bay.P[i] = vars.bay.P[i-1] - (π * D_i * τ_1p(i-1) / A_cross) * dl
        end

        # 二相流に関する計算
        if tp_flag == :tp # 二相流領域
            # println("tp!")
            vars.bay.x[i] = vars.bay.x[i-1] + (h_wall_f_2p(i-1) * π * D_i * (u.T_cc - vars.bay.T[i-1]) / (m * λ(i-1))) * dl
            vars.bay.P[i] = vars.bay.P[i-1] - ΔP_2p(i-1) * dl
            vars.bay.T[i] = calc_Tsat(vars.bay.P[i], FLUID)
            if vars.bay.x[i] > 1.0 # 直前の計算結果から二相流の出口を検知
                tp_flag = :v # 気相開始
                # println("vapor starts!")
                vars.bay.x[i] = 1.0
            end
        elseif tp_flag == :v # 気相領域
            vars.bay.x[i] = 1.0
        end
        # println("T=$(vars.bay.T[i]), P=$(vars.bay.P[i]), x=$(vars.bay.x[i])")
    end
end

function calc_residual!(f_res, u, vars, geo, condition, Q_nom)
    FLUID = condition.fluid
    h_ev = condition.h_ev
    T_cond = condition.T_cond
    Q_apply = condition.Q_apply

    # 独立変数
    T_ev = u.T_ev
    T_ec = u.T_ec

    f_res_t = zeros(3)

    m, Q_gr = calc_grv!(u, vars, geo, h_ev, FLUID) # グルーブでの温度計算と全体の質量流量出力
    calc_vl!(vars, m, geo, FLUID) # 蒸気管計算
    calc_cond!(vars, m, geo, T_cond, FLUID) # 凝縮器計算
    calc_ll!(vars, m, geo, FLUID) # 液管計算
    # calc_bay!(u, vars, m, geo) # バイヨネット管計算

    # 蒸発器の熱平衡
    A_ec_wick = geo.grv.A_ec_wick
    Q_ev = h_ev * A_ec_wick * (T_ec - T_ev)
    Q_leak = Q_apply - (Q_ev + Q_gr)

    # CCへの熱漏れ
    k_eff = calc_k_eff(T_ev, FLUID)
    L_wick = geo.grv.L
    A_ec_cc = (π / 4) * (geo.cc.D_o^2 - geo.cc.D_i^2)
    l_ec_cc = (geo.grv.L + geo.cc.L) / 2
    G_ec_cc = k_sus304(T_ec) * (A_ec_cc / l_ec_cc)

    # CCの温度バランス
    dTdP = calc_dT_dP_sat(T_ev, FLUID)
    ΔP_grv_cc = vars.grv.P[2] - vars.ll.P[geo.ll.N]
    T_cc = T_ev - dTdP * ΔP_grv_cc

    # 熱漏れとサブクールのバランス
    cp = PropsSI("C", "T", T_cc, "Q", 0.0, FLUID)

    # 残差方程式
    f_res_t[1] = Q_leak - (2 * π * k_eff * L_wick) * (T_ev - T_cc) / (log(geo.grv.D_o / geo.grv.D_i)) - (G_ec_cc * (T_ec - T_cc))# CCへの熱漏れ
    f_res_t[2] = Q_leak - m * cp * (T_cc - vars.ll.T[geo.ll.N]) # 熱漏れとサブクールのバランス
    println("u: T_ev=$(T_ev), T_ec=$(T_ec)")
    println("Q_ev=$(Q_ev), Q_leak=$(Q_leak), Q_gr=$(Q_gr), T_cc=$(T_cc), m=$(m)")

    # save_vars_to_csv(vars, geo) # 外部で保存できる
    save_residual_to_csv(f_res_t, condition)
    save_Q_T_m_to_csv(Q_ev, Q_gr, Q_leak, T_ev, T_ec, T_cc, m)
    
    f_res[1] = f_res_t[1] / Q_nom
    f_res[2] = f_res_t[2] / Q_nom
end

function f_nonlinear_solve!(f_res, u_scaled, p)
    u_physical = ComponentArray(
        T_ev    = u_scaled.T_ev * p.Q_nom,
        T_ec    = u_scaled.T_ec * p.Q_nom
    )

    calc_residual!(f_res, u_physical, p.vars, p.geo, p.condition, p.Q_nom)
end

function estimate_cc_void_ratio!(α, vars, geo, condition, mass_charged=0.03)
    # 各コンポーネントの質量を計算
    # グルーブ質量
    mass_grv    = sum(vars.grv.mass) * geo.grv.n_grv
    mass_vl     = sum(vars.vl.mass)
    mass_cond   = sum(vars.cond.mass)
    mass_ll     = sum(vars.ll.mass)

    mass_cc = mass_charged - (mass_grv + mass_vl + mass_cond + mass_ll) # [kg]
    T_cc = only(load_T_cc_only())
    println("α=$(α)")
    vars.cc.x[geo.cc.N] = calc_vapor_quality(T_cc, α, condition.fluid)
    vars.cc.α[geo.cc.N] = α
    vars.cc.T[geo.cc.N] = T_cc
    vars.cc.P[geo.cc.N] = calc_Psat(T_cc, condition.fluid)
    return geo.cc.A_cross * geo.cc.L - mass_cc / calc_ρ_2p(T_cc, α, condition.fluid)
end

function f_cc_void_solve(α, p)
    return estimate_cc_void_ratio!(α, p.vars, p.geo, p.condition, p.mass_charged)
end

function solve_steady_lhp!(vars, geo, condition)
    Q_nom = 373.15 # 正規化用

    # 初期条件
    initial_struct = opt_u_with_Q(condition.Q_apply) # Q_applyを用いて初期値を定義
    initial_u_scaled = ComponentArray(
        T_ev = initial_struct.T_ev / Q_nom,
        T_ec = initial_struct.T_ec / Q_nom,
    )

    # 制約条件
    # 下限値 (Lower Bounds) の設定
    lb = ComponentArray(
        T_ev   = 253.15,    
        T_ec   = 253.15,
    )

    # 上限値 (Upper Bounds) の設定
    # 特に上限がない変数は「Inf（無限大）」にしておけばOK
    ub = ComponentArray(
        T_ev   = 360.0,       
        T_ec   = 373.15
    )

    lb_scaled = lb ./ ComponentArray(T_ev=Q_nom, T_ec=Q_nom)
    ub_scaled = ub ./ ComponentArray(T_ev=Q_nom, T_ec=Q_nom)

    params = (vars=vars, geo=geo, condition=condition, Q_nom=Q_nom)

    prob = NonlinearProblem(f_nonlinear_solve!, initial_u_scaled, params; lb=lb_scaled, ub=ub_scaled)
    sol = solve(prob, TrustRegion(autodiff=AutoFiniteDiff()); show_trace = Val(true), abstol=1.0e-5)

    alpha_interval = (0.0, 1.0)
    params_CC = (vars=vars, geo=geo, condition=condition, mass_charged=0.03)
    prob_CC = IntervalNonlinearProblem(f_cc_void_solve, alpha_interval, params_CC)
    sol_CC = solve(prob_CC, Bisection(); show_trace = Val(true), abstol=1.0e-6)
end

