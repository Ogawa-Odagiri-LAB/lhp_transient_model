function wall_ec(m_evap_old, var_old, T_wall_old, geo, condition)
    Q_apply = condition.Q_apply

    # Q_ev
    Q_ev = calc_Q_ev(var_old, T_wall_old, geo, condition)

    # Q_grv
    Q_grv = calc_Q_grv(m_evap_old, var_old, T_wall_old, geo, condition)

    # Q_ec_cc
    Q_ec_cc = calc_Q_ec_cc(T_wall_old, geo)

    # Q_ec_vl
    Q_ec_vl = calc_Q_ec_vl(var_old, T_wall_old, geo, condition)

    # specific heat
    C = ((geo.grv.D_o^2 - geo.grv.D_i^2) * 0.25 * geo.grv.L) * ρ_sus304(T_wall_old.evap) * c_sus304(T_wall_old.evap) 

    dTdt = (Q_apply - Q_ev - Q_grv - Q_ec_cc - Q_ec_vl) / C
    return dTdt
end

function wall_common(var_old, T_wall_old, geo, condition, comp1, comp, comp2, T_a, T_z)
    # 引数
    # wall_common(var_old, T_wall_old, geo, condition, comp1（1個前のコンポーネントのSymbol）, comp（計算したいコンポーネントのSymbol）, comp2（1個次のコンポーネントのSymbol）， T_a（1個前の壁面温度）, T_z（1個次の壁面温度）)
    FLUID = condition.fluid
    if comp == :cond
        T_amb = condition.T_cond
    else
        T_amb = condition.T_amb
    end
    T_1 = 0.0
    T_2 = 0.0
    Q_amb = 0.0

    geo_part_1 = getproperty(geo, comp1)
    geo_part = getproperty(geo, comp)
    geo_part_2 = getproperty(geo, comp2)
    T_wall_old_part = getproperty(T_wall_old, comp)
    var_old_part = getproperty(var_old, comp)

    dTdt = zeros(geo_part.N)

    V = (geo_part.D_o^2 - geo_part.D_i^2) * 0.25 * geo_part.dl

    for i in 1:geo_part.N
        T = T_wall_old_part[i]
        T_f = var_old_part.T[i]
        P = var_old_part.P[i]
        α = var_old_part.α[i]
        m = var_old_part.m[i]

        # Q_cond
        k = k_sus304(T)
        G_1 = calc_G_pipe(k, geo_part.D_o, geo_part.D_i, geo_part.dl)
        G   = calc_G_pipe(k, geo_part.D_o, geo_part.D_i, geo_part.dl)
        G_2 = calc_G_pipe(k, geo_part.D_o, geo_part.D_i, geo_part.dl)
        
        if i == 1
            k_a = k_sus304(T_a)
            G_1 = calc_G_pipe(k_a, geo_part_1.D_o, geo_part_1.D_i, geo_part_1.dl)
            T_1 = T_a
        elseif i == geo_part.N
            k_z = k_sus304(T_z)
            G_2 = calc_G_pipe(k_z, geo_part_2.D_o, geo_part_2.D_i, geo_part_2.dl)
            T_2 = T_z
        end

        T_1 = (i == 1) ? T_a : T_wall_old_part[i-1]
        T_2 = (i == geo_part.N) ? T_z : T_wall_old_part[i+1]

        if comp == :cond
            Q_amb = calc_Q_wall_amb(k, geo_part.D_o, geo_part.D_i, T_amb, T)
        else
            Q_amb = 0.0
        end
        Q_cond = calc_Q_wall_cond(G_1, G, T_1, T) + calc_Q_wall_cond(G_2, G, T_2, T) +  Q_amb # 管壁内の熱伝導と凝縮器または周囲の温度の影響（熱伝導）で考慮

        # Q_conv
        if α == 0.0 || α == 1.0
            h = calc_h_wall_f_1p(T_f, P, m, geo_part.D_i, geo_part.A_cross, FLUID)
        else
            h = calc_h_condense(T_f, P, α, m, geo_part.D_i, geo_part.A_cross, FLUID)
        end

        Q_conv = calc_Q_conv(h, geo_part.S_inner, T_f, T)
        
        # specific heat
        C = V * ρ_sus304(T_wall_old_part[i]) * c_sus304(T_wall_old_part[i])

        dTdt[i] = (Q_cond + Q_conv) / C
    end
    return dTdt
end

function wall_cc(var_old, T_wall_old, geo, condition)
    T_cc = T_wall_old.cc
    k = k_sus304(T_cc)

    # Q_ec_cc
    Q_ec_cc = calc_Q_ec_cc(T_wall_old, geo)

    # Q_cc_ll
    Q_cc_ll = calc_Q_cc_ll(T_wall_old, geo)

    # Q_amb
    Q_amb = 0.0

    # Q_fluid
    Q_fluid = calc_Q_fluid(k, geo.cc.D_o, geo.cc.D_i, var_old.cc.T[1], T_cc)

    # specific heat
    C = ((geo.cc.D_o^2 - geo.cc.D_i^2) * 0.25 * geo.cc.L) * ρ_sus304(T_wall_old.cc) * c_sus304(T_wall_old.cc) 

    dTdt = (Q_ec_cc + Q_cc_ll + Q_amb + Q_fluid) / C
    return dTdt
end

function Euler_wall!(m_evap_old, var_old, T_wall, T_wall_old, geo, condition) # 分割Nに対してdtが大きすぎると計算が破綻する．
    dt = condition.dt

    dTdt_evap = wall_ec(m_evap_old, var_old, T_wall_old, geo, condition)
    dTdt_vl   = wall_common(var_old, T_wall_old, geo, condition, :grv,   :vl,    :cond, T_wall_old.evap,       T_wall_old.cond[1])
    dTdt_cond = wall_common(var_old, T_wall_old, geo, condition, :vl,    :cond,  :ll,   T_wall_old.vl[end],    T_wall_old.ll[1])
    dTdt_ll   = wall_common(var_old, T_wall_old, geo, condition, :cond,  :ll,    :cc,   T_wall_old.cond[end],  T_wall_old.cc)
    dTdt_cc   = wall_cc(var_old, T_wall_old, geo, condition)

    T_wall.evap  = T_wall_old.evap  + dt  * dTdt_evap
    T_wall.vl   .= T_wall_old.vl   .+ dt .* dTdt_vl
    T_wall.cond .= T_wall_old.cond .+ dt .* dTdt_cond
    T_wall.ll   .= T_wall_old.ll   .+ dt .* dTdt_ll
    T_wall.cc    = T_wall_old.cc    + dt  * dTdt_cc
end

function calc_dTdt_wall(m_evap_old, var_old, T_wall, geo, condition)
    # T_wall_now（現在の暫定温度）をベースに各部の微分値を計算
    dTdt_evap = wall_ec(m_evap_old, var_old, T_wall, geo, condition)
    dTdt_vl   = wall_common(var_old, T_wall, geo, condition, :grv,   :vl,    :cond, T_wall.evap,       T_wall.cond[1])
    dTdt_cond = wall_common(var_old, T_wall, geo, condition, :vl,    :cond,  :ll,   T_wall.vl[end],    T_wall.ll[1])
    dTdt_ll   = wall_common(var_old, T_wall, geo, condition, :cond,  :ll,    :cc,   T_wall.cond[end],  T_wall.cc)
    dTdt_cc   = wall_cc(var_old, T_wall, geo, condition)

    # T_wall と同じ構造の ComponentVector として微分値をまとめて返す
    # (内部配列を新しく確保して返します)
    return ComponentVector(
        evap = dTdt_evap,
        vl   = dTdt_vl,
        cond = dTdt_cond,
        ll   = dTdt_ll,
        cc   = dTdt_cc
    )
end

# function RK4_wall!(m_evap_old, var_old, T_wall, T_wall_old, geo, condition)
#     dt = condition.dt

#     # k1
#     T_tmp1 = T_wall_old
#     k1 = calc_dTdt_wall(m_evap_old, var_old, T_tmp1, geo, condition)

#     # k2
#     T_tmp2 = T_wall_old .+ (dt / 2.0) .* k1
#     k2 = calc_dTdt_wall(m_evap_old, var_old, T_tmp2, geo, condition)

#     # k3
#     T_tmp3 = T_wall_old .+ (dt / 2.0) .* k2
#     k3 = calc_dTdt_wall(m_evap_old, var_old, T_tmp3, geo, condition)

#     # k4
#     T_tmp4 = T_wall_old .+ dt .* k3
#     k4 = calc_dTdt_wall(m_evap_old, var_old, T_tmp4, geo, condition)

#     T_wall .= T_wall_old .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
# end

function RK4_wall!(m_evap_old, var_old, T_wall, T_wall_old, geo, condition)
    dt_main = condition.dt
    
    # 凝縮器の巨大な熱伝達率に耐えるため、内部だけで例えば 100 分割する
    N_sub = 100
    dt_sub = dt_main / N_sub
    
    T_current = copy(T_wall_old)
    
    for i in 1:N_sub
        k1 = calc_dTdt_wall(m_evap_old, var_old, T_current, geo, condition)
        
        T_tmp2 = T_current .+ (dt_sub / 2.0) .* k1
        k2 = calc_dTdt_wall(m_evap_old, var_old, T_tmp2, geo, condition)
        
        T_tmp3 = T_current .+ (dt_sub / 2.0) .* k2
        k3 = calc_dTdt_wall(m_evap_old, var_old, T_tmp3, geo, condition)
        
        T_tmp4 = T_current .+ dt_sub .* k3
        k4 = calc_dTdt_wall(m_evap_old, var_old, T_tmp4, geo, condition)
        
        T_current .= T_current .+ (dt_sub / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
    end
    
    T_wall .= T_current
end