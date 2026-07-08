function save_T_wall_to_csv(T_wall, geo, varname)
    res_dir = "./res/vars_transient/"
    if !ispath(res_dir)
        mkpath(res_dir)
    end
    combined_df = DataFrame()
    wall_components = (:evap, :vl, :cond, :ll, :cc)
    current_position = 0.0
    
    for component in wall_components
        comp_data = getproperty(T_wall, component)
        N_nodes = length(comp_data)
        if component == :evap
            dl = getproperty(geo, :grv).L
        elseif component == :cc
            dl = getproperty(geo, :cc).L
        else
            dl = getproperty(geo, component).dl
        end
        positions = current_position .+ [(i - 0.5) * dl for i in 1:N_nodes]

        # 各部位の DataFrame を作成
        df = DataFrame(
            Component = fill(string(component), N_nodes),
            Position = positions,
            Node    = 1:N_nodes,
            T_wall  = comp_data
        )        
        append!(combined_df, df)

        current_position += dl * N_nodes
    end

    filename = res_dir * varname * ".csv"
    
    CSV.write(filename, combined_df)
end

function save_var_to_csv(m_evap, var, geo, varname)
    res_dir = "./res/vars_transient/"
    if !ispath(res_dir)
        mkpath(res_dir)
    end
    combined_df = DataFrame()
    fluid_components = (:vl, :cond, :ll)
    current_position = 0.0

    # wick
    comp_data = getproperty(var, :wick)
    N_nodes = length(comp_data)
    dl = getproperty(geo, :wick).dD
    positions = current_position .+ [(i - 0.5) * dl for i in 1:N_nodes]

    df = DataFrame(
        Component = fill("wick", N_nodes),
        Position = positions,
        Node = 1:N_nodes,
        T = comp_data,
        P = fill(NaN, N_nodes),
        α = fill(NaN, N_nodes),
        m = m_evap

    )
    append!(combined_df, df)

    # grv
    comp_data = getproperty(var, :grv)
    N_nodes = length(comp_data.T)
    dl = getproperty(geo, :grv).L
    positions = current_position .+ [(i - 0.5) * dl for i in 1:N_nodes]

    df = DataFrame(
        Component = fill("grv", N_nodes),
        Position = positions,
        Node = 1:N_nodes,
        T = comp_data.T,
        P = comp_data.P,
        α = comp_data.α,
        m = comp_data.m
    )
    append!(combined_df, df)
    
    for component in fluid_components
        comp_data = getproperty(var, component)
        N_nodes = length(comp_data.T)
        dl = getproperty(geo, component).dl
        positions = current_position .+ [(i - 0.5) * dl for i in 1:N_nodes]

        # 各部位の DataFrame を作成
        df = DataFrame(
            Component   = fill(string(component), N_nodes),
            Position    = positions,
            Node        = 1:N_nodes,
            T           = comp_data.T,
            P           = comp_data.P,
            α           = comp_data.α,
            m           = comp_data.m
        )        
        append!(combined_df, df)

        current_position += dl * N_nodes
    end

    # CC
    comp_data = getproperty(var, :cc)
    N_nodes = length(comp_data.T)
    dl = getproperty(geo, :cc).L
    positions = current_position .+ [(i - 0.5) * dl for i in 1:N_nodes]

    df = DataFrame(
        Component = fill("cc", N_nodes),
        Position = positions,
        Node = 1:N_nodes,
        T = comp_data.T,
        P = comp_data.P,
        α = comp_data.α,
        m = comp_data.m
    )
    append!(combined_df, df)

    filename = res_dir * varname * ".csv"
    
    CSV.write(filename, combined_df)
end

function save_f_to_csv(f, geo, varname)
    res_dir = "./res/vars_transient/"
    if !ispath(res_dir)
        mkpath(res_dir)
    end
    combined_df = DataFrame()
    fluid_components = (:vl, :cond, :ll)
    current_position = 0.0

    # wick
    comp_data = getproperty(f, :wick)
    N_nodes = length(comp_data)
    dl = getproperty(geo, :wick).dD
    positions = current_position .+ [(i - 0.5) * dl for i in 1:N_nodes]

    df = DataFrame(
        Component = fill("wick", N_nodes),
        Position = positions,
        Node = 1:N_nodes,
        f_mass = fill(NaN, N_nodes),
        f_momentum = fill(NaN, N_nodes),
        f_energy = comp_data,
        f_state = fill(NaN, N_nodes)

    )
    append!(combined_df, df)

    # grv
    comp_data = getproperty(f, :grv)
    N_nodes = length(comp_data.T)
    dl = getproperty(geo, :grv).L
    positions = current_position .+ [(i - 0.5) * dl for i in 1:N_nodes]

    df = DataFrame(
        Component = fill("grv", N_nodes),
        Position = positions,
        Node = 1:N_nodes,
        f_mass = comp_data.α,
        f_momentum = comp_data.m,
        f_energy = comp_data.T,
        f_state = comp_data.P
        # T = comp_data.T,
        # P = comp_data.P,
        # α = comp_data.α,
        # m = comp_data.m
    )
    append!(combined_df, df)
    
    for component in fluid_components
        comp_data = getproperty(f, component)
        N_nodes = length(comp_data.T)
        dl = getproperty(geo, component).dl
        positions = current_position .+ [(i - 0.5) * dl for i in 1:N_nodes]

        # 各部位の DataFrame を作成
        df = DataFrame(
            Component   = fill(string(component), N_nodes),
            Position    = positions,
            Node        = 1:N_nodes,
            f_mass = comp_data.α,
            f_momentum = comp_data.m,
            f_energy = comp_data.T,
            f_state = comp_data.P
        )        
        append!(combined_df, df)

        current_position += dl * N_nodes
    end

    # CC
    comp_data = getproperty(f, :cc)
    N_nodes = length(comp_data.T)
    dl = getproperty(geo, :cc).L
    positions = current_position .+ [(i - 0.5) * dl for i in 1:N_nodes]

    df = DataFrame(
        Component = fill("cc", N_nodes),
        Position = positions,
        Node = 1:N_nodes,
        f_mass = comp_data.α,
        f_momentum = comp_data.m,
        f_energy = comp_data.T,
        f_state = comp_data.P
    )
    append!(combined_df, df)

    filename = res_dir * varname * ".csv"
    
    CSV.write(filename, combined_df)
end

# steady
function save_steady_var_to_csv(vars, geo, condition)
    res_dir = "./res/vars_steady/"
    if !ispath(res_dir)
        mkpath(res_dir)
    end
    combined_df = DataFrame()
    
    current_position = 0 # 累積距離

    for component in keys(vars)
        comp_data = vars[component]
        N_nodes = length(comp_data.T)
        comp_geo = getproperty(geo, component)
        dl = comp_geo.dl
        
        nodes_positions = current_position .+ ((1:N_nodes) .- 0.5) .* dl
        # 各部位の DataFrame を作成
        df = DataFrame(
            Position  = nodes_positions,
            Component = fill(string(component), N_nodes),
            Node      = 1:N_nodes,
            T         = comp_data.T,
            P         = comp_data.P,
            # x（蒸気クオリティ）を持たない部位は NaN（欠損値）で埋める
            x         = haskey(comp_data, :x) ? comp_data.x : fill(NaN, N_nodes),
            α         = comp_data.α
        )
        
        # メインの DataFrame に縦結合
        append!(combined_df, df)

        # 次のコンポーネントの始点として、この部位全体の物理長さを累積距離に加算
        current_position += dl * N_nodes
    end

    filename = res_dir * "Q_apply_$(round(Int, condition.Q_apply))W.csv"
    
    CSV.write(filename, combined_df)
end

function save_residual_to_csv(f, condition)
    res_f_dir = "./res/residuals/"
    if !ispath(res_f_dir)
        mkpath(res_f_dir)
    end
    combined_df = DataFrame()

    for i in 1:2
        df = DataFrame(
            num = i,
            f     = f[i]
        )
        append!(combined_df, df)
    end

    filename = res_f_dir * "Q_apply_$(round(Int, condition.Q_apply))W_residuals.csv"
    CSV.write(filename, combined_df)
end

function save_Q_T_m_to_csv(Q_ev, Q_gr, Q_leak, T_ev, T_ec, T_cc, m)
    res_Q_T_dir = "./res/Q_Ts/"
    if !ispath(res_Q_T_dir)
        mkpath(res_Q_T_dir)
    end
    combined_df = DataFrame(Q_ev = Q_ev,
                            Q_gr = Q_gr,
                            Q_leak = Q_leak,
                            T_ev = T_ev,
                            T_ec = T_ec,
                            T_cc = T_cc,
                            m    = m
                        )
    filename = res_Q_T_dir * "Q_Ts.csv"
    CSV.write(filename, combined_df)
end