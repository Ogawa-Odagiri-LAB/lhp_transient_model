function save_vars_to_csv(vars, geo, condition)
    res_dir = "./res/vars/"
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