function initialize_transient_var(geo::LHPGeometry)
    return ComponentArray(
        wick    = zeros(geo.wick.N),
        grv     = (T=0.0,               P=0.0,                  α=0.0,                  m=0.0),
        vl      = (T=zeros(geo.vl.N),   P=zeros(geo.vl.N),      α=ones(geo.vl.N),       m=zeros(geo.vl.N)),
        cond    = (T=zeros(geo.cond.N), P=zeros(geo.cond.N),    α=ones(geo.cond.N),     m=zeros(geo.cond.N)),
        ll      = (T=zeros(geo.ll.N),   P=zeros(geo.ll.N),      α=zeros(geo.ll.N),      m=zeros(geo.ll.N)),
        cc      = (T=0.0,               P=0.0,                  α=0.0,                  m=0.0)
    )
end

function set_values_on_var!(var, condition)
    FLUID = condition.fluid
    
    # T
    var.wick   .= 320.0 # 303.0
    var.wick[1] -= 0.5
    var.grv.T   = 300.0
    var.vl.T   .= 300.0
    var.cond.T .= 290.0
    var.ll.T   .= 300.0
    var.cc.T    = 300.0

    # P
    # var.wick
    var.grv.P   = calc_Psat(320.0, FLUID)
    var.vl.P   .= calc_Psat(300.0, FLUID)
    var.cond.P .= calc_Psat(288.0, FLUID)
    var.ll.P   .= calc_Psat(300.0, FLUID)
    var.cc.P    = calc_Psat(300.0, FLUID)

    # α
    # var.wick
    var.grv.α   = 1.0
    var.vl.α   .= 1.0
    var.cond.α .= 0.5
    var.ll.α   .= 0.0
    var.cc.α    = 0.3

    # m
    # var.wick
    var.grv.m   = 0.00015
    var.vl.m   .= 0.00015
    var.cond.m .= 0.00015
    var.ll.m   .= 0.00015
    var.cc.m    = 0.00015    
end

function set_smooth_values_on_var(geo, condition)
    FLUID = condition.fluid

    wick_T = fill(23.0 + 273.15, geo.wick.N)

    grv_T  = 23.0 + 273.15
    grv_P  = calc_Psat(grv_T, FLUID)
    grv_α  = 1.0
    grv_m  = 0.00015

    vl_T   = collect(range(23.0 + 273.15, 20.0 + 273.15, geo.vl.N)) 
    vl_P   = collect(range(calc_Psat(23.0 + 273.15, FLUID), calc_Psat(20.0 + 273.15, FLUID), geo.vl.N))
    vl_α   = fill(1.0 - 1.0e-4, geo.vl.N)
    vl_m   = fill(0.00015, geo.vl.N)

    cond_T = fill(20.0 + 273.15, geo.cond.N)
    cond_P = fill(calc_Psat(20.0 + 273.15, FLUID), geo.cond.N)
    cond_α = collect(range(1.0, 0.0, geo.cond.N))
    cond_m = fill(0.00015, geo.cond.N)

    ll_T   = collect(range(20.0 + 273.15, 23.0 + 273.15, geo.ll.N))
    ll_P   = collect(range(calc_Psat(20.0 + 273.15, FLUID), calc_Psat(23.0 + 273.15, FLUID), geo.ll.N))
    ll_α   = fill(0.0 + 1.0e-4, geo.ll.N)
    ll_m   = fill(0.00015, geo.ll.N)
    
    cc_T   = 23.0 + 273.15
    cc_P   = calc_Psat(cc_T, FLUID)
    cc_α   = 0.5
    cc_m  = 0.00015

    return ComponentArray(
            wick    = 
                wick_T,
            grv     = (
                T= grv_T,
                P= grv_P,
                α= grv_α,
                m= grv_m,
            ),
            vl      = (
                T= vl_T,
                P= vl_P,
                α= vl_α,
                m= vl_m,
            ),
            cond      = (
                T= cond_T,
                P= cond_P,
                α= cond_α,
                m= cond_m,
            ),
            ll      = (
                T= ll_T,
                P= ll_P,
                α= ll_α,
                m= ll_m,
            ),
            cc      = (
                T= cc_T,
                P= cc_P,
                α= cc_α,
                m= cc_m,
            ),
        )
end

function set_values_on_lb!(var)    
    # T
    var.wick   .= 273.15 # 303.0
    var.grv.T  .= 273.15
    var.vl.T   .= 273.15
    var.cond.T .= 273.15
    var.ll.T   .= 273.15
    var.cc.T   .= 273.15

    # P
    # var.wick
    var.grv.P  .= 0.0
    var.vl.P   .= 0.0
    var.cond.P .= 0.0
    var.ll.P   .= 0.0
    var.cc.P   .= 0.0

    # α
    # var.wick
    var.grv.α  .= -Inf
    var.vl.α   .= -Inf
    var.cond.α .= -Inf
    var.ll.α   .= -Inf
    var.cc.α   .= -Inf

    # m
    # var.wick
    var.grv.m  .= 0.0
    var.vl.m   .= 0.0
    var.cond.m .= 0.0
    var.ll.m   .= 0.0
    var.cc.m   .= 0.0    
end

function set_values_on_ub!(var)    
    # T
    var.wick   .= Inf # 303.0
    var.grv.T  .= Inf
    var.vl.T   .= Inf
    var.cond.T .= Inf
    var.ll.T   .= Inf
    var.cc.T   .= Inf

    # P
    # var.wick
    var.grv.P  .= Inf
    var.vl.P   .= Inf
    var.cond.P .= Inf
    var.ll.P   .= Inf
    var.cc.P   .= Inf

    # α
    # var.wick
    var.grv.α  .= Inf
    var.vl.α   .= Inf
    var.cond.α .= Inf
    var.ll.α   .= Inf
    var.cc.α   .= Inf

    # m
    # var.wick
    var.grv.m  .= Inf
    var.vl.m   .= Inf
    var.cond.m .= Inf
    var.ll.m   .= Inf
    var.cc.m   .= Inf    
end

function initialize_T_wall(geo::LHPGeometry)
    return ComponentArray(
        evap = 0.0,
        vl   = zeros(geo.vl.N),
        cond = zeros(geo.cond.N),
        ll   = zeros(geo.ll.N),
        cc   = 0.0      
    )
end

function set_values_on_T_wall!(T_wall, geo, condition)
    T_wall.evap  = 23.0 + 273.15 # 318.0
    T_wall.vl    = collect(range(23.0 + 273.15, 20.0 + 273.15, geo.vl.N))
    T_wall.cond .= 20.0 + 273.15
    T_wall.ll    = collect(range(20.0 + 273.15, 23.0 + 273.15, geo.ll.N))
    T_wall.cc    = 23.0 + 273.15
end

function initialize_f(geo::LHPGeometry)
    return ComponentArray(
        wick    = zeros(geo.wick.N),
        grv     = (T=0.0,               P=0.0,                  α=0.0,                  m=0.0),
        vl      = (T=zeros(geo.vl.N),   P=zeros(geo.vl.N),      α=zeros(geo.vl.N),      m=zeros(geo.vl.N)),
        cond    = (T=zeros(geo.cond.N), P=zeros(geo.cond.N),    α=zeros(geo.cond.N),    m=zeros(geo.cond.N)),
        ll      = (T=zeros(geo.ll.N),   P=zeros(geo.ll.N),      α=zeros(geo.ll.N),      m=zeros(geo.ll.N)),
        cc      = (T=0.0,               P=0.0,                  α=0.0,                  m=0.0)
    )
end

function λ_plotter()

    plot(
        title = "latent heat",
        xlabel = "Temp[K]",
        ylabel = "latent heat[J]",
        legend = :outerright, # 凡例をグラフの外側右に配置
        grid = true,
        size = (800, 500)     # グラフのサイズ (幅, 高さ)
    )

    T = 300.0:1.0:350.0
    λ = calc_latent_heat.(T, "propylene")

    plot!(
        T,
        λ,
        marker = :circle,
        markersize = 4,
        linewidth = 2
    )

    current()
    fig_dir = "./res/vars_transient/figs/"
    if !ispath(fig_dir)
        mkpath(fig_dir)
    end
    fig_file = fig_dir * "latent_heat" * ".png"
    savefig(fig_file)
end

function show_sat()
    @show FLUID_TABLES["Propylene"].Psat_T(300.0)
    @show calc_Psat(300.0, "Propylene")
end

function initialize_steady_vars(geo::LHPSteadyGeometry)
    return ComponentArray(
        grv     = (T=zeros(geo.grv.N),  P=zeros(geo.grv.N),     x=ones(geo.grv.N),      α=ones(geo.grv.N),  mass=zeros(geo.grv.N)),
        vl      = (T=zeros(geo.vl.N),   P=zeros(geo.vl.N),      x=ones(geo.vl.N),       α=ones(geo.vl.N),   mass=zeros(geo.vl.N)),
        cond    = (T=zeros(geo.cond.N), P=zeros(geo.cond.N),    x=ones(geo.cond.N),     α=ones(geo.cond.N), mass=zeros(geo.cond.N)),
        ll      = (T=zeros(geo.ll.N),   P=zeros(geo.ll.N),      x=zeros(geo.ll.N),      α=zeros(geo.ll.N),  mass=zeros(geo.ll.N)),
        bay     = (T=zeros(geo.bay.N),  P=zeros(geo.bay.N),     x=zeros(geo.bay.N),     α=zeros(geo.bay.N), mass=zeros(geo.bay.N)),
        cc      = (T=zeros(geo.cc.N),   P=zeros(geo.cc.N),      x=zeros(geo.cc.N),      α=zeros(geo.cc.N),  mass=zeros(geo.cc.N))
    )
end

function set_values_on_var_from_steady(geo, Q_apply)
    steady_var_file = "./res/vars_steady/" * "Q_apply_$(Q_apply)W" * ".csv" 
    df_steady_var = CSV.read(steady_var_file, DataFrame)

    Q_Ts_file = "./res/Q_Ts/" * "Q_Ts.csv"
    df_Q_Ts = CSV.read(Q_Ts_file, DataFrame)

    # Q_Tsからの値取得
    T_ec = df_Q_Ts.T_ec[1]
    T_ev = df_Q_Ts.T_ev[1]
    T_cc = df_Q_Ts.T_cc[1]
    m_st    = df_Q_Ts.m[1]

    # 4. 各コンポーネントのデータをグループ分けして辞書化（修正版）
    comp_groups = groupby(df_steady_var, :Component)

    comp_dfs = Dict(
        Symbol(k.Component) => DataFrame(sdf) 
        for (k, sdf) in pairs(comp_groups)
    )

    # データの整合性をチェックして取得するヘルパー関数
    function get_comp_data(comp_name::Symbol, field::Symbol, expected_N::Int, default_val=0.0)
        if !haskey(comp_dfs, comp_name)
            return fill(default_val, expected_N)
        end
        
        df_sub = comp_dfs[comp_name]
        actual_N = nrow(df_sub) # 定常解析の実際のノード数
        
        # パターンA: ノード数が完全に一致している場合
        if actual_N == expected_N
            return df_sub[:, field]
            
        # パターンB: 定常の方がノード数が多い場合（数個飛びでサンプリング）
        elseif actual_N > expected_N
            @warn "$comp_name: Steady nodes ($actual_N) > Transient nodes ($expected_N). Subsampling grid data."
            
            # 1番目から actual_N 番目までを、expected_N 個に等間隔分割するインデックスを計算
            # 例: actual_N=100, expected_N=10 のとき、[1, 12, 23, 34, 45, 56, 67, 78, 89, 100] のようになる

            if expected_N == 1
                # 1つだけ抜き出してサイズ1の配列（ベクトル）にする
                # ※ df_sub[1:1, field] と書くことで、スカラーではなく1要素の配列として取得できます
                return df_sub[1:1, field] 
            end
            chosen_indices = round.(Int, range(1, actual_N, length=expected_N))
            
            return df_sub[chosen_indices, field]
            
        # パターンC: 定常の方がノード数が少ない場合（やむを得ず平均値で補間）
        else
            @warn "$comp_name: Steady nodes ($actual_N) < Transient nodes ($expected_N). Falling back to mean value."
            # 前に導入した skipmissing を使って安全に平均を計算
            ave_val = mean(skipmissing(df_sub[:, field]))
            return fill(ave_val, expected_N)
        end
    end

    # 5. 各コンポーネントの生データをCSVから抽出
    grv_T  = get_comp_data(:grv,  :T, geo.grv.N)
    grv_P  = get_comp_data(:grv,  :P, geo.grv.N)
    grv_α  = get_comp_data(:grv,  :α, geo.grv.N)
    grv_m  = fill(m_st, geo.grv.N)

    vl_T   = get_comp_data(:vl,   :T, geo.vl.N)
    vl_P   = get_comp_data(:vl,   :P, geo.vl.N)
    vl_α   = get_comp_data(:vl,   :α, geo.vl.N)

    cond_T = get_comp_data(:cond, :T, geo.cond.N)
    cond_P = get_comp_data(:cond, :P, geo.cond.N)
    cond_α = get_comp_data(:cond, :α, geo.cond.N)

    ll_T   = get_comp_data(:ll,   :T, geo.ll.N)
    ll_P   = get_comp_data(:ll,   :P, geo.ll.N)
    ll_α   = get_comp_data(:ll,   :α, geo.ll.N)
    
    cc_T   = get_comp_data(:cc,   :T, geo.cc.N)
    cc_P   = get_comp_data(:cc,   :P, geo.cc.N)
    cc_α   = get_comp_data(:cc,   :α, geo.cc.N)
    cc_m   = fill(m_st, geo.cc.N)

    return ComponentArray(
        wick    = fill(T_ev, geo.wick.N),
        grv     = (
            T= grv_T[1],
            P= grv_P[1],
            α= grv_α[1],
            m= grv_m[1]
        ),
        vl      = (
            T= vl_T,   
            P= vl_P,      
            α= vl_α,      
            m= fill(m_st, geo.vl.N)
        ),
        cond      = (
            T= cond_T,   
            P= cond_P,      
            α= cond_α,      
            m= fill(m_st, geo.cond.N)
        ),
        ll      = (
            T= ll_T,   
            P= ll_P,      
            α= ll_α,      
            m= fill(m_st, geo.ll.N)
        ),
        cc      = (
            T= cc_T[1],   
            P= cc_P[1],      
            α= cc_α[1],      
            m= cc_m[1]
        ),
    )
    println("var is imported from steady analysis")
end

function set_values_on_T_wall_from_steady(geo, condition)
    Q_Ts_file = "./res/Q_Ts/" * "Q_Ts.csv"
    df_Q_Ts = CSV.read(Q_Ts_file, DataFrame)

    # Q_Tsからの値取得
    T_ec = df_Q_Ts.T_ec[1]
    T_cc = df_Q_Ts.T_cc[1]

    T_amb = condition.T_amb

    return ComponentArray(
        evap = T_ec,
        vl   = fill(T_amb, geo.vl.N),
        cond = fill(T_amb, geo.cond.N),
        ll   = fill(T_amb, geo.ll.N),
        cc   = T_cc
    )
    println("T_wall is imported from steady analysis")
end