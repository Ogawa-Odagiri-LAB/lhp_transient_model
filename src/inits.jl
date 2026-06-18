function initialize_steady_vars(geo::LHPGeometry)
    return ComponentArray(
        grv     = (T=zeros(geo.grv.N),  P=zeros(geo.grv.N),     x=ones(geo.grv.N),      α=ones(geo.grv.N),  mass=zeros(geo.grv.N)),
        vl      = (T=zeros(geo.vl.N),   P=zeros(geo.vl.N),      x=ones(geo.vl.N),       α=ones(geo.vl.N),   mass=zeros(geo.vl.N)),
        cond    = (T=zeros(geo.cond.N), P=zeros(geo.cond.N),    x=ones(geo.cond.N),     α=ones(geo.cond.N), mass=zeros(geo.cond.N)),
        ll      = (T=zeros(geo.ll.N),   P=zeros(geo.ll.N),      x=zeros(geo.ll.N),      α=zeros(geo.ll.N),  mass=zeros(geo.ll.N)),
        bay     = (T=zeros(geo.bay.N),  P=zeros(geo.bay.N),     x=zeros(geo.bay.N),     α=zeros(geo.bay.N), mass=zeros(geo.bay.N)),
        cc      = (T=zeros(geo.cc.N),   P=zeros(geo.cc.N),      x=zeros(geo.cc.N),      α=zeros(geo.cc.N),  mass=zeros(geo.cc.N))
    )
end

function initialize_transient_vars(geo::LHPGeometry)
    return ComponentArray(
        wick    = (T_in, m_wick),
        grv     = (T=zeros(geo.grv.N),  P=zeros(geo.grv.N),     α=ones(geo.grv.N),      m=zeros(geo.grv.N),     T_wall=zeros(geo.grv.N)),
        vl      = (T=zeros(geo.vl.N),   P=zeros(geo.vl.N),      α=ones(geo.vl.N),       m=zeros(geo.vl.N),      T_wall=zeros(geo.vl.N)),
        cond    = (T=zeros(geo.cond.N), P=zeros(geo.cond.N),    α=ones(geo.cond.N),     m=zeros(geo.cond.N),    T_wall=zeros(geo.cond.N)),
        ll      = (T=zeros(geo.ll.N),   P=zeros(geo.ll.N),      α=zeros(geo.ll.N),      m=zeros(geo.ll.N),      T_wall=zeros(geo.ll.N)),
        bay     = (T=zeros(geo.bay.N),  P=zeros(geo.bay.N),     α=zeros(geo.bay.N),     m=zeros(geo.bay.N),     T_wall=zeros(geo.bay.N)),
        cc      = (T_cc, α_cc, T_wall)
    )
end

function load_initial_vars_from_csv(geo::LHPGeometry, transient_condition)
    # 1. フォルダパスとファイル名の設定
    res_dir = "./res/vars/"
    res_Q_T_dir = "./res/Q_Ts/"
    
    vars_filename = res_dir * "Q_apply_$(round(Int, transient_condition.Q_apply_init))W.csv"
    qt_filename = res_Q_T_dir * "Q_Ts.csv"
    
    # ファイルの存在チェック
    if !isfile(vars_filename) || !isfile(qt_filename)
        error("targeted CSV not found \n$vars_filename\n$qt_filename")
    end

    # 2. CSVデータの読み込み
    df_vars = CSV.read(vars_filename, DataFrame)
    df_qt = CSV.read(qt_filename, DataFrame)

    # 3. 単一値 m の取得
    init_m = df_qt.m[1]

    # 4. 各コンポーネントのデータをグループ分けして辞書化（修正版）
    comp_groups = groupby(df_vars, :Component)
    
    # pairs(comp_groups) を使うことで、k (キー情報) と sdf (実際の DataFrame) が正しく取れます
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
        if nrow(df_sub) != expected_N
            @warn "$comp_name のCSVノード数($(nrow(df_sub)))とGeometryのノード数($expected_N)が一致しません。CSVのデータを使用します。"
        end
        return df_sub[:, field]
    end

    # 5. 各コンポーネントの生データをCSVから抽出
    # (wick と bay の特殊な初期化で ll や grv の値を参照するため、先に各配列を確定させます)
    ll_T   = get_comp_data(:ll,   :T, geo.ll.N)
    ll_P   = get_comp_data(:ll,   :P, geo.ll.N)
    ll_α   = get_comp_data(:ll,   :α, geo.ll.N)

    grv_T  = get_comp_data(:grv,  :T, geo.grv.N)
    grv_P  = get_comp_data(:grv,  :P, geo.grv.N)
    grv_α  = get_comp_data(:grv,  :α, geo.grv.N)

    vl_T   = get_comp_data(:vl,   :T, geo.vl.N)
    vl_P   = get_comp_data(:vl,   :P, geo.vl.N)
    vl_α   = get_comp_data(:vl,   :α, geo.vl.N)

    cond_T = get_comp_data(:cond, :T, geo.cond.N)
    cond_P = get_comp_data(:cond, :P, geo.cond.N)
    cond_α = get_comp_data(:cond, :α, geo.cond.N)

    cc_T   = get_comp_data(:cc,   :T, geo.cc.N)
    cc_P   = get_comp_data(:cc,   :P, geo.cc.N)
    cc_α   = get_comp_data(:cc,   :α, geo.cc.N)

    # 6. ComponentArray の構築と初期値代入
    return ComponentArray(
        # ウィック: 1番目をll[1]、2番目をgrv[1]の値にする
        wick = (
            T_in = ll_T[geo.ll.N], 
            m_wick = init_m
        ),
        grv = (
            T = grv_T,  
            P = grv_P,   
            α = grv_α,       
            m = fill(init_m, geo.grv.N),     
            T_wall = fill(grv_T[1], geo.grv.N)
        ),
        vl = (
            T = vl_T,   
            P = vl_P,      
            α = vl_α,       
            m = fill(init_m, geo.vl.N),      
            T_wall = fill(vl_T[1], geo.vl.N)
        ),
        cond = (
            T = cond_T, 
            P = cond_P,    
            α = cond_α,     
            m = fill(init_m, geo.cond.N),    
            T_wall = fill(cond_T[1], geo.cond.N)
        ),
        ll = (
            T = ll_T,   
            P = ll_P,      
            α = ll_α,       
            m = fill(init_m, geo.ll.N),      
            T_wall = fill(ll_T[1], geo.ll.N)
        ),
        # バイヨネット管: llの最後の値をすべてのノードに代入
        bay = (
            T = fill(ll_T[geo.ll.N], geo.bay.N),  
            P = fill(ll_P[geo.ll.N], geo.bay.N),     
            α = fill(ll_α[geo.ll.N], geo.bay.N),       
            m = fill(init_m, geo.bay.N),     
            T_wall = fill(ll_T[geo.ll.N], geo.bay.N) # T_wall も同様に ll[1] から流用
        ),
        cc = (
            T_cc = cc_T[1],   
            α_cc = cc_α[1], 
            T_wall = cc_T[1]
        )
    )
end

function prepare_f_res(geo::LHPGeometry)
    return ComponentArray(
        wick  = zeros(2),
        cc    = zeros(2),
        mass  = (grv=zeros(geo.grv.N), vl=zeros(geo.vl.N), cond=zeros(geo.cond.N), ll=zeros(geo.ll.N), bay=zeros(geo.bay.N)),
        mom   = (grv=zeros(geo.grv.N), vl=zeros(geo.vl.N), cond=zeros(geo.cond.N), ll=zeros(geo.ll.N), bay=zeros(geo.bay.N)),
        ene   = (grv=zeros(geo.grv.N), vl=zeros(geo.vl.N), cond=zeros(geo.cond.N), ll=zeros(geo.ll.N), bay=zeros(geo.bay.N)),
        state = (grv=zeros(geo.grv.N), vl=zeros(geo.vl.N), cond=zeros(geo.cond.N), ll=zeros(geo.ll.N), bay=zeros(geo.bay.N))
    )
end