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

function save_transient_vars_to_csv(vars, T_wall, geo, condition, t)
    res_dir = "./res/vars_transient/"
    if !ispath(res_dir)
        mkpath(res_dir)
    end
    combined_df = DataFrame()
    
    current_position = 0.0 # 累積距離

    # 1. 通常の流路コンポーネントの処理
    flow_components = (:grv, :vl, :cond, :ll, :bay)
    
    for component in flow_components
        comp_data = getproperty(vars, component)
        comp_wall = getproperty(T_wall, component)
        
        N_nodes = length(comp_data.T)
        comp_geo = getproperty(geo, component)
        dl = comp_geo.dl
        
        # ノードの中心位置を計算
        nodes_positions = current_position .+ ((1:N_nodes) .- 0.5) .* dl
        
        # 各部位の DataFrame を作成
        df = DataFrame(
            Time      = fill(t, N_nodes),             # 時刻を記録
            Component = fill(string(component), N_nodes),
            Node      = 1:N_nodes,
            Position  = nodes_positions,
            T         = comp_data.T,
            P         = comp_data.P,
            α         = comp_data.α,
            m         = comp_data.m,
            T_wall    = comp_wall,
            T_in      = fill(NaN, N_nodes),
            m_wick    = fill(NaN, N_nodes)
        )
        
        append!(combined_df, df)

        # 次のコンポーネントの始点として累積距離を加算
        current_position += dl * N_nodes
    end

    # 2. 特殊コンポーネント (wick) の処理 (1ノードとして扱う)
    df_wick = DataFrame(
        Time      = [t],
        Component = ["wick"],
        Node      = [1],
        Position  = [NaN],
        T         = [NaN],
        P         = [NaN],
        α         = [NaN],
        m         = [NaN],
        T_wall    = [NaN],
        T_in      = [vars.wick.T_in],
        m_wick    = [vars.wick.m_wick]
    )
    append!(combined_df, df_wick)

    # 3. 特殊コンポーネント (cc) の処理 (1ノードとして扱う)
    df_cc = DataFrame(
        Time      = [t],
        Component = ["cc"],
        Node      = [1],
        Position  = [NaN],
        T         = [vars.cc.T_cc],
        P         = [NaN],
        α         = [vars.cc.α_cc],
        m         = [NaN],
        T_wall    = [T_wall.cc],
        T_in      = [NaN],
        m_wick    = [NaN]
    )
    append!(combined_df, df_cc)

    # 出力ファイル名に電力と時刻（小数点以下2桁の固定長など）を組み込む
    # 例: Q_apply_100W_t_00010.00s.csv
    time_str = lpad(string(round(t, digits=2)), 8, '0')
    filename = res_dir * "Q_apply_$(round(Int, condition.Q_apply))W_t_$(time_str)s.csv"
    
    CSV.write(filename, combined_df)
end

# ==========================================
# 1. 各物理量のスケーリング基準値（正規化因子）
# ==========================================
const SCALE_T = 100.0     # 温度[K]のスケール (300K -> 3.0)
const SCALE_P = 1.0e5     # 圧力[Pa]のスケール (100,000Pa -> 1.0)
const SCALE_M = 1.0e-3    # 質量流量[kg/s]のスケール (0.001kg/s -> 1.0)
const SCALE_A = 1.0       # ボイド率[-]はもともと0~1なので1.0

# ＝方程式ごとの残差の正規化用＝
const SCALE_RES_MASS = 1.0e-3  # 質量保存の残差スケール
const SCALE_RES_MOM  = 1.0     # 運動量保存の残差スケール
const SCALE_RES_ENE  = 1.0e3   # エネルギー保存の残差スケール (Joule単位は大きくなりやすいため)
const SCALE_RES_STAT = 1.0     # 状態方程式の残差スケール


# ==========================================
# 2. 構造体 var と ソルバー用ベクトル u の相互変換
# ==========================================

# 物理構造体 var から 正規化された1次元配列 u0 を作成する関数
function pack_variables(var, geo)
    u = Float64[]
    
    # 1. wick
    push!(u, var.wick.T_in / SCALE_T)
    push!(u, var.wick.m_wick / SCALE_M)

    # 2. grv (各ノードごと)
    for i in 1:geo.grv.N
        push!(u, var.grv.T[i] / SCALE_T)
        push!(u, var.grv.P[i] / SCALE_P)
        push!(u, var.grv.α[i] / SCALE_A)
        push!(u, var.grv.m[i] / SCALE_M)
    end
    
    # 3. vl (各ノードごと)
    for i in 1:geo.vl.N
        push!(u, var.vl.T[i] / SCALE_T)
        push!(u, var.vl.P[i] / SCALE_P)
        push!(u, var.vl.α[i] / SCALE_A)
        push!(u, var.vl.m[i] / SCALE_M)
    end

    # 4. cond (各ノードごと)
    for i in 1:geo.cond.N
        push!(u, var.cond.T[i] / SCALE_T)
        push!(u, var.cond.P[i] / SCALE_P)
        push!(u, var.cond.α[i] / SCALE_A)
        push!(u, var.cond.m[i] / SCALE_M)
    end

    # 5. ll (各ノードごと)
    for i in 1:geo.ll.N
        push!(u, var.ll.T[i] / SCALE_T)
        push!(u, var.ll.P[i] / SCALE_P)
        push!(u, var.ll.α[i] / SCALE_A)
        push!(u, var.ll.m[i] / SCALE_M)
    end
    
    # 6. bay (各ノードごと)
    for i in 1:geo.bay.N
        push!(u, var.bay.T[i] / SCALE_T)
        push!(u, var.bay.P[i] / SCALE_P)
        push!(u, var.bay.α[i] / SCALE_A)
        push!(u, var.bay.m[i] / SCALE_M)
    end
    
    # 7. cc
    push!(u, var.cc.T_cc / SCALE_T)
    push!(u, var.cc.α_cc / SCALE_A)
    
    
    return u
end

# 正規化されたベクトル u から 物理構造体 var へ値を復元（書き込み）する関数
function unpack_variables!(var, u, geo)
    idx = 1
    
    # 1. wick
    var.wick.T_in   = u[idx] * SCALE_T; idx += 1
    var.wick.m_wick = u[idx] * SCALE_M; idx += 1

    # 2. grv
    for i in 1:geo.grv.N
        var.grv.T[i] = u[idx] * SCALE_T; idx += 1
        var.grv.P[i] = u[idx] * SCALE_P; idx += 1
        var.grv.α[i] = u[idx] * SCALE_A; idx += 1
        var.grv.m[i] = u[idx] * SCALE_M; idx += 1
    end

    # 3. vl
    for i in 1:geo.vl.N
        var.vl.T[i] = u[idx] * SCALE_T; idx += 1
        var.vl.P[i] = u[idx] * SCALE_P; idx += 1
        var.vl.α[i] = u[idx] * SCALE_A; idx += 1
        var.vl.m[i] = u[idx] * SCALE_M; idx += 1
    end
    
    # 4. cond
    for i in 1:geo.cond.N
        var.cond.T[i] = u[idx] * SCALE_T; idx += 1
        var.cond.P[i] = u[idx] * SCALE_P; idx += 1
        var.cond.α[i] = u[idx] * SCALE_A; idx += 1
        var.cond.m[i] = u[idx] * SCALE_M; idx += 1
    end

    # 5. ll
    for i in 1:geo.ll.N
        var.ll.T[i] = u[idx] * SCALE_T; idx += 1
        var.ll.P[i] = u[idx] * SCALE_P; idx += 1
        var.ll.α[i] = u[idx] * SCALE_A; idx += 1
        var.ll.m[i] = u[idx] * SCALE_M; idx += 1
    end
    
    # 6. bay
    for i in 1:geo.bay.N
        var.bay.T[i] = u[idx] * SCALE_T; idx += 1
        var.bay.P[i] = u[idx] * SCALE_P; idx += 1
        var.bay.α[i] = u[idx] * SCALE_A; idx += 1
        var.bay.m[i] = u[idx] * SCALE_M; idx += 1
    end
    
    # 7. cc
    var.cc.T_cc     = u[idx] * SCALE_T; idx += 1
    var.cc.α_cc     = u[idx] * SCALE_A; idx += 1
    
end

# ==========================================
# 3. ソルバー用の正規化された残差計算関数
# ==========================================
function normalized_lhp_residual!(f_scaled, u_scaled, p)
    # p の中身を展開 (事前に NonlinearProblem の params に渡しておく)
    var         = p.var
    var_old     = p.vars_old
    geo         = p.geo
    condition   = p.transient_condition

    # ① ソルバーが振ってきた「無次元値 u_scaled」を、本来の「物理単位」に直して var にセット
    unpack_variables!(var, u_scaled, geo)

    # ② 各コンポーネント用の「物理単位の残差配列」をローカルで確保
    f_wick  = zeros(eltype(u_scaled), 2)

    f_grv_mass  = zeros(eltype(u_scaled), geo.grv.N)
    f_grv_mom   = zeros(eltype(u_scaled), geo.grv.N)
    f_grv_ene   = zeros(eltype(u_scaled), geo.grv.N)
    f_grv_state = zeros(eltype(u_scaled), geo.grv.N)

    f_vl_mass   = zeros(eltype(u_scaled), geo.vl.N)
    f_vl_mom    = zeros(eltype(u_scaled), geo.vl.N)
    f_vl_ene    = zeros(eltype(u_scaled), geo.vl.N)
    f_vl_state  = zeros(eltype(u_scaled), geo.vl.N)

    f_cond_mass   = zeros(eltype(u_scaled), geo.cond.N)
    f_cond_mom    = zeros(eltype(u_scaled), geo.cond.N)
    f_cond_ene    = zeros(eltype(u_scaled), geo.cond.N)
    f_cond_state  = zeros(eltype(u_scaled), geo.cond.N)

    f_ll_mass   = zeros(eltype(u_scaled), geo.ll.N)
    f_ll_mom    = zeros(eltype(u_scaled), geo.ll.N)
    f_ll_ene    = zeros(eltype(u_scaled), geo.ll.N)
    f_ll_state  = zeros(eltype(u_scaled), geo.ll.N)
        
    f_bay_mass  = zeros(eltype(u_scaled), geo.bay.N)
    f_bay_mom   = zeros(eltype(u_scaled), geo.bay.N)
    f_bay_ene   = zeros(eltype(u_scaled), geo.bay.N)
    f_bay_state = zeros(eltype(u_scaled), geo.bay.N)

    f_cc    = zeros(eltype(u_scaled), 2)

    # ③ 既存の物理計算関数をそのまま実行（アップロードいただいた関数群です）
    calc_transient_wick!(f_wick, var, var_old, geo, condition)
    calc_transient_grv!(f_grv_mass, f_grv_mom, f_grv_ene, f_grv_state, var, var_old, geo, condition)
    calc_transient_vl!(f_vl_mass, f_vl_mom, f_vl_ene, f_vl_state, var, var_old, geo, condition)
    calc_transient_cond!(f_cond_mass, f_cond_mom, f_cond_ene, f_cond_state, var, var_old, geo, condition)
    calc_transient_ll!(f_ll_mass, f_ll_mom, f_ll_ene, f_ll_state, var, var_old, geo, condition)
    calc_transient_bay!(f_bay_mass, f_bay_mom, f_bay_ene, f_bay_state, var, var_old, geo, condition)
    calc_transient_cc!(f_cc, var, var_old, geo, condition)

    # ④ 得られた物理的な残差も、オーダーを揃えるために正規化因子で割って f_scaled に集約
    f_idx = 1
    
    # wick (f[1]はm_wick用[mass系], f[2]はT_in用[energy系])
    f_scaled[f_idx] = f_wick[1] / SCALE_RES_MASS; f_idx += 1
    f_scaled[f_idx] = f_wick[2] / SCALE_RES_ENE;  f_idx += 1

    # grv
    for i in 1:geo.grv.N
        f_scaled[f_idx] = f_grv_mass[i]  / SCALE_RES_MASS;  f_idx += 1
        f_scaled[f_idx] = f_grv_mom[i]   / SCALE_RES_MOM;   f_idx += 1
        f_scaled[f_idx] = f_grv_ene[i]   / SCALE_RES_ENE;   f_idx += 1
        f_scaled[f_idx] = f_grv_state[i] / SCALE_RES_STAT;  f_idx += 1
    end
    
    # vl
    for i in 1:geo.vl.N
        f_scaled[f_idx] = f_vl_mass[i]   / SCALE_RES_MASS;  f_idx += 1
        f_scaled[f_idx] = f_vl_mom[i]    / SCALE_RES_MOM;   f_idx += 1
        f_scaled[f_idx] = f_vl_ene[i]    / SCALE_RES_ENE;   f_idx += 1
        f_scaled[f_idx] = f_vl_state[i]  / SCALE_RES_STAT;  f_idx += 1
    end

    # cond
    for i in 1:geo.cond.N
        f_scaled[f_idx] = f_cond_mass[i]   / SCALE_RES_MASS;  f_idx += 1
        f_scaled[f_idx] = f_cond_mom[i]    / SCALE_RES_MOM;   f_idx += 1
        f_scaled[f_idx] = f_cond_ene[i]    / SCALE_RES_ENE;   f_idx += 1
        f_scaled[f_idx] = f_cond_state[i]  / SCALE_RES_STAT;  f_idx += 1
    end

    # ll
    for i in 1:geo.ll.N
        f_scaled[f_idx] = f_ll_mass[i]   / SCALE_RES_MASS;  f_idx += 1
        f_scaled[f_idx] = f_ll_mom[i]    / SCALE_RES_MOM;   f_idx += 1
        f_scaled[f_idx] = f_ll_ene[i]    / SCALE_RES_ENE;   f_idx += 1
        f_scaled[f_idx] = f_ll_state[i]  / SCALE_RES_STAT;  f_idx += 1
    end

    # bay
    for i in 1:geo.bay.N
        f_scaled[f_idx] = f_bay_mass[i]  / SCALE_RES_MASS;  f_idx += 1
        f_scaled[f_idx] = f_bay_mom[i]   / SCALE_RES_MOM;   f_idx += 1
        f_scaled[f_idx] = f_bay_ene[i]   / SCALE_RES_ENE;   f_idx += 1
        f_scaled[f_idx] = f_bay_state[i] / SCALE_RES_STAT;  f_idx += 1
    end
    
    # cc (f[1]はT_cc用[energy系], f[2]はα_cc用[mass系])
    f_scaled[f_idx] = f_cc[1] / SCALE_RES_ENE;  f_idx += 1
    f_scaled[f_idx] = f_cc[2] / SCALE_RES_MASS; f_idx += 1
    
    return nothing
end