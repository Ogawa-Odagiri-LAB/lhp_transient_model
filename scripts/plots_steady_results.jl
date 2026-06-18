using Plots, DataFrames, CSV

function auto_plot_physical_length_lhp(dir_path="./res/var/", geo=nothing)
    if geo === nothing
        error("'geo' not found")
    end

    gr() 
    
    # 1. フォルダ内のCSVファイルを自動検出
    all_files = readdir(dir_path)
    csv_files = filter(f -> endswith(f, ".csv"), all_files)
    
    if isempty(csv_files)
        println("no CSV found")
        return nothing
    end
    
    # 【追加】保存先となる png フォルダのパスを設定し、自動作成する
    output_dir = joinpath(dir_path, "png")
    mkpath(output_dir) # フォルダが存在しない場合は自動作成、既にある場合は何もしない
    
    # LHPの物理的な流路の順番を定義
    flow_order = ["grv", "vl", "cond", "ll"]
    
    # geoオブジェクトから各部位の全長さLとノード数Nをマッピング
    # ※最初の定義（GrvGeometryにNやLがない場合）に備え、プログラム側のフィールド名（N_grv, L_grvなど）に合わせて適宜調整してください
    geo_map = Dict(
        "grv"  => (L = geo.grv.L,  N = geo.grv.N),
        "vl"   => (L = geo.vl.L,   N = geo.vl.N),
        "cond" => (L = geo.cond.L, N = geo.cond.N),
        "ll"   => (L = geo.ll.L,   N = geo.ll.N),
        # "bay"  => (L = geo.bay.L,  N = geo.bay.N)
    )
    
    println("start plotting")
    
    for file in csv_files
        full_path = joinpath(dir_path, file)
        df = CSV.read(full_path, DataFrame)
        
        if !all(col -> col in names(df), ["Component", "Node", "T", "P", "x"])
            continue
        end
        
        ordered_df = DataFrame()
        current_position = 0.0  # LHPの始点からの累積距離 [m]
        
        # 境界線の位置情報を記録する配列
        boundaries = Float64[]
        boundary_labels = String[]
        
        for comp in flow_order
            df_comp = filter(row -> row.Component == comp, df)
            if isempty(df_comp) continue end
            
            sort!(df_comp, :Node)
            
            # 対応するコンポーネントの幾何情報を取得
            comp_geo = geo_map[comp]
            dl = comp_geo.L / comp_geo.N  # 1ノードあたりの長さ [m]
            
            # 実際の物理位置を計算
            df_comp[!, :Position] = current_position .+ ((df_comp.Node .- 0.5) .* dl)
            
            # このコンポーネントの総長さを累積位置に足す
            current_position += comp_geo.L
            
            append!(ordered_df, df_comp)
            
            # 部位の境界位置 [m] を記録
            push!(boundaries, current_position)
            push!(boundary_labels, comp)
        end
        
        # グラフ全体の基本設定
        plot_opts = (
            xlabel = "LHP Flow Path Length [m]", 
            linewidth = 2, 
            legend = false, 
            xlims = (0.0, current_position)
        )
        
        cond_title = replace(file, ".csv" => "")
        
        # X軸を物理位置にして描画
        p_T = plot(ordered_df.Position, ordered_df.T, title="Temperature Profile ($(cond_title))", ylabel="Temperature [K]"; plot_opts...)
        p_P = plot(ordered_df.Position, ordered_df.P/10^6, title="Pressure Profile ($(cond_title))", ylabel="Pressure [MPa]"; plot_opts...)
        p_x = plot(ordered_df.Position, ordered_df.x, title="Vapor Quality Profile ($(cond_title))", ylabel="Vapor Quality [-]", ylims=(-0.05, 1.05); plot_opts...)
        
        # 各部位の境界に縦の点線を追加
        for i in 1:(length(boundaries)-1)
            pos = boundaries[i]
            vline!(p_T, [pos], line=(:gray, :dash, 1))
            vline!(p_P, [pos], line=(:gray, :dash, 1))
            vline!(p_x, [pos], line=(:gray, :dash, 1))
        end
        
        # 各エリアの中央にコンポーネント名のテキストを配置
        prev_b = 0.0
        for i in 1:length(boundaries)
            b = boundaries[i]
            mid_point = (prev_b + b) / 2.0
            
            y_max = maximum(ordered_df.T)
            y_min = minimum(ordered_df.T)
            annotate!(p_T, mid_point, y_max - (y_max-y_min)*0.08, text(boundary_labels[i], :center, 9, :blue, :bold))
            
            prev_b = b
        end
        
        # 縦に3段に結合
        final_plot = plot(p_T, p_P, p_x, layout=(3, 1), size=(800, 950), margin=5Plots.mm)
        
        # 【変更】保存先を新しく作った `output_dir` (dir_path/png/) の中に指定
        output_png = joinpath(output_dir, replace(file, ".csv" => ".png"))
        savefig(final_plot, output_png)
        println("images saved: $(output_png)")
    end
    println("all CSVs plotted")
end

# 改良版
function auto_plot_physical_length_lhp_2(dir_path="./res/vars/")
    gr() 
    
    # 1. フォルダ内のCSVファイルを自動検出
    all_files = readdir(dir_path)
    csv_files = filter(f -> endswith(f, ".csv"), all_files)
    
    if isempty(csv_files)
        println("no CSV found")
        return nothing
    end
    
    # 保存先となる png フォルダの自動作成
    output_dir = joinpath(dir_path, "png")
    mkpath(output_dir)
    
    # 縦の境界線を描画するために必要なコンポーネント順（流路順）
    flow_order = ["grv", "vl", "cond", "ll"]
    
    println("start plotting")
    
    for file in csv_files
        full_path = joinpath(dir_path, file)
        df = CSV.read(full_path, DataFrame)
        
        # 必須列のチェック（Position 列が最左列にあることを前提にチェック対象に追加）
        if !all(col -> col in names(df), ["Position", "Component", "Node", "T", "P", "x"])
            continue
        end
        
        # 流路の順序通りにデータをソートして並び替える
        ordered_df = DataFrame()
        boundaries = Float64[]
        boundary_labels = String[]
        
        for comp in flow_order
            df_comp = filter(row -> row.Component == comp, df)
            if isempty(df_comp) continue end
            
            sort!(df_comp, :Node)
            append!(ordered_df, df_comp)
            
            # 【改良】各部位の末尾ノードのPositionを、そのまま境界位置として記録
            # 1ノードの中心が Position = x だとすると、コンポーネントの末尾の物理端は 
            # 「中心位置 + 0.5*dl」となるため、それを正確に割り出します
            dl = (df_comp.Node[end] > 1) ? (df_comp.Position[2] - df_comp.Position[1]) : (df_comp.Position[1] * 2)
            push!(boundaries, df_comp.Position[end] + 0.5 * dl)
            push!(boundary_labels, comp)
        end
        
        # グラフ全体の最大長さを取得
        total_length = maximum(ordered_df.Position)
        
        # グラフ全体の基本設定
        plot_opts = (
            xlabel = "LHP Flow Path Length [m]", 
            linewidth = 2, 
            legend = false, 
            xlims = (0.0, total_length)
        )
        
        cond_title = replace(file, ".csv" => "")
        
        # 【改良】CSVから読み取った Position 列をそのままX軸（第1引数）に指定
        p_T = plot(ordered_df.Position, ordered_df.T, title="Temperature Profile ($(cond_title))", ylabel="Temperature [K]"; plot_opts...)
        p_P = plot(ordered_df.Position, ordered_df.P/10^6, title="Pressure Profile ($(cond_title))", ylabel="Pressure [MPa]"; plot_opts...)
        p_x = plot(ordered_df.Position, ordered_df.x, title="Vapor Quality Profile ($(cond_title))", ylabel="Vapor Quality [-]", ylims=(-0.05, 1.05); plot_opts...)
        
        # 各部位の境界に縦の点線を追加
        for i in 1:(length(boundaries)-1)
            pos = boundaries[i]
            vline!(p_T, [pos], line=(:gray, :dash, 1))
            vline!(p_P, [pos], line=(:gray, :dash, 1))
            vline!(p_x, [pos], line=(:gray, :dash, 1))
        end
        
        # 各エリアの中央にコンポーネント名のテキストを配置
        prev_b = 0.0
        for i in 1:length(boundaries)
            b = boundaries[i]
            mid_point = (prev_b + b) / 2.0
            
            y_max = maximum(ordered_df.T)
            y_min = minimum(ordered_df.T)
            annotate!(p_T, mid_point, y_max - (y_max-y_min)*0.08, text(boundary_labels[i], :center, 9, :blue, :bold))
            
            prev_b = b
        end
        
        # 縦に3段に結合
        final_plot = plot(p_T, p_P, p_x, layout=(3, 1), size=(800, 950), margin=5Plots.mm)
        
        output_png = joinpath(output_dir, replace(file, ".csv" => ".png"))
        savefig(final_plot, output_png)
        println("images saved: $(output_png)")
    end
    println("all CSVs plotted")
end

auto_plot_physical_length_lhp_2("./res/vars/")