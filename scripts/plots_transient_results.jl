function plot_node_time_series(res_dir::String, component_name::String, node_idx::Int)
    # 1. ディレクトリ内のCSVファイルをすべて取得
    all_files = readdir(res_dir, join=true)
    csv_files = filter(f -> endswith(f, ".csv"), all_files)
    
    if isempty(csv_files)
        error("指定されたディレクトリにCSVファイルが見つかりません: $res_dir")
    end
    
    # データの格納用ベクトル
    times = Float64[]
    temperatures = Float64[]
    
    # 2. 各CSVから該当ノードのデータを抽出
    for file in csv_files
        df = CSV.read(file, DataFrame)
        
        # 条件に合う行をフィルター (Component と Node が一致するもの)
        row = df[(df.Component .== component_name) .& (df.Node .== node_idx), :]
        
        if !isempty(row)
            push!(times, row.Time[1])
            push!(temperatures, row.T[1])
        end
    end
    
    # 3. 念のため時間（Time）の昇順でソート（ファイル読み込み順が前後しても大丈夫なように）
    p = sortperm(times)
    times = times[p]
    temperatures = temperatures[p]
    
    # 4. プロットの作成
    plot(
        times, 
        temperatures, 
        marker = :circle, 
        linewidth = 2,
        title = "Time Series: $component_name (Node $node_idx)",
        xlabel = "Time [s]",
        ylabel = "Temperature [K]", # 必要に応じて単位（℃など）に変えてください
        label = "Node $node_idx",
        legend = :best,
        grid = true
    )
end

# 保存先ディレクトリを指定
res_dir = "./res/vars_transient/"

# grv の ノード1 の温度履歴をプロット
p = plot_node_time_series(res_dir, "grv", 1)