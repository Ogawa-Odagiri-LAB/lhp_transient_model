using Plots, DataFrames, CSV

function plot_T_wall(varname)
    res_file = "./res/vars_transient/" * varname * ".csv" 
    df_plot = CSV.read(res_file, DataFrame)

    # 2. グラフの土台を初期化（背景色やフォントサイズなどの設定）
    plot(
        title = "Wall Temperature Distribution",
        xlabel = "Position from Component Inlet [m]",
        ylabel = "Wall Temperature [K]",
        legend = :outerright, # 凡例をグラフの外側右に配置
        grid = true,
        size = (800, 500)     # グラフのサイズ (幅, 高さ)
    )

    # 3. コンポーネントごとにデータを抽出して自動で重ね書きする
    # CSVに保存されている順番 (:evap, :vl, :cond, :ll, :cc) にループ
    wall_components = ["evap", "vl", "cond", "ll", "cc"]

    for component in wall_components
        # 現在のコンポーネントの行だけをフィルタリング
        df_sub = filter(row -> row.Component == component, df_plot)
        
        # データが空でなければプロットに追加
        if !isempty(df_sub)
            plot!(
                df_sub.Position,  # 横軸：位置 [m]
                df_sub.T_wall,    # 縦軸：壁面温度 [K]
                label = component,     # 凡例に表示する名前
                marker = :circle, # ノードの位置に丸いマーカーをつける
                markersize = 4,
                linewidth = 2     # 線を少し太くする
            )
        end
    end

    # 4. 💡 グラフを画面に表示、または画像として保存
    fig_dir = "./res/vars_transient/figs/"
    if !ispath(fig_dir)
        mkpath(fig_dir)
    end
    fig_file = fig_dir * varname * "_profile.png"
    savefig(fig_file)
end

function plot_var(varname)
    res_file = "./res/vars_transient/" * varname * "_1.0s" * ".csv" 
    df_plot = CSV.read(res_file, DataFrame)

    # Temperature
    plot(
        title = "Independent Variables Plotted Along Distance",
        xlabel = "Position from Component Inlet [m]",
        ylabel = "Temperature [K]",
        legend = :outerright, # 凡例をグラフの外側右に配置
        grid = true,
        size = (800, 500)     # グラフのサイズ (幅, 高さ)
    )

    # 3. コンポーネントごとにデータを抽出して自動で重ね書きする
    # CSVに保存されている順番 (:evap, :vl, :cond, :ll, :cc) にループ
    wall_components = ["wick", "grv", "vl", "cond", "ll", "cc"]

    for component in wall_components
        # 現在のコンポーネントの行だけをフィルタリング
        df_sub = filter(row -> row.Component == component, df_plot)
        
        # データが空でなければプロットに追加
        if !isempty(df_sub)
            plot!(
                df_sub.Position,  # 横軸：位置 [m]
                df_sub.T,    # 縦軸：壁面温度 [K]
                label = component,     # 凡例に表示する名前
                marker = :circle, # ノードの位置に丸いマーカーをつける
                markersize = 4,
                linewidth = 2     # 線を少し太くする
            )
        end
    end

    # 4. グラフを画面に表示、または画像として保存
    fig_dir = "./res/vars_transient/figs/"
    if !ispath(fig_dir)
        mkpath(fig_dir)
    end
    fig_file = fig_dir * varname * "_profile.png"
    savefig(fig_file)
end

using Plots
using CSV
using DataFrames

function plot_transient_component_nodes(varname, dir_path="./res/vars_transient/")
    # 1. フォルダ内の「var_数字s.csv」にマッチするファイルをすべて探索
    if !ispath(dir_path)
        error("ディレクトリ $dir_path が存在しません。")
    end
    all_files = readdir(dir_path)
    csv_files = filter(f -> occursin(r"^var_\d+(\.\d+)?s\.csv$", f), all_files)
    
    if isempty(csv_files)
        error("該当するCSVファイル（var_*s.csv）が見つかりません。")
    end

    # 2. ファイル名から時間（秒）を抽出して、時系列順にソート
    parse_time(f) = parse(Float64, match(r"^var_\d+(\.\d+)?s\.csv$", f).captures[1])
    file_time_pairs = [(file=f, time=parse_time(f)) for f in csv_files]
    sort!(file_time_pairs, by = x -> x.time)

    times = [x.time for x in file_time_pairs]

    # 時系列データを格納するコンテナ（各コンポーネントのNode=1用）
    t_wick = Float64[]
    t_grv  = Float64[]
    t_vl   = Float64[]
    t_cond = Float64[]
    t_ll   = Float64[]
    t_cc   = Float64[]

    # 3. 各ファイルをループして Node==1 の温度データを抽出
    for pair in file_time_pairs
        file_full_path = joinpath(dir_path, pair.file)
        df = CSV.read(file_full_path, DataFrame)
        
        # 各コンポーネント名、かつ Node == 1 の行をフィルタリングして温度 T を取得
        # ※確実に値を取得するため、該当行の最初の1要素（[1]）を抽出
        push!(t_wick, df[(df.Component .== "wick") .& (df.Node .== 1), :T][1])
        push!(t_grv,  df[(df.Component .== "grv")  .& (df.Node .== 1), :T][1])
        push!(t_vl,   df[(df.Component .== "vl")   .& (df.Node .== 1), :T][1])
        push!(t_cond, df[(df.Component .== "cond") .& (df.Node .== 1), :T][1])
        push!(t_ll,   df[(df.Component .== "ll")   .& (df.Node .== 1), :T][1])
        push!(t_cc,   df[(df.Component .== "cc")   .& (df.Node .== 1), :T][1])
    end

    # 4. 1つのウィンドウに重ねてプロット
    gr() # GRバックエンドを使用
    
    plot(title = "Temperature",
         xlabel = "Time [s]",
         ylabel = "Temperature [K]",
         lw = 2,
         grid = true,
         legend = :outerright) # 凡例を外側右に配置

    # 各コンポーネントのデータをプロットに追加
    plot!(times, t_wick, label="Wick", marker=:circle)
    plot!(times, t_grv,  label="grv", marker=:utriangle)
    plot!(times, t_vl,   label="vl", marker=:square)
    plot!(times, t_cond, label="cond", marker=:triangle)
    plot!(times, t_ll,   label="ll", marker=:diamond)
    plot!(times, t_cc,   label="cc", marker=:star7)

    # 5. グラフの表示と保存
    display(current())

    fig_dir = "./res/vars_transient/figs/"
    if !ispath(fig_dir)
        mkpath(fig_dir)
    end
    fig_file = fig_dir * varname * "_transient_profile.png"
    savefig(fig_file)
end

function plot_transient_T_wall_nodes(varname, dir_path="./res/vars_transient/")
    # 1. フォルダ内の「var_数字s.csv」にマッチするファイルをすべて探索
    if !ispath(dir_path)
        error("ディレクトリ $dir_path が存在しません。")
    end
    all_files = readdir(dir_path)
    csv_files = filter(f -> occursin(r"^T_wall_\d+(\.\d+)?s\.csv$", f), all_files)
    
    if isempty(csv_files)
        error("該当するCSVファイル（T_wall_*s.csv）が見つかりません。")
    end

    # 2. ファイル名から時間（秒）を抽出して、時系列順にソート
    parse_time(f) = parse(Float64, match(r"^T_wall_\d+(\.\d+)?s\.csv$", f).captures[1])
    file_time_pairs = [(file=f, time=parse_time(f)) for f in csv_files]
    sort!(file_time_pairs, by = x -> x.time)

    times = [x.time for x in file_time_pairs]

    # 時系列データを格納するコンテナ（各コンポーネントのNode=1用）
    t_evap = Float64[]
    t_vl   = Float64[]
    t_cond = Float64[]
    t_ll   = Float64[]
    t_cc   = Float64[]

    # 3. 各ファイルをループして Node==1 の温度データを抽出
    for pair in file_time_pairs
        file_full_path = joinpath(dir_path, pair.file)
        df = CSV.read(file_full_path, DataFrame)
        
        # 💡 各コンポーネント名、かつ Node == 1 の行をフィルタリングして温度 T を取得
        # ※確実に値を取得するため、該当行の最初の1要素（[1]）を抽出
        push!(t_evap, df[(df.Component .== "evap") .& (df.Node .== 1), :T_wall][1])
        push!(t_vl,   df[(df.Component .== "vl")   .& (df.Node .== 1), :T_wall][1])
        push!(t_cond, df[(df.Component .== "cond") .& (df.Node .== 1), :T_wall][1])
        push!(t_ll,   df[(df.Component .== "ll")   .& (df.Node .== 1), :T_wall][1])
        push!(t_cc,   df[(df.Component .== "cc")   .& (df.Node .== 1), :T_wall][1])
    end

    # 4. 1つのウィンドウに重ねてプロット
    gr() # GRバックエンドを使用
    
    plot(title = "Temperature",
         xlabel = "Time [s]",
         ylabel = "Temperature [K]",
         lw = 2,
         grid = true,
         legend = :outerright) # 凡例を外側右に配置

    # 各コンポーネントのデータをプロットに追加
    plot!(times, t_evap, label="evap", marker=:circle)
    plot!(times, t_vl,   label="vl", marker=:square)
    plot!(times, t_cond, label="cond", marker=:triangle)
    plot!(times, t_ll,   label="ll", marker=:diamond)
    plot!(times, t_cc,   label="cc", marker=:star7)

    # 5. グラフの表示と保存
    display(current())
    fig_dir = "./res/vars_transient/figs/"
    if !ispath(fig_dir)
        mkpath(fig_dir)
    end
    fig_file = fig_dir * varname * "_transient_profile.png"
    savefig(fig_file)
end

plot_transient_component_nodes("var")
# plot_transient_T_wall_nodes("T_wall")
# plot_T_wall("T_wall")
# plot_T_wall("T_wall_old")
plot_var("var")