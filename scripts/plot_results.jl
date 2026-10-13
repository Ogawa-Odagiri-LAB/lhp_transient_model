# =============================================================================
# plot_results.jl (保存済みの結果を眺めるスクリプト)
#
# save_results が書き出した CSV を読んで、プロットする。計算 (transient.jl) をやり直さずに、
# いつでも、何度でも、前回の結果を見られる。LHPTransientModel (Sundials など) は読み込まない
# ので、起動も速い。
#
# 使い方:
#   julia --project plot_results.jl        … 全部の図を作って、PNG に保存する
#   REPL で  include("plot_results.jl")    … そのあと、df_prof などの変数が残るので、
#                                            任意の時刻を、その場で何度でも指定して見られる:
#     plot_distribution(df_prof, df_wick, lengths, 75.0; summary = df_summary)
#     plot_distribution(df_prof, df_wick, lengths, [10.0, 60.0, 120.0]; summary = df_summary)
#     plot_component_series(df_prof, :T; comps = ("vl", "cond"), agg = :max)
#     available_times(df_prof)           # 保存されている時刻の一覧
# =============================================================================

using Plots, DataFrames, CSV, CoolProp
include(joinpath(@__DIR__, "..", "src", "plot.jl"))

result_name = "transient_flat_ammonia"

# ---------------------------- 設定 ----------------------------
# (REPL で、値を書き換えて何度も include できるよう、const にはしていない)
# 結果のフォルダ (transient.jl の save_results の出力先)
RES_DIR = joinpath(@__DIR__, "..", "res", result_name)
# 分布を見たい時刻 [s] (保存されている時刻のうち、最も近いものを使う)
TIMES = [10.0, 50.0, 100.0, 200.0]
# 作動流体 (蒸発による熱 mdot × h_fg を計算するのに使う)
FLUID = "ammonia"
# 基準線 (壁温度の図に引く)。transient.jl の運転条件と同じ値にする
T_COND = 15.0 + 273.15      # 凝縮器プレートの温度 [K]
T_AMB = 25.0 + 273.15       # 周囲温度 [K]
# geometry.csv が無い古い結果 (最初の 200 s の実行) のときに使う、各成分の長さ [m]。
# transient.jl の幾何と同じ値にする (セル数は、CSV から数える)。
LENGTHS_FALLBACK = (grv = 0.035, vl = 2.0, cond = 1.0, ll = 2.0, cc = 0.07)
# 成分内のセルの代表値 (時系列): :mean, :min, :max, :first (最上流), :last (最下流)
AGG = :mean
# 画像を保存するか (RES_DIR/overview.png など)。REPL で試行錯誤するときは false でもよい
SAVE_FIGURES = true
# コマンドラインで実行するとき、ウィンドウが閉じないよう、Enter を待つか
WAIT_FOR_ENTER = false
# --------------------------------------------------------------

# 日本語のラベルにしたいとき (日本語フォントが必要):
# set_plot_language!(:ja; fontfamily = "Meiryo")

df_summary, df_prof, df_wick = load_results(RES_DIR)
lengths = isfile(joinpath(RES_DIR, "geometry.csv")) ? load_lengths(RES_DIR) : LENGTHS_FALLBACK
isfile(joinpath(RES_DIR, "geometry.csv")) ||
    println("geometry.csv が無いので、LENGTHS_FALLBACK を使います: ", LENGTHS_FALLBACK)

avail = available_times(df_prof)
println("結果: ", RES_DIR)
println("  保存されている時刻: ", avail[1], " 〜 ", avail[end], " s (", length(avail), " 点)")
println("  成分ごとのセル数: ", Dict(c => maximum(df_prof.cell[df_prof.comp .== c]) for c in PLOT_FLUID_COMPS))

# 古い結果は、蒸発による熱 Q_evap の列を持たないので、h_fg(T_out) を、CoolProp で計算する
h_fg = T -> CoolProp.PropsSI("H", "T", T, "Q", 1.0, FLUID) - CoolProp.PropsSI("H", "T", T, "Q", 0.0, FLUID)

# --- 時系列の一覧 (成分ごとに色・凡例) ---
plt_overview = plot_overview(df_summary, df_prof; agg = AGG, T_cond = T_COND, T_amb = T_AMB, h_fg = h_fg)

# --- 任意の時刻の、ループ全体の分布 (温度、乾き度、ウィック) ---
times = filter(t -> avail[1] - 1e-9 <= t <= avail[end] + 1e-9, TIMES)       # 範囲外の時刻は除く
isempty(times) && (times = [avail[end]])
plt_dist = plot_distribution(df_prof, df_wick, lengths, times; summary = df_summary)          # 重ね描き
plt_dist_last = plot_distribution(df_prof, df_wick, lengths, times[end]; summary = df_summary)  # 最後の時刻 (成分ごとに色分け)

if SAVE_FIGURES
    savefig(plt_overview, joinpath(RES_DIR, "overview.png"))
    savefig(plt_dist, joinpath(RES_DIR, "distribution.png"))
    savefig(plt_dist_last, joinpath(RES_DIR, "distribution_t$(round(Int, times[end])).png"))
    println("画像を保存しました: ", abspath(RES_DIR), " (overview.png, distribution.png, distribution_t*.png)")
end

display(plt_overview)
display(plt_dist)
display(plt_dist_last)

if WAIT_FOR_ENTER
    println("Enter で終了")
    readline()
end