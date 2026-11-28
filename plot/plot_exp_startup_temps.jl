# =============================================================================
# plot_exp_startup_temps.jl  (実験データだけのプロット: 実験開始から加熱開始 + 1 分までの温度)
#
# 実験 CSV の EVA1 (= EC1)、RES1・RES2、VL1〜3、CON1〜5 を、1 枚の図に重ねて描く。
# 加熱開始の時刻は Heat Load (Q) から自動で読み取る (plot_startup_compare_heattop.jl と同じ方法)。
#
# 実験 CSV は次のどちらでも読める (列名は自動で読み替える):
#   (a) 実験装置の出力 (タブ区切り。ヘッダ: 日時, t, Time, V, Heat Load, 101_HB, 102_EC1, ..., 116_RES_2, ...)
#   (b) 整形済み (カンマ区切り。列: t, Q, HB, EC1, EC2, VL1..3, CON1..5, LL1..3, RES1, RES2)
#
# 使い方:
#   julia --project plot/plot_exp_startup_temps.jl                 … 下の設定を使う
#   julia --project plot/plot_exp_startup_temps.jl EXP_CSV [出力.png]
# =============================================================================
using Pkg; Pkg.activate(".")
using Plots, DataFrames, CSV, Statistics

# CSV.jl の警告を抑えるオプションは、バージョンで名前が違う (0.10 系: silencewarnings、1.0 以降: on_error)
const CSV_QUIET = pkgversion(CSV) >= v"1.0.0" ? (on_error = :collect,) : (silencewarnings = true,)

# ---------------------------- 設定 ----------------------------
EXP_CSV = joinpath(@__DIR__, "..", "data", "exp_5w_startup_heattop.csv")   # 下面加熱なら "exp_5w_startup.csv"
OUT_PNG = nothing             # nothing なら EXP_CSV と同じフォルダに <CSV名>_temps_to_onset+60s.png
T_AFTER = 60.0                # 加熱開始から何秒後まで描くか [s]
T_START = nothing             # 描き始めの時刻 [s] (実験 CSV の t)。nothing なら実験の最初の時刻
# 加熱開始: nothing なら Heat Load から読み取る (Q が Q_ON を初めて超えたサンプルと、その直前のサンプルの中点)。数値を入れるとその時刻 [s]
T_HEAT_EXP = nothing
Q_ON = nothing                # nothing なら max(0.5 W, 加熱後の Q の中央値 × 0.2)
# 横軸: :absolute = 実験 CSV の時刻 t [s]、:heating = 加熱開始からの時間 τ [s] (加熱前は負)
XAXIS = :absolute
RES_MEAN = false              # true なら RES1 と RES2 の平均を 1 本で描く (false は 2 本)
YLIMS = nothing               # 縦軸の範囲 [°C]。例: (19.0, 22.0)。nothing なら自動
SAVE_FIGURE = true
# --------------------------------------------------------------

length(ARGS) >= 1 && (EXP_CSV = ARGS[1])
length(ARGS) >= 2 && (OUT_PNG = ARGS[2])
OUT_PNG === nothing && (OUT_PNG = joinpath(dirname(EXP_CSV), splitext(basename(EXP_CSV))[1] * "_temps_to_onset+" * string(round(Int, T_AFTER)) * "s.png"))

# 実験 CSV の読み込み (plot_startup_compare.jl と同じ)。区切りはタブかカンマ (1 行目で判定)。
# 列名は先頭の "101_" などの番号と空白を除き、"Heat Load" → Q、"RES_2" → RES2 に読み替える
function load_exp(path)
    first_line = open(readline, path)
    delim = count(==('\t'), first_line) >= count(==(','), first_line) ? '\t' : ','
    d = CSV.read(path, DataFrame; delim = delim, normalizenames = false, CSV_QUIET...)
    function norm_name(s)
        s = strip(String(s))
        s == "Heat Load" && return "Q"
        m = match(r"^\d+_(.+)$", s)
        m === nothing && return s
        n = m.captures[1]
        return n == "RES_2" ? "RES2" : n
    end
    rename!(d, [norm_name(n) for n in names(d)]; makeunique = true)
    for need in ("t", "Q", "EC1", "VL1", "VL2", "VL3", "CON1", "CON2", "CON3", "CON4", "CON5", "RES1", "RES2")
        need in names(d) || error("実験 CSV に列 $(need) がありません。列: " * join(names(d), ", "))
    end
    dropmissing!(d, ["t", "Q"])
    return d
end
df = load_exp(EXP_CSV)

# 加熱開始の読み取り
function exp_onset(t, Q)
    on = Q .> 0.5
    any(on) || error("Heat Load (Q) が 0.5 W を超えません。T_HEAT_EXP を数値で入れてください")
    q_on = Q_ON === nothing ? max(0.5, 0.2 * median(Q[on])) : Q_ON
    i = findfirst(>(q_on), Q)
    t_on = i > 1 ? 0.5 * (t[i-1] + t[i]) : t[i]
    println("加熱開始: t = ", round(t_on, digits = 2), " s (Heat Load が ", round(q_on, digits = 2), " W を超えたサンプル t = ", t[i],
            " s と直前のサンプルの中点。加熱後の Q の中央値 ", round(median(Q[on]), digits = 2), " W)")
    return t_on
end
t0 = T_HEAT_EXP === nothing ? exp_onset(df.t, df.Q) : T_HEAT_EXP

# 描く範囲
t_lo = T_START === nothing ? minimum(df.t) : T_START
t_hi = t0 + T_AFTER
d = df[(df.t .>= t_lo) .& (df.t .<= t_hi), :]
println("描く範囲: t = ", round(t_lo, digits = 2), "〜", round(t_hi, digits = 2), " s (", nrow(d), " 点)")
x = XAXIS === :heating ? d.t .- t0 : d.t
xlabel = XAXIS === :heating ? "time from heating onset [s]" : "time [s]"
x_on = XAXIS === :heating ? 0.0 : t0

# 系列: (列名, 凡例, 色, 線種)。色は既存の図 (HB = Red、VL1 = Orange、CON1〜3 = MediumBlue / DodgerBlue / DeepSkyBlue、
# LL2 = Indigo、RES1 = DarkGreen) に合わせる。そこに無いセンサは同じ部位の色の系統から選ぶ:
#   EVA1 = Red (蒸発器。HB と同じ赤系。この図には HB を描かないので Red)、RES2 = ForestGreen (RES の緑系)、
#   VL2 = DarkOrange・VL3 = Goldenrod (VL のオレンジ系)、CON4 = SteelBlue・CON5 = LightSkyBlue (CON の青系)
series = Any[("EC1", "EVA1", :Red, :solid)]
if RES_MEAN
    push!(series, (["RES1", "RES2"], "RES (mean)", :DarkGreen, :solid))
else
    push!(series, ("RES1", "RES1", :DarkGreen, :solid), ("RES2", "RES2", :ForestGreen, :solid))
end
append!(series, [("VL1", "VL1", :Orange, :solid), ("VL2", "VL2", :DarkOrange, :solid), ("VL3", "VL3", :Goldenrod, :solid),
                 ("CON1", "CON1", :MediumBlue, :solid), ("CON2", "CON2", :DodgerBlue, :solid), ("CON3", "CON3", :DeepSkyBlue, :solid),
                 ("CON4", "CON4", :SteelBlue, :solid), ("CON5", "CON5", :LightSkyBlue, :solid)])

col_values(c::String) = Float64.(coalesce.(d[!, c], NaN))
col_values(c::Vector{String}) = sum(col_values(ci) for ci in c) ./ length(c)

fig = plot(; title = "Experiment: " * basename(EXP_CSV) * " (start → heating onset + " * string(round(Int, T_AFTER)) * " s)",
           xlabel = xlabel, ylabel = "T [°C]", legend = :outerright, grid = true, framestyle = :box,
           size = (1100, 550), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, titlefontsize = 10)
for (c, lab, colr, ls) in series
    plot!(fig, x, col_values(c); label = lab, color = colr, ls = ls, linewidth = 1.5)
end
vline!(fig, [x_on]; color = :gray, ls = :dot, lw = 1.2, label = "heating onset")
YLIMS === nothing || ylims!(fig, YLIMS)

# 加熱開始直前 (12 s) の平均と、終了時刻の値 [°C] を表示 (図の読み取りの補助)
pre = (d.t .>= t0 - 12.0) .& (d.t .< t0)
println("\n", rpad("センサ", 8), lpad("加熱前12s平均", 14), lpad("終了時刻", 10), lpad("変化", 8))
for (c, lab, _, _) in series
    v = col_values(c)
    v_pre = any(pre) ? mean(filter(!isnan, v[pre])) : NaN
    v_end = v[findlast(!isnan, v)]
    println(rpad(lab, 8), lpad(string(round(v_pre, digits = 3)), 14), lpad(string(round(v_end, digits = 3)), 10),
            lpad(string(round(v_end - v_pre, digits = 3)), 8))
end

if SAVE_FIGURE
    savefig(fig, OUT_PNG)
    println("\n保存しました: ", OUT_PNG)
end
display(fig)