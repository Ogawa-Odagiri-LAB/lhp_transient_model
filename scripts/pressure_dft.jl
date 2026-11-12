# pressure_dft.jl — 実験の圧力計 (P_RES, P_EVA) の生データの波形と DFT (パワースペクトル) を見る (約 2 時間分)
#
# 【目的】 ノイズの性質 (白色ノイズ / 特定周波数のピーク / A/D の量子化 / スパイク) を見て、
#   ローパスフィルタのカットオフ (と、ノッチが必要か) を決める材料にする。フィルタ自体はこのファイルにはまだ入れていない (次の段階)。
# 【実行】 プロジェクトのルートで julia scripts/pressure_dft.jl  (設定はこのファイルの冒頭を編集。引数なし)
#   DSP.jl が無ければ最初に Pkg.add する。
# 【出力】 res/pressure_dft/ に  pressure_waveform.png (波形)、pressure_psd.png (パワースペクトル)、pressure_psd.csv (スペクトルの数値)
#   と、画面に: サンプリング間隔、量子化の刻み、スパイクの数、区間ごとのノイズの大きさ (rms)、ホワイトノイズの床、スペクトルのピーク。

# ---------------------------- 設定 ----------------------------
pd_file      = "exp_5W_startup.csv"                 # data/ の下の実験 csv (2 時間分)
pd_cols      = ["P_RES", "P_EVA"]              # 圧力の列。生の列名 ("122_P_RES" など) でも "P_RES" でもよい (先頭の番号と "_" を除いて照合)
pd_scale     = 1000.0                          # 単位の換算: 実験の圧力が MPa なら 1000 で kPa になる
pd_t_heat    = 303.8                           # 加熱開始の時刻 [s] (実験の csv の t)。τ = t − pd_t_heat
# スペクトルを出す区間 (τ [s]; -Inf / Inf は端まで)。短い区間は、低い周波数の分解能が粗くなる
pd_segments  = [("before heating", -Inf, -2.0),
                ("startup 0-120 s", 0.0, 120.0),
                ("heating 120-1800 s", 120.0, 1800.0),
                ("heating 1800 s-end", 1800.0, Inf)]
pd_nperseg   = 256                             # Welch 法の 1 区間のサンプル数 (大きいほど周波数分解能が良く、平均の回数が減る)。区間が短いと自動で小さくする
pd_zoom_tau  = (-60.0, 120.0)                  # 波形の拡大図の範囲 τ [s]
pd_out_name  = "pressure_dft"                  # 出力フォルダ (res/<名前>/)
pd_spike_k   = 6.0                             # スパイクの判定: 3 点メディアンからの差が k × 1.4826 × MAD を超えたもの

# フィルタ (スペクトルの後に、DSP.jl のゼロ位相ローパスをかける)
pf_cutoffs   = [0.04, 0.06, 0.08]              # 比べるカットオフ周波数 [Hz] (ナイキスト 0.25 Hz 未満)。pf_use の値を含めること
pf_use       = 0.06                            # 保存する csv の「本命」のカットオフ [Hz]
pf_order     = 4                               # Butterworth の次数 (前向き + 後ろ向きにかけるので、実効は 2 倍の急峻さ)
pf_pad       = 100                             # 端の処理: 両端を奇対称に折り返して付ける点数 (ゼロ位相フィルタの端の暴れを抑える)
pf_zoom_tau  = (-20.0, 80.0)                   # 起動の拡大図の範囲 τ [s]
pf_levels    = (abs = 2.0, frac = 0.5)         # 立ち上がり時刻の定義: ΔP が 2 kPa を超える最初の τ / ΔP が (τ=0〜60 s の最大の 50 %) を超える最初の τ

# ---------------------------- 準備 ----------------------------
using Pkg
Pkg.activate(".")
Base.find_package("DSP") === nothing && Pkg.add("DSP")
using CSV, DataFrames, Statistics, Printf
using Plots
using DSP

pd_normh(x) = replace(replace(strip(replace(string(x), "﻿" => "")), r"^\d+_" => ""), "_" => "")
pd_tofloat(v) = [x isa Number ? Float64(x) : (y = tryparse(Float64, strip(string(x))); y === nothing ? NaN : y) for x in v]

pd_path = joinpath(@__DIR__, "..", "data", pd_file)
isfile(pd_path) || error("実験の csv が見つかりません: " * pd_path)
df = CSV.read(pd_path, DataFrame)
hdr = names(df)
function pd_col(nm)
    j = findfirst(==(nm), hdr)
    j === nothing && (j = findfirst(x -> pd_normh(x) == pd_normh(nm), hdr))
    j === nothing && error("列 \"" * nm * "\" がありません。列名: " * string(hdr))
    return pd_tofloat(df[!, j])
end
t_raw = pd_col("t")
good_t = isfinite.(t_raw)
t_raw = t_raw[good_t]
println("読み込み: ", pd_path, "  (", length(t_raw), " 行、t = ", round(t_raw[1], digits = 2), " 〜 ", round(t_raw[end], digits = 2), " s、約 ", round((t_raw[end] - t_raw[1]) / 60, digits = 1), " 分)")

# 一様なサンプリングにそろえる (測定の間隔が一定ならそのまま)
dts = diff(t_raw)
dt = median(dts)
uniform = maximum(abs.(dts .- dt)) <= 0.05 * dt
println("サンプリング間隔: 中央値 ", round(dt, digits = 4), " s (fs = ", round(1 / dt, digits = 4), " Hz、ナイキスト ", round(0.5 / dt, digits = 4), " Hz)。",
        uniform ? "間隔は一定。" : "間隔が不揃い (最小 " * string(round(minimum(dts), digits = 3)) * "、最大 " * string(round(maximum(dts), digits = 3)) * " s) → 一様な格子に線形補間する。")
lin(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))

function pd_fill(t, y)          # 欠損 (NaN) を線形補間で埋める。埋めた数も返す
    ok = isfinite.(y)
    all(ok) && return copy(y), 0
    return [ok[i] ? y[i] : lin(t[ok], y[ok], t[i]) for i in eachindex(y)], count(!, ok)
end

sig = Dict{String,Vector{Float64}}()
nfill = Dict{String,Int}()
tu = uniform ? t_raw : collect(t_raw[1]:dt:t_raw[end])
for nm in pd_cols
    y = pd_col(nm)[good_t] .* pd_scale
    y, nf = pd_fill(t_raw, y)
    nfill[nm] = nf
    sig[nm] = uniform ? y : [lin(t_raw, y, t) for t in tu]
end
τ = tu .- pd_t_heat
fs = 1 / dt

# ---------------------------- 基本統計 (量子化・スパイク) ----------------------------
med3(y) = [median(y[max(1, i - 1):min(length(y), i + 1)]) for i in eachindex(y)]
println("\n=== 基本統計 ===")
for nm in pd_cols
    y = sig[nm]
    u = sort(unique(round.(y, digits = 9)))
    steps = diff(u)
    q = isempty(steps) ? NaN : minimum(steps[steps .> 0])
    r = y .- med3(y)
    s = 1.4826 * median(abs.(r .- median(r)))
    nsp = s > 0 ? count(abs.(r) .> pd_spike_k * s) : 0
    println(rpad(nm, 8), "欠損を補間した数 ", nfill[nm], "   平均 ", round(mean(y), digits = 3), "  最小 ", round(minimum(y), digits = 3), "  最大 ", round(maximum(y), digits = 3),
            "   値の最小の刻み ≈ ", @sprintf("%.4g", q), " (圧力計/AD の分解能の目安)   スパイク (3 点メディアンから ", pd_spike_k, "σ 超) ", nsp, " 点")
end

# ---------------------------- 波形のプロット ----------------------------
isdir(joinpath(@__DIR__, "..", "res")) || mkpath(joinpath(@__DIR__, "..", "res"))
pd_out = joinpath(@__DIR__, "..", "res", pd_out_name)
mkpath(pd_out)
pw = []
for nm in pd_cols
    p1 = plot(τ ./ 60, sig[nm], label = nm, lw = 0.8, xlabel = "τ from heating start [min]", ylabel = nm * " [kPa]", title = nm * " raw (full record)", legend = false)
    push!(pw, p1)
end
for nm in pd_cols
    m = (τ .>= pd_zoom_tau[1]) .& (τ .<= pd_zoom_tau[2])
    p2 = plot(τ[m], sig[nm][m], label = nm, lw = 1.0, marker = :circle, ms = 2, xlabel = "τ [s]", ylabel = nm * " [kPa]", title = nm * " raw (zoom)", legend = false)
    push!(pw, p2)
end
savefig(plot(pw..., layout = (2, length(pd_cols)), size = (560 * length(pd_cols), 760)), joinpath(pd_out, "pressure_waveform.png"))

# ---------------------------- パワースペクトル (Welch 法) ----------------------------
function seg_data(nm, a, b)
    m = (τ .>= a) .& (τ .< b)
    return sig[nm][m], tu[m]
end
function detrend_lin(t, y)      # 直線 (傾向) を引く
    length(y) < 3 && return y .- mean(y)
    A = [ones(length(t)) (t .- mean(t))]
    c = A \ y
    return y .- A * c
end

psd_rows = DataFrame(sensor = String[], segment = String[], freq_Hz = Float64[], psd_kPa2_per_Hz = Float64[])
colors = [:black, :red, :blue, :green, :orange, :purple]
println("\n=== 区間ごとのノイズ (直線の傾向を引いた後) ===")
println(rpad("センサ", 8), rpad("区間", 22), lpad("点数", 7), lpad("rms [kPa]", 12), lpad("白色の床 [kPa²/Hz]", 20), lpad("(床から見た rms 相当 [kPa])", 30))
pp = []
for nm in pd_cols
    pl = plot(xscale = :log10, yscale = :log10, xlabel = "frequency [Hz]", ylabel = "PSD [kPa²/Hz]", title = nm * " power spectrum (Welch, Hann)", legend = :bottomleft)
    peaks_report = String[]
    for (ks, (lab, a, b)) in enumerate(pd_segments)
        y, ts = seg_data(nm, a, b)
        if length(y) < 16
            println(rpad(nm, 8), rpad(lab, 22), lpad(string(length(y)), 7), "   (点数が少なすぎるので飛ばす)")
            continue
        end
        yd = detrend_lin(ts, y)
        nper = min(pd_nperseg, max(8, length(yd) ÷ 2))
        pg = welch_pgram(yd, nper, nper ÷ 2; fs = fs, window = hanning)
        f = collect(freq(pg)); P = collect(power(pg))
        keep = f .> 0
        f = f[keep]; P = P[keep]
        append!(psd_rows, DataFrame(sensor = fill(nm, length(f)), segment = fill(lab, length(f)), freq_Hz = f, psd_kPa2_per_Hz = P))
        plot!(pl, f, P, label = lab * " (" * string(length(yd)) * " pts)", lw = 1.4, color = colors[mod1(ks, length(colors))])
        # 床: ナイキストに近い上位 25 % の周波数の中央値
        hi = f .>= 0.75 * maximum(f)
        floor_psd = median(P[hi])
        rms = std(yd)
        println(rpad(nm, 8), rpad(lab, 22), lpad(string(length(y)), 7), lpad(@sprintf("%.4g", rms), 12), lpad(@sprintf("%.3g", floor_psd), 20),
                lpad(@sprintf("%.4g", sqrt(floor_psd * 0.5 * fs)), 30))
        # ピーク: 周りの 11 点の中央値の 5 倍を超える局所最大 (上位 4 つ)
        runmed(P, i) = median(P[max(1, i - 5):min(length(P), i + 5)])
        cand = [i for i in 2:length(P)-1 if P[i] > P[i-1] && P[i] >= P[i+1] && P[i] > 5 * runmed(P, i) && f[i] > 2 / (length(yd) * dt)]
        sort!(cand, by = i -> -P[i])
        for i in cand[1:min(4, length(cand))]
            push!(peaks_report, @sprintf("  %-22s ピーク %.4f Hz (周期 %.1f s)  PSD %.3g (周りの中央値の %.0f 倍)", lab, f[i], 1 / f[i], P[i], P[i] / runmed(P, i)))
        end
    end
    push!(pp, pl)
    println(rpad(nm, 8), "スペクトルのピーク (なければ表示なし):")
    isempty(peaks_report) ? println("  (目立つピークなし)") : foreach(println, peaks_report)
end
savefig(plot(pp..., layout = (length(pd_cols), 1), size = (900, 480 * length(pd_cols))), joinpath(pd_out, "pressure_psd.png"))
CSV.write(joinpath(pd_out, "pressure_psd.csv"), psd_rows)
println("\n保存しました: ", pd_out, " (pressure_waveform.png, pressure_psd.png, pressure_psd.csv)")
println("""
見方:
  ・白色の床がほぼ平ら (低周波から高周波まで同じ高さ) なら、ノイズは白色ノイズ。ローパスで落とせる。
  ・特定の周波数に鋭いピークがあれば、チラー・ポンプ・電源などの周期的な外乱。ローパスの前にノッチで除くとよい。
  ・「値の最小の刻み」がノイズの rms と同程度なら、ノイズの正体は A/D の量子化。ローパスよりも平滑化 (移動平均など) で足りる。
  ・信号 (立ち上がり) の帯域は 'startup 0-120 s' のスペクトルが、静かな区間 ('before heating') の床から立ち上がる周波数まで。
    カットオフはその床と立ち上がりの交わる周波数より少し上に置く。""")

# =============================================================================
# フィルタ: DSP.jl の Butterworth ローパスを、前向きと後ろ向きに 1 回ずつかける (filtfilt = ゼロ位相。立ち上がりの時刻がずれない)
# =============================================================================
println("\n=== ローパスフィルタ (Butterworth ", pf_order, " 次、filtfilt、カットオフ ", pf_cutoffs, " Hz) ===")
all(c -> 0 < c < fs / 2, pf_cutoffs) || error("カットオフはナイキスト (", fs / 2, " Hz) 未満にしてください")
pf_use in pf_cutoffs || error("pf_use は pf_cutoffs の中の値にしてください")
function pf_lowpass(y, fc)
    n = length(y)
    pad = min(pf_pad, n - 1)
    ext = [2 * y[1] .- y[pad+1:-1:2]; y; 2 * y[end] .- y[end-1:-1:end-pad]]      # 奇対称に折り返して端を延ばす
    h = digitalfilter(Lowpass(fc / (fs / 2)), Butterworth(pf_order))      # カットオフはナイキスト周波数に対する比 (DSP のバージョンによらず動く書き方)
    z = filtfilt(h, ext)
    return z[pad+1:pad+n]
end
function pf_cross(τs, y, level, from)     # τ ≥ from で y が level を (下から) 最初に超える τ (線形補間)。なければ NaN
    for i in 1:length(y)-1
        τs[i+1] < from && continue
        if y[i] < level <= y[i+1]
            return τs[i] + (level - y[i]) / (y[i+1] - y[i]) * (τs[i+1] - τs[i])
        end
    end
    return NaN
end
filt = Dict{Tuple{String,Float64},Vector{Float64}}()
for nm in pd_cols, fc in pf_cutoffs
    filt[(nm, fc)] = pf_lowpass(sig[nm], fc)
end
pre = (τ .>= -280.0) .& (τ .< -2.0)
st  = (τ .>= 0.0) .& (τ .<= 60.0)
println("\n--- 起動の立ち上がり時刻 (ΔP = フィルタ後 − 加熱前の平均) と、ノイズの減り方 ---")
println(rpad("センサ", 8), rpad("カットオフ", 12), lpad("ΔP>" * string(pf_levels.abs) * " kPa の τ [s]", 22), lpad(string(Int(100 * pf_levels.frac)) * "% の τ [s]", 14), lpad("ΔP の最大 [kPa]", 17),
        lpad("加熱前の rms (残り) [kPa]", 27), lpad("(生) [kPa]", 12), lpad("差 (生−フィルタ) の rms [kPa]", 30))
for nm in pd_cols
    base = mean(sig[nm][pre])
    for fc in pf_cutoffs
        z = filt[(nm, fc)]
        d = z .- base
        pk = maximum(d[st])
        t_abs = pf_cross(τ[st], d[st], pf_levels.abs, 0.0)
        t_frc = pf_cross(τ[st], d[st], pf_levels.frac * pk, 0.0)
        println(rpad(nm, 8), rpad(string(fc) * " Hz", 12), lpad(@sprintf("%.2f", t_abs), 22), lpad(@sprintf("%.2f", t_frc), 14), lpad(@sprintf("%.2f", pk), 17),
                lpad(@sprintf("%.3f", std(detrend_lin(tu[pre], z[pre]))), 27), lpad(@sprintf("%.3f", std(detrend_lin(tu[pre], sig[nm][pre]))), 12),
                lpad(@sprintf("%.3f", std(sig[nm] .- z)), 30))
    end
end
println("(立ち上がり時刻がカットオフを変えてもほとんど動かなければ、フィルタが時刻をずらしていない。0.5 s 以上動くなら、カットオフを上げる)")

# 図: 起動の拡大 (生 + フィルタ後)、全体
zm = (τ .>= pf_zoom_tau[1]) .& (τ .<= pf_zoom_tau[2])
pz = []
for nm in pd_cols
    q = scatter(τ[zm], sig[nm][zm], label = "raw", ms = 2.5, color = :gray, xlabel = "τ [s]", ylabel = nm * " [kPa]", title = nm * " startup: raw vs low-pass", legend = :topleft)
    for (k, fc) in enumerate(pf_cutoffs)
        plot!(q, τ[zm], filt[(nm, fc)][zm], label = string(fc) * " Hz", lw = 2, color = colors[k + 1])
    end
    push!(pz, q)
end
savefig(plot(pz..., layout = (length(pd_cols), 1), size = (900, 420 * length(pd_cols))), joinpath(pd_out, "pressure_filter_startup.png"))
pf_full = []
for nm in pd_cols
    q = plot(τ ./ 60, sig[nm], label = "raw", lw = 0.5, color = :lightgray, xlabel = "τ [min]", ylabel = nm * " [kPa]", title = nm * " full record: raw vs " * string(pf_use) * " Hz low-pass", legend = :bottomright)
    plot!(q, τ ./ 60, filt[(nm, pf_use)], label = string(pf_use) * " Hz", lw = 1.5, color = :red)
    push!(pf_full, q)
end
savefig(plot(pf_full..., layout = (length(pd_cols), 1), size = (900, 420 * length(pd_cols))), joinpath(pd_out, "pressure_filter_full.png"))

# csv: 一様な格子 (dt の間隔)。t, τ [s]、各センサの生と、各カットオフのフィルタ後 [kPa]
out = DataFrame(t = tu, tau = τ)
for nm in pd_cols
    out[!, Symbol(nm * "_raw")] = sig[nm]
    for fc in pf_cutoffs
        out[!, Symbol(nm * "_lp" * string(fc))] = filt[(nm, fc)]
    end
end
CSV.write(joinpath(pd_out, "pressure_filtered.csv"), out)
println("\n保存しました: ", pd_out, " (pressure_filter_startup.png, pressure_filter_full.png, pressure_filtered.csv)")