# =============================================================================
# plot_startup_compare.jl  (実験とモデルの起動の比較プロット)
#
# 入力:
#   ・モデルの結果フォルダの summary.csv と profiles.csv (save_results の出力)
#   ・実験の CSV。次のどちらでも読める (列名は自動で読み替える):
#       (a) 実験装置の出力 (タブ区切り。ヘッダ: 日時, t, Time, V, Heat Load, 101_HB, 102_EC1, 103_EC2, 104_VL1, ...,
#           115_RES1, 116_RES_2, 117_cond, 118_coldplate, ..., 121_P_EVA, 122_P_RES。t は [s]、温度は [°C]、圧力は [MPa])
#       (b) 整形済み (カンマ区切り。列: t, Q, HB, EC1, EC2, VL1..3, CON1..5, LL1..3, RES1, RES2。温度は °C)
# 出力:
#   時系列の図 (HB/EC/RES、VL1-3、CON1-5、LL1-3、圧力と印加熱) を 1 枚の PNG に保存する。
#
# 使い方:
#   julia --project plot/plot_startup_compare.jl                          … 下の設定を使う
#   julia --project plot/plot_startup_compare.jl RES_DIR EXP_CSV [出力.png]
#
# モデルのセンサ位置 (transient_flat.jl と同じ。管の全長に対する比 s):
#   VL1/2/3 = 1/6, 1/2, 5/6 (蒸発器側から)、CON1..5 = 0, 0.25, 0.5, 0.75, 1 (入口から)、
#   LL1/2/3 = 0.10/1.15, 0.65/1.15, 1.11/1.15 (凝縮器側から)。壁温度 (profiles.csv の T_wall) を、
#   セル中心 (i−0.5)/N の間で線形補間して使う。HB・EC は集中容量の壁 (hb, evap)、RES は CC の壁 (cc)。
# =============================================================================

using Plots, DataFrames, CSV, Statistics, CoolProp
using DSP        # 圧力のローパスフィルタ (P_FILTER = true のとき)

# ---------------------------- 設定 ----------------------------
RES_DIR = joinpath(@__DIR__, "..", "res", "transient_flat_ammonia_2h_v9a_hcon_geo2_foot")
EXP_CSV = joinpath(@__DIR__, "..", "data", "exp_5W_startup.csv")
OUT_PNG = nothing             # nothing なら RES_DIR/startup_compare.png
# 時間軸:
#   :absolute … 実験・モデルとも、それぞれの時刻 t をそのまま使う (実験の t=0 とモデルの t=0 が同じ計算のとき)
#   :heating  … 各自の加熱開始 (Q>0 になる最初の時刻) からの時間 τ にそろえる (:short の計算など、加熱開始が違うとき)
ALIGN = :heating                 # :auto … モデルの計算が 300 s 未満 (:short など) なら :heating、それ以外は :absolute
# 加熱開始の時刻 [s]。計算ログの τ = 0 と同じ基準 (加熱を指令した時刻) にそろえる。nothing なら Q>0.5 W になる最初の時刻から自動で求める
# (自動だと、モデルのヒーター遅れ 3.5 s のため、モデルは指令より約 4 s、実験は約 1 s 遅い時刻になり、図で 3 s 前後ずれる)。
T_HEAT_MODEL = 10.0           # transient_flat.jl の t_heat
T_HEAT_EXP   = 303.8
XRANGE_MIN = nothing          # 横軸の範囲 [分]。nothing なら自動 (:absolute は 0〜12 分、:heating は加熱前 0.2 分〜モデルの終了時刻)
SHOW_EXP_EC_EACH = false      # true なら EC1, EC2 / RES1, RES2 を別々に描く (false は平均)
SAVE_FIGURE = true
# 縦軸の範囲 [°C や kPa や K]。パネルの番号: 1 = HB/EC/RES、2 = VL、3 = CON、4 = LL、5 = 圧力、6 = ΔLL3 と ΔT_sat。
# 例: YLIMS = Dict(1 => (22.0, 25.0), 3 => (20.0, 23.0))。書かないパネルは、表示範囲 (横軸) の中のデータから自動 (外れ値は両端 0.2 % を除く)。
YLIMS = Dict{Int,Tuple{Float64,Float64}}()
# 圧力: 実験は P_RES (CC の圧力計。単位 MPa と想定)。圧力計の絶対値はモデルとずれる (加熱前 0.864 MPa に対しモデル 0.931 MPa) ので、
# P_REL = true なら、加熱開始前 P_BASE_S 秒の平均からの変化 ΔP を描く。false なら絶対値 [kPa]。
P_REL = true
P_BASE_S = 30.0
P_EXP_COL = "P_RES"           # "P_EVA" にすると蒸発器側の圧力計
P_EXP_TO_KPA = 1000.0         # 実験の圧力の単位が MPa のとき 1000
# 圧力のローパスフィルタ (圧力のパネル 5 と 6 だけ。温度には使わない)。実験とモデルの両方に同じフィルタをかける (公平な比較のため)。
#   DSP.jl の Butterworth (P_FILTER_ORDER 次) を前向きと後ろ向きに 1 回ずつ (filtfilt = ゼロ位相。立ち上がりの時刻がずれない)。両端は奇対称に折り返す。
#   間隔が不揃いな実験・モデルのデータは、P_FILTER_DT [s] の一様な格子に線形補間してからフィルタをかけ、元の時刻に戻す。
#   0.04 Hz: 実験の圧力のノイズ (白色、rms 約 2 kPa) を約 1/2 にする。山の高さは 0.08 Hz より 3〜8 % 低く出る (実験・モデルとも同じ)。
P_FILTER = true               # false なら生のデータを描く
P_FILTER_HZ = 0.04            # カットオフ周波数 [Hz] (P_FILTER_DT の逆数の半分 = ナイキスト未満)
P_FILTER_ORDER = 4
P_FILTER_DT = 2.0             # [s]
# センサ位置 (管の全長に対する比)
S_VL  = [1/6, 1/2, 5/6]
S_CON = [0.0, 0.25, 0.5, 0.75, 1.0]
S_LL  = [0.10 / 1.15, 0.65 / 1.15, 1.11 / 1.15]
# --------------------------------------------------------------

length(ARGS) >= 1 && (RES_DIR = ARGS[1])
length(ARGS) >= 2 && (EXP_CSV = ARGS[2])
length(ARGS) >= 3 && (OUT_PNG = ARGS[3])
OUT_PNG === nothing && (OUT_PNG = joinpath(RES_DIR, "startup_compare.png"))

K = 273.15

# ---------- 読み込み ----------
df_s = CSV.read(joinpath(RES_DIR, "summary.csv"), DataFrame)
df_p = CSV.read(joinpath(RES_DIR, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :T_wall])

# 実験 CSV の読み込み。区切りはタブかカンマ (1 行目で判定)。列名は、先頭の "101_" などの番号と空白を除き、
# "Heat Load" → Q、"RES_2" → RES2 に読み替える。データ行の先頭が空 (日時の列) でも読める。
function load_exp(path)
    first_line = open(readline, path)
    delim = count(==('\t'), first_line) >= count(==(','), first_line) ? '\t' : ','
    d = CSV.read(path, DataFrame; delim = delim, normalizenames = false, silencewarnings = true)
    function norm_name(s)
        s = strip(String(s))
        s == "Heat Load" && return "Q"
        m = match(r"^\d+_(.+)$", s)
        m === nothing && return s
        n = m.captures[1]
        return n == "RES_2" ? "RES2" : n
    end
    rename!(d, [norm_name(n) for n in names(d)]; makeunique = true)
    for need in ("t", "Q", "HB", "EC1", "EC2", "VL1", "VL2", "VL3", "CON1", "CON2", "CON3", "CON4", "CON5",
                 "LL1", "LL2", "LL3", "RES1", "RES2")
        need in names(d) || error("実験 CSV に列 $(need) がありません。列: " * join(names(d), ", "))
    end
    dropmissing!(d, ["t", "Q", "HB"])
    return d
end
df_e = load_exp(EXP_CSV)

# ---------- 加熱開始時刻 ----------
first_heat(t, Q) = (i = findfirst(>(0.5), Q); i === nothing ? t[1] : t[i])
t0_model = T_HEAT_MODEL === nothing ? first_heat(df_s.t, df_s.Q_apply) : T_HEAT_MODEL
t0_exp   = T_HEAT_EXP === nothing ? first_heat(df_e.t, df_e.Q) : T_HEAT_EXP
println("加熱開始: モデル t = ", round(t0_model, digits = 2), " s,  実験 t = ", round(t0_exp, digits = 2), " s")
ALIGN_USED = ALIGN === :auto ? (maximum(df_s.t) < 300.0 ? :heating : :absolute) : ALIGN
println("時間軸: ", ALIGN_USED)
shift_m = ALIGN_USED === :heating ? t0_model : 0.0
shift_e = ALIGN_USED === :heating ? t0_exp : 0.0
tm(t, shift) = (t .- shift) ./ 60.0          # [分]
xlabel = ALIGN_USED === :heating ? "time from heating onset [min]" : "time [min]"

# ---------- モデルの壁温度 (profiles.csv) を、(時刻 × セル) の行列にする ----------
function wall_matrix(df, comp)
    d = sort(filter(r -> r.comp == comp, df), [:t, :cell])
    nrow(d) == 0 && return (Float64[], zeros(0, 0))
    ts = unique(d.t)
    N = maximum(d.cell)
    M = permutedims(reshape(d.T_wall, N, length(ts))) .- K      # (時刻, セル) [°C]
    return ts, M
end

# 管の全長に対する比 s の位置の温度 (セル中心 (i−0.5)/N の間の線形補間。端は外挿せず端のセル)
function at_position(M, s)
    N = size(M, 2)
    x = (collect(1:N) .- 0.5) ./ N
    s <= x[1] && return M[:, 1]
    s >= x[N] && return M[:, N]
    i = searchsortedlast(x, s)
    w = (s - x[i]) / (x[i+1] - x[i])
    return (1 - w) .* M[:, i] .+ w .* M[:, i+1]
end

t_vl, M_vl = wall_matrix(df_p, "vl")
t_cd, M_cd = wall_matrix(df_p, "cond")
t_ll, M_ll = wall_matrix(df_p, "ll")
t_hb, M_hb = wall_matrix(df_p, "hb")
t_ev, M_ev = wall_matrix(df_p, "evap")
t_cc, M_cc = wall_matrix(df_p, "cc")

# ---------- 描画 ----------
xlim = XRANGE_MIN !== nothing ? XRANGE_MIN :
       (ALIGN_USED === :heating ? (-0.2, (maximum(df_s.t) - t0_model) / 60.0) : (0.0, 12.0))
cols = palette(:tab10)
exp_style = (lw = 1.6, alpha = 0.85)
mod_style = (lw = 2.0, ls = :dash)

# plot! と同じだが、横軸の表示範囲内のデータを記録して、縦軸の範囲の自動設定に使う
VISIBLE = Dict{UInt,Vector{Float64}}()
function addline!(p, x, y; kw...)
    plot!(p, x, y; kw...)
    m = (x .>= xlim[1]) .& (x .<= xlim[2]) .& isfinite.(y)
    append!(get!(VISIBLE, objectid(p), Float64[]), y[m])
    return p
end

function base_plot(ttl, ylab)
    plot(; title = ttl, xlabel = xlabel, ylabel = ylab, xlims = xlim, legend = :outerright,
         grid = true, framestyle = :box, size = (900, 300))
end

# 1) HB / EC / RES
p1 = base_plot("HB / EC / RES (wall temperature)", "T [°C]")
addline!(p1, tm(df_e.t, shift_e), df_e.HB; color = cols[1], label = "HB exp", exp_style...)
addline!(p1, tm(t_hb, shift_m), M_hb[:, 1]; color = cols[1], label = "HB model", mod_style...)
if SHOW_EXP_EC_EACH
    addline!(p1, tm(df_e.t, shift_e), df_e.EC1; color = cols[2], label = "EC1 exp", exp_style...)
    addline!(p1, tm(df_e.t, shift_e), df_e.EC2; color = cols[2], ls = :dot, label = "EC2 exp", lw = 1.6)
    addline!(p1, tm(df_e.t, shift_e), df_e.RES1; color = cols[3], label = "RES1 exp", exp_style...)
    addline!(p1, tm(df_e.t, shift_e), df_e.RES2; color = cols[3], ls = :dot, label = "RES2 exp", lw = 1.6)
else
    addline!(p1, tm(df_e.t, shift_e), (df_e.EC1 .+ df_e.EC2) ./ 2; color = cols[2], label = "EC exp (mean)", exp_style...)
    addline!(p1, tm(df_e.t, shift_e), (df_e.RES1 .+ df_e.RES2) ./ 2; color = cols[3], label = "RES exp (mean)", exp_style...)
end
addline!(p1, tm(t_ev, shift_m), M_ev[:, 1]; color = cols[2], label = "EC model", mod_style...)
addline!(p1, tm(t_cc, shift_m), M_cc[:, 1]; color = cols[3], label = "RES model", mod_style...)

# 2) VL1-3
p2 = base_plot("Vapor line (VL1-3, wall temperature)", "T [°C]")
for (k, name) in enumerate(("VL1", "VL2", "VL3"))
    addline!(p2, tm(df_e.t, shift_e), df_e[!, name]; color = cols[k], label = name * " exp", exp_style...)
    addline!(p2, tm(t_vl, shift_m), at_position(M_vl, S_VL[k]); color = cols[k], label = name * " model", mod_style...)
end

# 3) CON1-5
p3 = base_plot("Condenser (CON1-5, wall temperature)", "T [°C]")
for (k, name) in enumerate(("CON1", "CON2", "CON3", "CON4", "CON5"))
    addline!(p3, tm(df_e.t, shift_e), df_e[!, name]; color = cols[k], label = name * " exp", exp_style...)
    addline!(p3, tm(t_cd, shift_m), at_position(M_cd, S_CON[k]); color = cols[k], label = name * " model", mod_style...)
end

# 4) LL1-3
p4 = base_plot("Liquid line (LL1-3, wall temperature)", "T [°C]")
for (k, name) in enumerate(("LL1", "LL2", "LL3"))
    addline!(p4, tm(df_e.t, shift_e), df_e[!, name]; color = cols[k], label = name * " exp", exp_style...)
    addline!(p4, tm(t_ll, shift_m), at_position(M_ll, S_LL[k]); color = cols[k], label = name * " model", mod_style...)
end

# 5) 圧力 (モデルの P_cc と、実験の圧力計)
function lp_filter(t, y)
    P_FILTER || return y
    ok = isfinite.(t) .& isfinite.(y)
    count(ok) < 8 && return y
    tt, yy = t[ok], y[ok]
    tg = collect(tt[1]:P_FILTER_DT:tt[end])
    length(tg) < 8 && return y
    interp(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
    yg = [interp(tt, yy, x) for x in tg]
    n = length(yg); pad = min(100, n - 1)
    ext = [2 * yg[1] .- yg[pad+1:-1:2]; yg; 2 * yg[end] .- yg[end-1:-1:end-pad]]
    h = digitalfilter(Lowpass(P_FILTER_HZ / (0.5 / P_FILTER_DT)), Butterworth(P_FILTER_ORDER))
    zg = filtfilt(h, ext)[pad+1:pad+n]
    z = copy(y)
    z[ok] = [interp(tg, zg, x) for x in tt]
    return z
end
P_TAG = P_FILTER ? " (LP " * string(P_FILTER_HZ) * " Hz)" : ""
exp_pressure() = lp_filter(df_e.t, df_e[!, P_EXP_COL] .* P_EXP_TO_KPA)
p5 = base_plot(P_REL ? "CC pressure change from pre-heating" : "CC pressure", P_REL ? "ΔP [kPa]" : "P [kPa]")
base(t, P, t0) = (m = (t .>= t0 - P_BASE_S) .& (t .< t0); any(m) ? mean(P[m]) : P[1])
P_m = lp_filter(df_s.t, df_s.P_cc ./ 1e3)
P_m_plot = P_REL ? P_m .- base(df_s.t, P_m, t0_model) : P_m
addline!(p5, tm(df_s.t, shift_m), P_m_plot; color = cols[4], label = "P_cc model" * P_TAG, lw = 2.0)
# 液相の最下層 (層 N、ウィックや液管につながる側 = バルクに近い液) の温度 T_lN から決まる飽和圧。P_cc との差が「界面の層が熱くなったことによる過剰圧力」(旧版の意味の ΔP_ne)。
if "T_lN" in names(df_s) && !all(isnan, df_s.T_lN)
    Psat_lN = lp_filter(df_s.t, [isnan(T) ? NaN : CoolProp.PropsSI("P", "T", T, "Q", 0, "Ammonia") / 1e3 for T in df_s.T_lN])   # [kPa]
    addline!(p5, tm(df_s.t, shift_m), P_REL ? Psat_lN .- base(df_s.t, P_m, t0_model) : Psat_lN; color = cols[2], label = "P_sat(T_lN) model (liquid bottom layer)" * P_TAG, ls = :dashdot, lw = 2.0)
end
if P_EXP_COL in names(df_e)
    Pe = exp_pressure()
    Pe_plot = P_REL ? Pe .- base(df_e.t, Pe, t0_exp) : Pe
    addline!(p5, tm(df_e.t, shift_e), Pe_plot; color = cols[4], label = P_EXP_COL * " exp" * P_TAG, exp_style...)
end
# 6) LL3・VL1 の加熱前からの変化と、圧力から求めた飽和温度の変化 ΔT_sat(ΔP)。
#    LL3 が ΔT_sat に沿うなら、LL3 の流体は飽和 (気液接触) に近い。沿わないなら、熱伝導や流れの影響。
#    ΔT_sat は、加熱前の RES 平均温度の飽和圧 P0 に ΔP を足して、T_sat(P0+ΔP) − T_sat(P0) (圧力計の絶対値のずれは使わない)。
p6 = base_plot("LL3, VL1 change vs saturation-temperature change from ΔP", "ΔT [K]")
function tsat_change(t, dP_kPa, t0, T0_C)
    P0 = CoolProp.PropsSI("P", "T", T0_C + K, "Q", 0, "Ammonia")
    T_ref = CoolProp.PropsSI("T", "P", P0, "Q", 0, "Ammonia")
    return [CoolProp.PropsSI("T", "P", P0 + 1e3 * x, "Q", 0, "Ammonia") - T_ref for x in dP_kPa]
end
pre_e = (df_e.t .>= t0_exp - P_BASE_S) .& (df_e.t .< t0_exp)
any(pre_e) || (pre_e = df_e.t .< t0_exp)
T_res0 = mean((df_e.RES1[pre_e] .+ df_e.RES2[pre_e]) ./ 2)
dLL3_e = df_e.LL3 .- mean(df_e.LL3[pre_e])
dVL1_e = df_e.VL1 .- mean(df_e.VL1[pre_e])
addline!(p6, tm(df_e.t, shift_e), dLL3_e; color = cols[3], label = "ΔLL3 exp", exp_style...)
addline!(p6, tm(df_e.t, shift_e), dVL1_e; color = cols[1], label = "ΔVL1 exp", exp_style...)
if P_EXP_COL in names(df_e)
    Pe = exp_pressure()
    dPe = Pe .- mean(Pe[pre_e])
    addline!(p6, tm(df_e.t, shift_e), tsat_change(df_e.t, dPe, t0_exp, T_res0); color = :black, label = "ΔT_sat from P exp", lw = 1.6)
end
pre_m = (df_s.t .>= t0_model - P_BASE_S) .& (df_s.t .< t0_model)
any(pre_m) || (pre_m = df_s.t .< t0_model)
dPm = P_m .- mean(P_m[pre_m])
addline!(p6, tm(df_s.t, shift_m), tsat_change(df_s.t, dPm, t0_model, T_res0); color = :black, ls = :dash, label = "ΔT_sat from P model", lw = 2.0)
addline!(p6, tm(t_ll, shift_m), at_position(M_ll, S_LL[3]) .- mean(at_position(M_ll, S_LL[3])[t_ll .< t0_model]); color = cols[3], label = "ΔLL3 model", mod_style...)

# 加熱開始の目印 (全図)
for p in (p1, p2, p3, p4, p5, p6)
    vline!(p, [ALIGN_USED === :heating ? 0.0 : t0_exp / 60.0]; color = :gray, ls = :dot, label = "")
end

# 縦軸の範囲
for (k, p) in enumerate((p1, p2, p3, p4, p5, p6))
    if haskey(YLIMS, k)
        ylims!(p, YLIMS[k])
    elseif haskey(VISIBLE, objectid(p)) && !isempty(VISIBLE[objectid(p)])
        v = VISIBLE[objectid(p)]
        lo, hi = quantile(v, 0.002), quantile(v, 0.998)
        pad = 0.06 * (hi - lo) + 1e-6
        ylims!(p, (lo - pad, hi + pad))
    end
end

fig = plot(p1, p2, p3, p4, p5, p6; layout = (6, 1), size = (1000, 1800), left_margin = 6Plots.mm, bottom_margin = 3Plots.mm)
if SAVE_FIGURE
    savefig(fig, OUT_PNG)
    println("保存しました: ", OUT_PNG)
end
display(fig)