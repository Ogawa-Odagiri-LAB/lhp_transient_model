# =============================================================================
# plot_Nsweep.jl  (CC の液の層数 N = 1, 2, 4, 8, 16 の :short の結果を重ねて描く。2026-11-05。未実行)
#   入力: res/transient_flat_ammonia_ic_C2_N<n>/summary.csv (ic_Nsweep.jl の出力)、実験 data/exp_5W_startup.csv
#   出力: res/ic_C2_Nsweep_plots/ に
#         Nsweep_timeseries.png … ΔP_cc (実験の ΔP_RES と重ねる)、T_l1 − T_lN、CC の頂部の壁温度の変化 ΔT_ccB (実験の RES と重ねる)、Q_nc_l、ṁ_cond (τ = −5〜90 s)
#         Nsweep_metrics.png / Nsweep_metrics.csv … N に対する、圧力の山の高さと時刻、max(T_l1−T_lN)、τ=100 s の ΔT_ccB (頂部の壁温度の変化)、層 1 の厚さ
#   使い方: julia --project plot/plot_Nsweep.jl          (プロジェクトのルートで)
#   時間: τ = t − 10 s (モデル)、τ = t − 303.8 s (実験)。ΔP は、モデルは t < 10 s の平均、実験は τ = −30〜0 s の平均からの変化。
#   圧力にフィルタはかけない (生のデータ)。実験の圧力計のノイズは rms 約 2 kPa。
# =============================================================================

using Plots, CSV, DataFrames, Statistics, Printf

# ---------------------------- 設定 ----------------------------
NS        = [1, 2, 4, 8, 16]
TAG_BASE  = "_ic_C2_N"
RES_ROOT  = joinpath(@__DIR__, "..", "res")
EXP_CSV   = joinpath(@__DIR__, "..", "data", "exp_5W_startup.csv")
OUT_DIR   = joinpath(RES_ROOT, "ic_C2_Nsweep_plots")      # 既にあれば、ファイルを上書きする (結果の res/ のフォルダは変えない)
T_HEAT_MODEL, T_HEAT_EXP = 10.0, 303.8
XLIM      = (-5.0, 90.0)
R_RATIO, H_LIQ_MM = 1.23, 9.0                              # 層 1 の厚さの見積もり用 (CC_LAYER_RATIO と h_liq ≈ 9 mm)
# ----------------------------------------------------------------

mkpath(OUT_DIR)
K = 273.15
result_dir(n) = joinpath(RES_ROOT, "transient_flat_ammonia" * TAG_BASE * string(n))
have = [n for n in NS if isfile(joinpath(result_dir(n), "summary.csv"))]
isempty(have) && error("結果が見つかりません: ", result_dir(NS[1]))
for n in setdiff(NS, have)
    println("注意: N = ", n, " の summary.csv がありません (", result_dir(n), ")。描きません")
end

# ---------- 実験 ----------
function load_exp(path)
    line = open(readline, path)
    delim = count(==('\t'), line) >= count(==(','), line) ? '\t' : ','
    d = CSV.read(path, DataFrame; delim = delim, normalizenames = false, silencewarnings = true)
    nm = strip.(String.(names(d)))
    rename!(d, nm)
    pick(rx) = (i = findfirst(n -> occursin(rx, n), nm); i === nothing ? error("実験 CSV に ", rx, " の列がありません。列: ", join(nm, ", ")) : nm[i])
    tcol = findfirst(==("t"), nm); tcol === nothing && error("実験 CSV に t の列がありません")
    out = DataFrame(t = Float64.(coalesce.(d[!, "t"], NaN)),
                    P = Float64.(coalesce.(d[!, pick(r"P_RES")], NaN)) .* 1000.0,     # MPa → kPa
                    RES1 = Float64.(coalesce.(d[!, pick(r"RES1")], NaN)), RES2 = Float64.(coalesce.(d[!, pick(r"RES_?2")], NaN)))
    return filter(r -> isfinite(r.t) && isfinite(r.P), out)
end
exp_df = load_exp(EXP_CSV)
tau_e = exp_df.t .- T_HEAT_EXP
pre_e = (tau_e .>= -30) .& (tau_e .< 0)
dP_e  = exp_df.P .- mean(exp_df.P[pre_e])
res_e = (exp_df.RES1 .+ exp_df.RES2) ./ 2
dRES_e = res_e .- mean(res_e[pre_e])
sel_e = (tau_e .>= XLIM[1]) .& (tau_e .<= XLIM[2])

# ---------- モデル ----------
struct Run; n::Int; tau::Vector{Float64}; dP::Vector{Float64}; dT12::Vector{Float64}; dTcc::Vector{Float64}; Qnc::Vector{Float64}; mc::Vector{Float64}; end   # dTcc = CC の頂部の壁温度 T_ccB の変化 (T_sat(P) ではない)
function load_run(n)
    df = CSV.read(joinpath(result_dir(n), "summary.csv"), DataFrame)
    tau = df.t .- T_HEAT_MODEL
    pre = tau .< 0
    any(pre) || (pre = tau .< 0.5)
    # CC ケースの頂部側の壁節点 T_ccB (= profiles.csv の comp = "cc" の T_wall = plot_startup_compare.jl の RES モデル)。T_ccA は蒸発器側の節点
    hasproperty(df, :T_ccB) || error("summary.csv に T_ccB の列がありません: ", result_dir(n))
    return Run(n, tau, (df.P_cc .- mean(df.P_cc[pre])) ./ 1e3, df.T_l1 .- df.T_lN, df.T_ccB .- mean(df.T_ccB[pre]),
               1e3 .* df.Q_nc_l, 1e6 .* df.mdot_cond)             # ΔP [kPa]、T_l1−T_lN [K]、ΔT_cc [K]、Q_nc_l [mW]、ṁ_cond [mg/s]
end
runs = [load_run(n) for n in have]
col(i) = cgrad(:viridis)[length(runs) == 1 ? 0.5 : (i - 1) / (length(runs) - 1)]
lab(n) = n == 1 ? "N = 1 (集中定数)" : "N = $(n)"

# ---------- 時系列 ----------
function base(ttl, yl; last = false)
    plot(; title = ttl, ylabel = yl, xlabel = last ? "time from heating onset, τ [s]" : "", xlims = XLIM, legend = :outerright,
         grid = true, framestyle = :box)
end
function addruns!(p, f; kw...)
    for (i, r) in enumerate(runs)
        m = (r.tau .>= XLIM[1]) .& (r.tau .<= XLIM[2])
        plot!(p, r.tau[m], f(r)[m]; color = col(i), lw = 2, label = lab(r.n), kw...)
    end
    vline!(p, [0.0]; color = :gray, ls = :dot, label = "")
    return p
end
p1 = base("CC pressure change ΔP (model P_cc, experiment P_RES)", "ΔP [kPa]")
addruns!(p1, r -> r.dP)
plot!(p1, tau_e[sel_e], dP_e[sel_e]; color = :black, lw = 1.2, alpha = 0.7, label = "P_RES exp (raw)")
p2 = base("Liquid layer 1 − layer N (model)", "T_l1 − T_lN [K]")
addruns!(p2, r -> r.dT12)
p3 = base("CC wall temperature change (model: top-side wall node T_ccB, experiment: RES mean)", "ΔT [K]")
addruns!(p3, r -> r.dTcc)
plot!(p3, tau_e[sel_e], dRES_e[sel_e]; color = :black, lw = 1.2, alpha = 0.7, label = "RES exp (mean of RES1, RES2)")
p4 = base("Natural-convection heat into the liquid (secondary wick → liquid)", "Q_nc_l [mW]")
addruns!(p4, r -> r.Qnc)
p5 = base("Interface condensation rate", "ṁ_cond [mg/s]"; last = true)
addruns!(p5, r -> r.mc)
fig = plot(p1, p2, p3, p4, p5; layout = (5, 1), size = (1000, 1400), left_margin = 6Plots.mm, bottom_margin = 2Plots.mm)
savefig(fig, joinpath(OUT_DIR, "Nsweep_timeseries.png"))
println("保存しました: ", joinpath(OUT_DIR, "Nsweep_timeseries.png"))

# ---------- N に対する指標 ----------
rows = NamedTuple[]
for r in runs
    m = (r.tau .>= 15) .& (r.tau .<= 100)
    ip = argmax(ifelse.(m, r.dP, -Inf)); it = argmax(ifelse.(m, r.dT12, -Inf)); i100 = argmin(abs.(r.tau .- 100))
    g1 = r.n == 1 ? 1.0 : (R_RATIO - 1) / (R_RATIO^r.n - 1)
    push!(rows, (N = r.n, layer1_mm = g1 * H_LIQ_MM, dP_peak_kPa = r.dP[ip], tau_dP_peak_s = r.tau[ip],
                 dT1N_max_K = r.dT12[it], tau_dT1N_max_s = r.tau[it], dTccB_100s_K = r.dTcc[i100], Qnc_max_mW = maximum(r.Qnc[m]),
                 mcond_max_mgps = maximum(r.mc[m])))
end
tab = DataFrame(rows)
# 実験の値 (比較用)
m_e = (tau_e .>= 15) .& (tau_e .<= 100)
ipe = argmax(ifelse.(m_e, dP_e, -Inf))
exp_peak, exp_tau = dP_e[ipe], tau_e[ipe]
CSV.write(joinpath(OUT_DIR, "Nsweep_metrics.csv"), tab)
show(stdout, "text/plain", tab); println()
@printf("実験: ΔP_RES の山 (τ = 15〜100 s の最大、生データ) %.1f kPa (τ = %.1f s)\n", exp_peak, exp_tau)

xt = (tab.N, string.(tab.N))
q1 = plot(tab.N, tab.dP_peak_kPa; marker = :circle, lw = 2, xscale = :log2, xticks = xt, label = "model", color = :blue,
          title = "Pressure peak (τ = 15–100 s)", ylabel = "ΔP peak [kPa]", framestyle = :box, grid = true, legend = :best)
hline!(q1, [exp_peak]; color = :black, ls = :dash, label = "P_RES exp (raw max)")
q2 = plot(tab.N, tab.dT1N_max_K; marker = :circle, lw = 2, xscale = :log2, xticks = xt, label = "", color = :red,
          title = "max(T_l1 − T_lN)", ylabel = "[K]", framestyle = :box, grid = true)
q3 = plot(tab.N, tab.dTccB_100s_K; marker = :circle, lw = 2, xscale = :log2, xticks = xt, label = "", color = :green,
          title = "ΔT_ccB (CC top wall) at τ = 100 s", ylabel = "[K]", xlabel = "number of liquid layers N", framestyle = :box, grid = true)
q4 = plot(tab.N, tab.layer1_mm; marker = :circle, lw = 2, xscale = :log2, yscale = :log10, xticks = xt, label = "", color = :gray,
          title = "Thickness of layer 1 (r = $(R_RATIO), h_liq ≈ $(H_LIQ_MM) mm)", ylabel = "[mm]", xlabel = "number of liquid layers N", framestyle = :box, grid = true)
figm = plot(q1, q2, q3, q4; layout = (2, 2), size = (1100, 800), left_margin = 6Plots.mm, bottom_margin = 4Plots.mm)
savefig(figm, joinpath(OUT_DIR, "Nsweep_metrics.png"))
println("保存しました: ", joinpath(OUT_DIR, "Nsweep_metrics.png"), "、", joinpath(OUT_DIR, "Nsweep_metrics.csv"))
display(fig)