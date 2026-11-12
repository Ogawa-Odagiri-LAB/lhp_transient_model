# plot_exp_eva_vl_p.jl — experiment: (EVA1 - VL1) and P_EVA versus time from heating onset, tau = 15-25 s
# Run from the project root:  julia scripts/plot_exp_eva_vl_p.jl   (edit the settings below; no arguments)
# Input : data/<ev_file>   Output: res/exp_plots/<ev_out>.png (+ .pdf/.svg)
# tau = 0 is the first time Heat Load >= ev_onset_W.

# ---------------------------- settings ----------------------------
ev_file    = "exp_5W_startup.csv"
ev_out     = "exp_EVA1-VL1_PEVA_15-25s"
ev_tcol    = "t"                    # time column [s] ("t" or "Time")
ev_heat    = "HeatLoad"             # heat load column
ev_eva     = "EVA1"                 # also tried: EC1
ev_vl      = "VL1"
ev_pcol    = "PEVA"                 # P_EVA column (raw name 121_P_EVA); assumed MPa
ev_pscale  = 1000.0                 # MPa -> kPa
ev_onset_W = 1.0
ev_xlim    = (15.0, 25.0)
ev_xtick   = 1.0
ev_dTlim   = :auto                  # e.g. (0, 10) or :auto  [K]
ev_Plim    = :auto                  # e.g. (800, 900) or :auto [kPa]
ev_P_rel   = false                  # true: plot P_EVA - P_EVA(tau=0) instead of the absolute value
ev_size    = (900, 600)
ev_dpi     = 300
ev_font_label, ev_font_tick, ev_font_legend = 14, 12, 12
ev_lw      = 3.0
ev_ms      = 5                      # marker size (samples are shown as markers)
ev_vector  = true

# ---------------------------- code ----------------------------
using Pkg
Pkg.activate(".")
using CSV, DataFrames, Plots, Statistics

path = isabspath(ev_file) ? ev_file : joinpath(@__DIR__, "..", "data", ev_file)
df = CSV.read(path, DataFrame; header = 1, silencewarnings = true)
normh(x) = lowercase(replace(replace(replace(strip(replace(string(x), r"[\"﻿]" => "")), r"^\d+_" => ""), "_" => ""), " " => ""))
nm = Dict(normh(n) => n for n in names(df))
function col(name)
    k = normh(name)
    haskey(nm, k) || error("column '" * name * "' not found. Available: " * join(names(df), ", "))
    return [v isa Number ? Float64(v) : (something(tryparse(Float64, string(v)), NaN)) for v in df[!, nm[k]]]
end
t, q, eva, vl, pe = col(ev_tcol), col(ev_heat), col(ev_eva), col(ev_vl), col(ev_pcol)
i0 = findfirst(v -> isfinite(v) && v >= ev_onset_W, q)
i0 === nothing && error("heating onset (Heat Load >= " * string(ev_onset_W) * " W) not found")
t0 = t[i0]
τ = t .- t0
println("heating onset: t = ", t0, " s (row ", i0, ")")

sel = findall(i -> ev_xlim[1] <= τ[i] <= ev_xlim[2] && isfinite(eva[i]) && isfinite(vl[i]), eachindex(τ))
isempty(sel) && error("no data in the window")
dT = eva[sel] .- vl[sel]
P  = pe[sel] .* ev_pscale
if ev_P_rel
    ip = findlast(i -> τ[i] <= 0 && isfinite(pe[i]), eachindex(τ))
    P .-= pe[ip] * ev_pscale
end
x = τ[sel]

p = plot(x, dT; label = "EVA1 - VL1", lw = ev_lw, marker = :circle, ms = ev_ms, color = :royalblue,
         xlabel = "Time from heating onset [s]", ylabel = "EVA1 - VL1 [K]", xlims = ev_xlim, xticks = ev_xlim[1]:ev_xtick:ev_xlim[2],
         size = ev_size, dpi = ev_dpi, legend = :topleft, guidefontsize = ev_font_label, tickfontsize = ev_font_tick, legendfontsize = ev_font_legend,
         right_margin = 12Plots.mm)
ev_dTlim === :auto || plot!(p, ylims = ev_dTlim)
p2 = twinx(p)
plot!(p2, x, P; label = ev_P_rel ? "P_EVA - P_EVA(0)" : "P_EVA", lw = ev_lw, marker = :square, ms = ev_ms, color = :crimson, ls = :dash,
      ylabel = ev_P_rel ? "ΔP_EVA [kPa]" : "P_EVA [kPa]", xlims = ev_xlim, legend = :topright, grid = false,
      guidefontsize = ev_font_label, tickfontsize = ev_font_tick, legendfontsize = ev_font_legend)
ev_Plim === :auto || plot!(p2, ylims = ev_Plim)

outdir = joinpath(@__DIR__, "..", "res", "exp_plots"); mkpath(outdir)
savefig(p, joinpath(outdir, ev_out * ".png"))
if ev_vector
    savefig(p, joinpath(outdir, ev_out * ".pdf")); savefig(p, joinpath(outdir, ev_out * ".svg"))
end
CSV.write(joinpath(outdir, ev_out * ".csv"), DataFrame(tau = x, EVA1_minus_VL1 = dT, P_EVA_kPa = P))
display(p)
println("Saved: ", joinpath(outdir, ev_out * ".png"))