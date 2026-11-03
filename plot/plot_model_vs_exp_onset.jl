# plot_model_vs_exp_onset.jl — experiment vs model from the heating onset to 8 min: temperatures (3 groups, one SVG each) and CC pressure (P_RES)
#
# Run from the project root:  julia scripts/plot_model_vs_exp_onset.jl   (edit the settings below; no arguments)
# Input : data/<mp_exp_file>               experiment (same format as plot_exp_onset.jl; the onset = first time the heat load exceeds mp_onset_W)
#         res/<result>/summary.csv, profiles.csv   model results of the :long run (816 s) — one entry per model in mp_models
#         mp_models: the 16-layer CC (default) and the lumped CC (CC_LUMPED = true in types.jl, result folder ends with "_lump").
#                    A model whose folder does not exist is skipped with a warning, so the script also works with the first model only.
# Output: res/<mp_outdir>/  temp_HB_EC_RES.svg, temp_VL1-3.svg, temp_CON1-3.svg, pressure_RES.svg (all models overlaid),
#         pressure_RES_<tag>.svg (one model each), and the same files as .png if mp_png = true.
# Line styles: experiment = solid (thin, with the sensor colour), model (16 layers) = dashed, model (lumped) = dash-dot. Same colour = same sensor.
# Model positions: HB = mean of hb, EC = evaporator-case side node (eca), RES = mean of the CC wall, VL1-3 = vapor line at 1/6, 1/2, 5/6,
#         CON1-3 = condenser at 0, 0.25, 0.5 of its length. Same definitions as param_sweep.jl.

# ---------------------------- settings ----------------------------
mp_exp_file  = "exp_5W_startup.csv"              # under data/
mp_models    = [("Model", "transient_flat_ammonia_v9a"),            # (label, folder under res/)
                ("Model (lumped CC)", "transient_flat_ammonia_v9a_lump")]
mp_t_heat    = 10.0                            # heating start in the model [s] (= t_heat in transient_flat.jl)
mp_onset_W   = 1.0                             # experiment: onset = first time the heat load exceeds this [W]
mp_t_onset   = :auto                           # or a number [s] to set the experimental onset by hand (e.g. 303.8)
mp_xlim      = (0.0, 8.0)                      # [min since the heating onset]
mp_xtick     = 1.0
mp_Tlim      = :auto                           # temperature axis [°C], e.g. (18, 28)
mp_dt        = 1.0                             # sampling of the model curves [s]
mp_outdir    = "model_vs_exp_onset"
mp_size      = (1100, 650)
mp_dpi       = 300
mp_font_label, mp_font_tick, mp_font_legend = 14, 12, 11
mp_lw        = 2.5
mp_png       = true                            # also save PNG
# --- pressure ---
mp_P_exp_col = "P_RES"                         # experiment pressure column (MPa)
mp_P_scale   = 1000.0                          # MPa -> kPa
mp_P_filter  = false                            # low-pass (zero-phase Butterworth, DSP.jl) for the experiment AND the model curve
mp_P_hz      = 0.04                            # cutoff [Hz]
mp_P_order   = 4
mp_P_dt      = 2.0                             # resampling step for the filter [s]
mp_P_rel     = true                            # true: plot the change P − P(τ = 0) (default); false: absolute pressure
mp_Plim      = :auto                           # [kPa] (or ΔP if mp_P_rel), e.g. (800, 1200)
# ----------------------------------------------------------------

using Pkg
Pkg.activate(".")
using CSV, DataFrames, Plots, Statistics, DSP

normh(x) = lowercase(replace(replace(strip(replace(string(x), "﻿" => "")), r"^\d+_" => ""), r"[ _\"]" => ""))
tofloat(v) = [x isa Number ? Float64(x) : (y = tryparse(Float64, strip(string(x))); y === nothing ? NaN : y) for x in v]
lin(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
const K0 = 273.15

# ---------------- experiment ----------------
epath = joinpath(@__DIR__, "..", "data", mp_exp_file)
isfile(epath) || error("File not found: " * epath)
dfe = CSV.read(epath, DataFrame)
hdr = names(dfe)
function ecol(name)
    j = findfirst(x -> normh(x) == normh(name), hdr)
    j === nothing && return nothing
    return tofloat(dfe[!, j])
end
te = ecol("t"); te === nothing && error("Column \"t\" not found. Columns: " * string(hdr))
qe = ecol("Heat Load"); qe === nothing && (qe = ecol("Q"))
t_on = if mp_t_onset === :auto
    qe === nothing && error("Heat load column not found; set mp_t_onset by hand")
    k = findfirst(i -> isfinite(qe[i]) && qe[i] > mp_onset_W, eachindex(qe))
    k === nothing && error("The heat load never exceeds mp_onset_W")
    te[k]
else
    Float64(mp_t_onset)
end
println("Experiment: heating onset t = ", round(t_on, digits = 2), " s")
function exp_mean(names)                  # mean of the available columns
    cs = [c for c in (ecol(n) for n in names) if c !== nothing]
    isempty(cs) && return nothing
    return [mean(filter(isfinite, [c[i] for c in cs])) for i in eachindex(te)]
end
ok_e = findall(isfinite, te)

# ---------------- model ----------------
mutable struct Res
    label::String; summ::DataFrame; prof::DataFrame; idx::Dict{Tuple{String,Float64},Vector{Int}}; tp::Vector{Float64}
end
models = Res[]
for (lab, folder) in mp_models
    d = joinpath(@__DIR__, "..", "res", folder)
    if !(isfile(joinpath(d, "summary.csv")) && isfile(joinpath(d, "profiles.csv")))
        @warn "Model results not found, skipped: " * d
        continue
    end
    summ = CSV.read(joinpath(d, "summary.csv"), DataFrame)
    prof = CSV.read(joinpath(d, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :T_wall])
    idx = Dict{Tuple{String,Float64},Vector{Int}}()
    for i in 1:nrow(prof)
        push!(get!(idx, (String(prof.comp[i]), Float64(prof.t[i])), Int[]), i)
    end
    push!(models, Res(lab, summ, prof, idx, sort(unique(prof.t))))
end
isempty(models) && error("No model results found. Run transient_flat.jl with stage = :long first.")
function cells(r::Res, comp, τ)
    t = r.tp[argmin(abs.(r.tp .- (mp_t_heat + τ)))]
    haskey(r.idx, (comp, Float64(t))) || return nothing
    return sort(r.prof[r.idx[(comp, Float64(t))], :], :cell)
end
function wall(r::Res, comp, τ)
    d = cells(r, comp, τ); d === nothing && return NaN
    return mean(d.T_wall) - K0
end
function wallpos(r::Res, comp, τ, pos)
    d = cells(r, comp, τ); d === nothing && return NaN
    N = nrow(d)
    return lin([(i - 0.5) / N for i in 1:N], Vector(d.T_wall), pos) - K0
end
function eva(r::Res, τ)
    v = wall(r, "eca", τ)
    return isnan(v) ? lin(r.summ.t, Vector(r.summ.T_ec), mp_t_heat + τ) - K0 : v
end
tmax(r::Res) = r.summ.t[end] - mp_t_heat

# sensors: (label, experiment columns (averaged), model function, colour)
G_T = [
  ("HB_EC_RES", [("HB", ["HB"], (r, τ) -> wall(r, "hb", τ), :red),
                 ("EC", ["EC1", "EC2"], (r, τ) -> eva(r, τ), :DarkOrange),
                 ("RES", ["RES1", "RES2"], (r, τ) -> wall(r, "cc", τ), :DarkGreen)]),
  ("VL1-3",    [("VL1", ["VL1"], (r, τ) -> wallpos(r, "vl", τ, 1 / 6), :MediumBlue),
                ("VL2", ["VL2"], (r, τ) -> wallpos(r, "vl", τ, 1 / 2), :DodgerBlue),
                ("VL3", ["VL3"], (r, τ) -> wallpos(r, "vl", τ, 5 / 6), :DeepSkyBlue)]),
  ("CON1-3",   [("CON1", ["CON1"], (r, τ) -> wallpos(r, "cond", τ, 0.0), :Indigo),
                ("CON2", ["CON2"], (r, τ) -> wallpos(r, "cond", τ, 0.25), :DeepPink),
                ("CON3", ["CON3"], (r, τ) -> wallpos(r, "cond", τ, 0.5), :Teal)]),
]
mstyle = [:dash, :dashdot, :dot, :dashdotdot]

outdir = joinpath(@__DIR__, "..", "res", mp_outdir); mkpath(outdir)
function newplot(ylabel)
    p = plot(; xlabel = "Time since heating onset [min]", ylabel = ylabel, xlims = mp_xlim, xticks = mp_xlim[1]:mp_xtick:mp_xlim[2],
             size = mp_size, dpi = mp_dpi, guidefontsize = mp_font_label, tickfontsize = mp_font_tick, legendfontsize = mp_font_legend,
             legend = :outertopright, framestyle = :box, grid = true, left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
    return p
end
function save(p, name)
    savefig(p, joinpath(outdir, name * ".svg"))
    mp_png && savefig(p, joinpath(outdir, name * ".png"))
    println("Saved: ", joinpath(outdir, name * ".svg"))
end

# ---------------- temperature figures ----------------
tau_ex = (te .- t_on)
for (gname, sens) in G_T
    p = newplot("Temperature [°C]")
    for (lab, ecols, f, col) in sens
        y = exp_mean(ecols)
        if y === nothing
            @warn "Experiment column(s) $(ecols) not found for $(lab); experiment skipped"
        else
            k = [i for i in ok_e if isfinite(y[i]) && mp_xlim[1] * 60 - 5 <= tau_ex[i] <= mp_xlim[2] * 60 + 5]
            plot!(p, tau_ex[k] ./ 60, y[k]; label = lab * " (exp.)", color = col, lw = mp_lw, ls = :solid)
        end
        for (m, r) in enumerate(models)
            τs = collect(0.0:mp_dt:min(mp_xlim[2] * 60, tmax(r)))
            plot!(p, τs ./ 60, [f(r, τ) for τ in τs]; label = lab * " (" * r.label * ")", color = col, lw = mp_lw, ls = mstyle[mod1(m, length(mstyle))])
        end
    end
    mp_Tlim === :auto || plot!(p, ylims = mp_Tlim)
    save(p, "temp_" * gname)
end

# ---------------- pressure (low-pass filtered) ----------------
function lp_filter(t, y)
    mp_P_filter || return y
    ok = isfinite.(t) .& isfinite.(y)
    count(ok) < 8 && return y
    tt, yy = t[ok], y[ok]
    tg = collect(tt[1]:mp_P_dt:tt[end])
    length(tg) < 8 && return y
    yg = [lin(tt, yy, x) for x in tg]
    n = length(yg); pad = min(100, n - 1)
    ext = [2 * yg[1] .- yg[pad+1:-1:2]; yg; 2 * yg[end] .- yg[end-1:-1:end-pad]]
    h = digitalfilter(Lowpass(mp_P_hz / (0.5 / mp_P_dt)), Butterworth(mp_P_order))
    zg = filtfilt(h, ext)[pad+1:pad+n]
    z = copy(y)
    z[ok] = [lin(tg, zg, x) for x in tt]
    return z
end
Pe_raw = ecol(mp_P_exp_col)
Pe = Pe_raw === nothing ? nothing : lp_filter(te, Pe_raw .* mp_P_scale)
Pe === nothing && @warn "Experiment pressure column '$(mp_P_exp_col)' not found; only the model is plotted"
ptag = mp_P_filter ? " (LP " * string(mp_P_hz) * " Hz)" : ""
ylab = (mp_P_rel ? "P − P(0)" : "Pressure") * " [kPa]"
function exp_curve()
    Pe === nothing && return nothing
    k = [i for i in ok_e if isfinite(Pe[i]) && mp_xlim[1] * 60 - 5 <= tau_ex[i] <= mp_xlim[2] * 60 + 5]
    y = Pe[k]
    if mp_P_rel
        i0 = argmin(abs.(tau_ex[ok_e]))
        y = y .- Pe[ok_e[i0]]
    end
    return tau_ex[k] ./ 60, y
end
function model_curve(r::Res)
    τs = collect(0.0:mp_dt:min(mp_xlim[2] * 60, tmax(r)))
    tt = r.summ.t
    P = lp_filter(tt, Vector(r.summ.P_cc) ./ 1e3)
    y = [lin(tt, P, mp_t_heat + τ) for τ in τs]
    mp_P_rel && (y = y .- lin(tt, P, mp_t_heat))
    return τs ./ 60, y
end
function pressure_plot(rs, name)
    p = newplot(ylab)
    ec = exp_curve()
    ec === nothing || plot!(p, ec[1], ec[2]; label = mp_P_exp_col * " exp." * ptag, color = :black, lw = mp_lw)
    for (m, r) in enumerate(rs)
        x, y = model_curve(r)
        plot!(p, x, y; label = "P_cc " * r.label * ptag, color = [:crimson, :royalblue, :seagreen, :purple][mod1(m, 4)], lw = mp_lw, ls = :dash)
    end
    mp_Plim === :auto || plot!(p, ylims = mp_Plim)
    save(p, name)
end
pressure_plot(models, "pressure_RES")
if length(models) > 1
    for (m, r) in enumerate(models)
        pressure_plot([r], "pressure_RES_" * replace(lowercase(r.label), r"[^a-z0-9]+" => "_"))
    end
end
println("Done. Output folder: ", outdir)