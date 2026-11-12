# plot_exp_overview.jl — plot raw experimental data over a long run (temperatures on the left axis, heat load on the right axis)
#
# Run from the project root:  julia scripts/plot_exp_overview.jl   (edit the settings below; no arguments)
# Input : data/<ep_file>  (same format as exp_5W_startup.csv: column "t" in seconds, the logger's column names such as
#         "101_HB", "102_EC1", ..., "118_coldplate", "Heat Load" are accepted; a leading number and "_" are ignored)
# Output: res/exp_overview/<ep_out>.png  (also shown on screen)

# ---------------------------- settings ----------------------------
ep_file    = "exp_all.csv"      # file under data/ (EDIT: the long-run file)
ep_out     = "exp_overview"            # output file name (without extension)
ep_t_zero  = 0.0                       # time origin [s] (subtracted from t). e.g. 303.8 to count from the start of heating
ep_xlim    = (0.0, 840.0)              # x range [min]
ep_Tlim    = (0.0, 60.0)               # left axis [°C]
ep_Qlim    = (0.0, 200.0)              # right axis [W]
ep_size    = (1200, 600)               # figure size [pt] (the pixel size of the PNG is size × ep_dpi / 100)
ep_dpi     = 300                       # PNG resolution (default is 100). 300 gives 3600 × 1800 px at size (1200, 600)
ep_font_label  = 14                    # axis label font size [pt]
ep_font_tick   = 12                    # tick label font size [pt]
ep_font_legend = 12                    # legend font size [pt]
ep_legend_cols = 3                     # number of columns of the temperature legend (3 columns x 3 rows for 9 series)
ep_legend_pos  = :topleft              # legend position (e.g. :topleft, :top, :outertop)
ep_legend_down = false                 # false: the legend fills across (HB, EVA1, VL1 in the first ROW). true: fills down (HB, EVA1, VL1 in the first COLUMN).
                                       #   If the order on the figure is not what you want, flip this flag.
ep_lw          = 1.2                   # line width of the temperatures
ep_lw_heat     = 2.0                   # line width of the heat load
ep_vector  = true                      # true: also save a vector file (PDF and SVG; always sharp, good for papers and slides)
# (label shown in the legend, [accepted column names, first match is used])
ep_series  = [("HB",         ["HB"]),
              ("EVA1",       ["EVA1", "EC1"]),
              ("VL1",        ["VL1"]),
              ("CON1",       ["CON1"]),
              ("CON2",       ["CON2"]),
              ("CON3",       ["CON3"]),
              ("LL1",        ["LL1"]),
              ("RES1",       ["RES1"]),
              ("Cold plate", ["coldplate", "cold_plate", "ColdPlate"])]
ep_heat    = ("Heat Load", ["Heat Load", "Q"])

# ---------------------------- code ----------------------------
using Pkg
Pkg.activate(".")
using CSV, DataFrames, Plots

normh(x) = lowercase(replace(replace(strip(replace(string(x), "﻿" => "")), r"^\d+_" => ""), r"[ _]" => ""))
tofloat(v) = [x isa Number ? Float64(x) : (y = tryparse(Float64, strip(string(x))); y === nothing ? NaN : y) for x in v]

path = joinpath(@__DIR__, "..", "data", ep_file)
isfile(path) || error("File not found: " * path)
df = CSV.read(path, DataFrame)
hdr = names(df)
function getcol(cands)
    for c in cands
        j = findfirst(x -> normh(x) == normh(c), hdr)
        j !== nothing && return tofloat(df[!, j]), hdr[j]
    end
    return nothing, nothing
end
t, _ = getcol(["t"])
t === nothing && error("Column \"t\" not found. Columns: " * string(hdr))
x = (t .- ep_t_zero) ./ 60                       # [min]
println("Read ", length(t), " rows: t = ", round(minimum(filter(isfinite, t)), digits = 1), " to ", round(maximum(filter(isfinite, t)), digits = 1),
        " s (", round((maximum(filter(isfinite, t)) - minimum(filter(isfinite, t))) / 60, digits = 1), " min)")
palette = [:red, :DeepPink, :orange, :MediumBlue, :DodgerBlue, :DeepSkyBlue, :Indigo, :DarkGreen, :gray]
p = plot(xlabel = "Time [min]", ylabel = "Temperature [°C]", xlims = ep_xlim, ylims = ep_Tlim, legend = ep_legend_pos, legend_columns = ep_legend_cols,
         size = ep_size, dpi = ep_dpi, guidefontsize = ep_font_label, tickfontsize = ep_font_tick, legendfontsize = ep_font_legend, right_margin = 18Plots.mm, left_margin = 6Plots.mm, bottom_margin = 6Plots.mm, grid = true)
nser = length(ep_series)
nrow_leg = cld(nser, ep_legend_cols)
# plotting order = legend order. ep_legend_down = true: arrange the entries so that the legend reads down each column
order = ep_legend_down ? [(c - 1) * nrow_leg + r for r in 1:nrow_leg for c in 1:ep_legend_cols if (c - 1) * nrow_leg + r <= nser] : collect(1:nser)
for k in order
    (lab, cands) = ep_series[k]
    y, used = getcol(cands)
    if y === nothing
        @warn "Column not found for $(lab) (tried $(cands)); skipped. Columns: $(hdr)"
        continue
    end
    plot!(p, x, y; label = lab, lw = ep_lw, color = palette[mod1(k, length(palette))])
end
q, usedq = getcol(ep_heat[2])
if q === nothing
    @warn "Heat load column not found (tried $(ep_heat[2]))"
else
    p2 = twinx(p)
    plot!(p2, x, q; label = ep_heat[1], lw = ep_lw_heat, color = :black, ls = :dash, ylabel = "Heat load [W]", ylims = ep_Qlim, xlims = ep_xlim, legend = :topright, grid = false, guidefontsize = ep_font_label, tickfontsize = ep_font_tick, legendfontsize = ep_font_legend)
end
out = joinpath(@__DIR__, "..", "res", "exp_overview")
mkpath(out)
savefig(p, joinpath(out, ep_out * ".png"))
if ep_vector
    savefig(p, joinpath(out, ep_out * ".pdf"))
    savefig(p, joinpath(out, ep_out * ".svg"))
end
display(p)
println("Saved: ", joinpath(out, ep_out * ".png"), ep_vector ? " (+ .pdf, .svg)" : "")