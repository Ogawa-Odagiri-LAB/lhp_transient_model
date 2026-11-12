# plot_heat_budget.jl — loop heat budget and condenser two-phase length of a finished run (SVG, one file per figure)
#
# Place in plot/ and run from the project root:  julia plot/plot_heat_budget.jl   (edit the settings below; no arguments)
#   plot/ must be a sibling of res/ (paths are resolved from this file as ../res)
# Input : res/<hb_result>/budget.csv (written by transient_flat.jl via heat_budget_calc.jl), profiles.csv, geometry.csv
# Output: res/<hb_outdir>/  heat_budget.svg        Q_apply, Q_plate, ambient loss (VL+LL), housing loss, storage rate, residual
#                           heat_storage.svg       storage rate by part (HB, EC case, CC case, tube walls, fluid, wick, plate/housing)
#                           cond_2ph_length.svg    two-phase length in the condenser vs time, with CON1-CON5 positions as horizontal lines
#                           cond_x_map.svg         quality x along the condenser (position vs time), CON1-CON5 marked
#                           housing_loss.svg       wall → insulation / PEEK heat of HB, EC case, CC case (from summary.csv)
#                           cond_amb_hyp.svg       hypothetical ambient heat gain of the condenser tube (not in the model) vs Q_apply
# CON positions (as in transient_flat.jl): CON1 = 0, CON2 = 0.25, CON3 = 0.5, CON4 = 0.75, CON5 = 1.0 of the condenser length (inlet -> outlet).
# The experiment has CON1-CON5 only (no CON6).

# ---------------------------- settings ----------------------------
hb_result   = "transient_flat_ammonia_v9a"          # folder under res/ (:long = 816 s run; for the :long2h run use "transient_flat_ammonia_2h_v9a")
hb_outdir   = "heat_budget"
hb_xlim     = :auto                                 # time axis [s since the heating start], e.g. (0, 7000) or (0, 300)
hb_logx     = false                                 # true: logarithmic time axis (shows both the start-up and the 2 h)
hb_font_label, hb_font_tick, hb_font_legend = 14, 12, 11
hb_size     = (1100, 650)
hb_lw       = 2.5
hb_dpi      = 300
hb_png      = true
hb_con_pos  = [("CON1", 0.0), ("CON2", 0.25), ("CON3", 0.5), ("CON4", 0.75), ("CON5", 1.0)]
# ------------------------------------------------------------------

using CSV, DataFrames, Plots
gr()
const HB_DIR = joinpath(@__DIR__, "..", "res", hb_result)
const HB_OUT = joinpath(@__DIR__, "..", "res", hb_outdir)
mkpath(HB_OUT)
isfile(joinpath(HB_DIR, "budget.csv")) || error("budget.csv が見つかりません: $(HB_DIR)\n  → transient_flat.jl (新しい版) をもう一度実行して budget.csv を作ってください")
b = CSV.read(joinpath(HB_DIR, "budget.csv"), DataFrame)
prof = CSV.read(joinpath(HB_DIR, "profiles.csv"), DataFrame)
geo = CSV.read(joinpath(HB_DIR, "geometry.csv"), DataFrame)
sm = CSV.read(joinpath(HB_DIR, "summary.csv"), DataFrame)
nrow(sm) == nrow(b) || error("summary.csv と budget.csv の行数が違います (別の計算の結果?)")
Lc = geo[geo.comp .== "cond", :L][1]
tau = b.tau
m = tau .>= 0
xs = hb_logx ? (v -> max.(v, 1.0)) : identity
xax = hb_logx ? :log10 : :identity
function base(; kw...)
    plot(; size = hb_size, dpi = hb_dpi, framestyle = :box, grid = true, gridalpha = 0.25,
         xguidefontsize = hb_font_label, yguidefontsize = hb_font_label, tickfontsize = hb_font_tick, legendfontsize = hb_font_legend,
         legend = :best, xscale = xax, left_margin = 6Plots.mm, bottom_margin = 5Plots.mm, kw...)
end
xl() = hb_xlim === :auto ? (hb_logx ? (1.0, maximum(tau)) : (0.0, maximum(tau))) : hb_xlim
function save(p, name)
    savefig(p, joinpath(HB_OUT, name * ".svg"))
    hb_png && savefig(p, joinpath(HB_OUT, name * ".png"))
    println("保存: ", joinpath(HB_OUT, name * ".svg"))
end

# 1. heat budget
p1 = base(xlabel = "Time since heating start [s]", ylabel = "Heat [W]", xlims = xl())
plot!(p1, xs(tau[m]), b.Q_apply[m], lw = hb_lw, color = :black, label = "Q_apply")
plot!(p1, xs(tau[m]), b.Q_plate[m], lw = hb_lw, color = :royalblue, label = "Q_plate (plate → cold plate)")
plot!(p1, xs(tau[m]), (b.Q_amb_vl .+ b.Q_amb_ll)[m], lw = hb_lw, color = :orange, label = "Q_amb (VL + LL)")
plot!(p1, xs(tau[m]), b.Q_house[m], lw = hb_lw, color = :seagreen, label = "Q_house (insulation, PEEK)")
plot!(p1, xs(tau[m]), b.dE_total[m], lw = hb_lw, color = :firebrick, label = "dE/dt (storage)")
plot!(p1, xs(tau[m]), b.resid[m], lw = 1.5, ls = :dash, color = :gray, label = "residual")
save(p1, "heat_budget")

# 2. storage by part
p2 = base(xlabel = "Time since heating start [s]", ylabel = "Storage rate [W]", xlims = xl())
for (col, lab, clr) in ((:dE_HB, "HB", :black), (:dE_EC, "EC case", :royalblue), (:dE_CC, "CC case", :orange), (:dE_wall, "tube walls", :seagreen),
                        (:dE_fluid, "fluid", :firebrick), (:dE_wick, "wick", :purple), (:dE_plate, "plate + cold plate + housing", :gray))
    plot!(p2, xs(tau[m]), b[m, col], lw = hb_lw, color = clr, label = lab)
end
save(p2, "heat_storage")

# 3. two-phase length with CON positions
p3 = base(xlabel = "Time since heating start [s]", ylabel = "Condenser length [m]", xlims = xl(), ylims = (0, Lc * 1.02), legend = :topright)
plot!(p3, xs(tau[m]), b.L_2ph[m], lw = hb_lw, color = :firebrick, label = "two-phase length (0.001 < x < 0.999)")
plot!(p3, xs(tau[m]), b.L_vap_eq[m], lw = hb_lw, color = :royalblue, ls = :dash, label = "equivalent vapor length (Σ dl·x)")
plot!(p3, xs(tau[m]), (b.x_front .* Lc)[m], lw = hb_lw, color = :black, ls = :dot, label = "vapor front (first cell with x < 0.5)")
xr = xl()
for (name, s) in hb_con_pos
    hline!(p3, [s * Lc], color = :gray, ls = :dashdot, lw = 1, label = "")
    annotate!(p3, xr[2], s * Lc, text(" " * name, hb_font_tick, :left, :gray30))
end
plot!(p3, right_margin = 12Plots.mm)
save(p3, "cond_2ph_length")

# 4. x map (position vs time)
c = prof[prof.comp .== "cond", :]
ts = sort(unique(c.t))
t_heat_guess = ts[1] - tau[1]            # t = tau + t_heat; budget.csv has the same time grid as profiles.csv
nC = maximum(c.cell)
X = fill(NaN, nC, length(ts))
let idx = Dict(t => j for (j, t) in enumerate(ts))
    for r in eachrow(c)
        X[r.cell, idx[r.t]] = clamp(r.x, 0.0, 1.0)
    end
end
tau_map = ts .- t_heat_guess
mm = tau_map .>= 0
pos = [(i - 0.5) / nC * Lc for i in 1:nC]
p4 = base(xlabel = "Time since heating start [s]", ylabel = "Position from condenser inlet [m]", xlims = xl(), ylims = (0, Lc), legend = false,
          right_margin = 22Plots.mm)
heatmap!(p4, xs(tau_map[mm]), pos, X[:, mm], c = :turbo, clims = (0, 1), colorbar_title = "quality x [-]")
for (name, s) in hb_con_pos
    hline!(p4, [s * Lc], color = :white, ls = :dash, lw = 1.2, label = "")
    annotate!(p4, xr[1], s * Lc, text(name * " ", hb_font_tick, :right, :black))
end
save(p4, "cond_x_map")

# 5. hypothetical ambient gain of the condenser tube
p5 = base(xlabel = "Time since heating start [s]", ylabel = "Heat [W]", xlims = xl())
plot!(p5, xs(tau[m]), b.Q_cond_amb_hyp[m], lw = hb_lw, color = :orange, label = "hypothetical ambient gain of the condenser tube (wrapped like VL)")
plot!(p5, xs(tau[m]), (0.01 .* b.Q_apply)[m], lw = 1.5, ls = :dash, color = :gray, label = "1 % of Q_apply")
save(p5, "cond_amb_hyp")

# 6. heat from the walls into the insulation / PEEK (+ = lost by the wall). Q_hi_cc = CC case → CC-side insulation node (curved surface only, G ≈ 4.5 mW/K in series with the outer surface).
p6 = base(xlabel = "Time since heating start [s]", ylabel = "Heat [W]", xlims = xl())
for (col, lab, clr) in ((:Q_hi_hb, "HB → insulation", :black), (:Q_peek, "HB → PEEK", :royalblue), (:Q_hi_ev, "EC case → insulation", :orange), (:Q_hi_cc, "CC case → insulation", :firebrick))
    plot!(p6, xs(tau[m]), sm[m, col], lw = hb_lw, color = clr, label = lab)
end
save(p6, "housing_loss")