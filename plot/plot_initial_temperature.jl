# plot_initial_temperature.jl — plot the initial liquid/vapor distribution (t = 0) on a schematic of the LHP
#
# Run from the project root:  julia scripts/plot_initial_phase.jl   (edit the settings below; no arguments)
# Input : res/<it_result>/profiles.csv (vapor quality x of every cell at the first saved time) and summary.csv (CC vapor fraction f_v_cc)
#         -> run transient_flat.jl first (any stage); the initial state is the same for the same initial-condition settings.
# Output: res/<it_result>/initial_phase_map.png / .pdf / .svg
# Drawing: vapor line (VL, 20 cells), condenser (CON, 50 cells, serpentine), liquid line (LL, 20 cells), groove, primary wick, CC (half cylinder).
#          Each tube cell is coloured by its vapor quality x (blue = liquid, red = vapor, in between = two-phase).
#          Cell numbers run from the evaporator side (VL1 ... ), from the vapor-line side (CON1 ...) and from the condenser side (LL1 ... LL3 near the CC).
#          This is a schematic: lengths are not to scale (every cell of a tube has the same drawn length).

# plot_initial_temperature.jl — heat map of the initial temperature distribution (t = 0) on a schematic of the LHP
#
# Run from the project root:  julia scripts/plot_initial_temperature.jl   (edit the settings below; no arguments)
# Input : res/<it_result>/profiles.csv (temperature of every cell at the first saved time) and summary.csv (CC vapor fraction f_v_cc)
#         -> run transient_flat.jl first (any stage); the initial state is the same for the same initial-condition settings.
# Output: res/<it_result>/initial_temperature_map.png / .pdf / .svg
# Drawing: same schematic as plot_initial_phase.jl (vapor line 20 cells, condenser 50 cells, liquid line 20 cells, HB, groove, primary wick, CC).
#          Each cell is coloured by its temperature. it_field = :T_wall (wall temperature, what the thermocouples read) or :T (fluid).
#          Components whose temperature is not in profiles.csv are drawn grey. Lengths are not to scale.

# ---------------------------- settings ----------------------------
it_result  = "transient_flat_ammonia_v9a"      # result folder under res/ (the initial state is read from its first saved time)
it_out     = "initial_temperature_map"         # output file name (without extension)
it_field   = :T_wall                           # :T_wall (wall) or :T (fluid)
it_Tlim    = :auto                             # colour range [°C], e.g. (20.0, 23.0); :auto = min to max of the plotted cells
it_cmap    = :turbo                            # colour map (e.g. :turbo, :thermal, :viridis, :coolwarm)
it_sensors = true                              # true: mark the thermocouple positions (VL1-3, CON1-5, LL1-3)
it_values  = true                              # true: write the temperature [°C] next to each sensor label (1 decimal)
it_size    = (1200, 1200)                      # square figure
it_dpi     = 300
it_font    = 20                                # base font size [pt] (all labels scale with it; sensor labels = base-2, small labels = base-1)
it_lw      = 14                                # drawn tube width [pt]
it_vector  = true                              # also save PDF and SVG

# ---------------------------- code ----------------------------
using Pkg
Pkg.activate(".")
using CSV, DataFrames, Plots

resdir = isabspath(it_result) ? it_result : joinpath(@__DIR__, "..", "res", it_result)
isfile(joinpath(resdir, "profiles.csv")) || error("profiles.csv not found in " * resdir * "  (run transient_flat.jl first)")
prof = CSV.read(joinpath(resdir, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :T, :T_wall])
summ = CSV.read(joinpath(resdir, "summary.csv"), DataFrame)
t0 = minimum(prof.t)
println("Initial state read at t = ", t0, " s from ", resdir)

const K0 = 273.15
function cell_T(comp, N)            # temperature [°C] of the cells of one component (NaN where not available)
    d = prof[(prof.t .== t0) .& (String.(prof.comp) .== comp), :]
    T = fill(NaN, N)
    for r in eachrow(d)
        v = r[it_field]
        (ismissing(v) || !isfinite(v)) && (v = r[it_field === :T ? :T_wall : :T])
        (!ismissing(v) && isfinite(v) && 1 <= r.cell <= N) && (T[r.cell] = v - K0)
    end
    return T
end
meanT(comp) = (T = cell_T(comp, 200); v = filter(isfinite, T); isempty(v) ? NaN : sum(v) / length(v))
x_vl, x_cond, x_ll = cell_T("vl", 20), cell_T("cond", 50), cell_T("ll", 20)
T_hb, T_grv, T_wick, T_cc = meanT("hb"), meanT("grv"), meanT("wick"), meanT("cc")
f_v_cc = "f_v_cc" in names(summ) && isfinite(summ.f_v_cc[1]) ? summ.f_v_cc[1] : NaN
isnan(f_v_cc) && @warn "f_v_cc not found in summary.csv: the CC is drawn with f_v = 0.43 (default)"
isnan(f_v_cc) && (f_v_cc = 0.43)

allT = filter(isfinite, vcat(x_vl, x_cond, x_ll, [T_hb, T_grv, T_wick, T_cc]))
isempty(allT) && error("No temperatures found in profiles.csv at t = " * string(t0))
Tmin, Tmax = it_Tlim === :auto ? (minimum(allT), maximum(allT)) : (Float64(it_Tlim[1]), Float64(it_Tlim[2]))
Tmax == Tmin && (Tmax = Tmin + 1e-3)
cg = cgrad(it_cmap)
xcol(T) = isnan(T) ? RGB(0.7, 0.7, 0.7) : cg[clamp((T - Tmin) / (Tmax - Tmin), 0.0, 1.0)]
println("Temperature range [°C]: ", round(minimum(allT), digits = 3), " to ", round(maximum(allT), digits = 3), "  (colour scale ", round(Tmin, digits = 3), " to ", round(Tmax, digits = 3), ")")

# ---- tube paths (waypoints) and helpers ----
wp_vl = [(3.0, 2.25), (4.0, 2.25), (4.0, 6.2), (11.5, 6.2)]
ys = [5.4, 4.6, 3.8, 3.0, 2.2]
wp_cond = Tuple{Float64,Float64}[(11.5, 6.2)]
for (k, y) in enumerate(ys)
    xa, xb = isodd(k) ? (11.5, 8.0) : (8.0, 11.5)
    push!(wp_cond, (xa, y)); push!(wp_cond, (xb, y))
end
wp_ll = [(8.0, 2.2), (8.0, 0.8), (0.3, 0.8), (0.3, 3.3), (1.0, 3.3)]

seglen(a, b) = hypot(b[1] - a[1], b[2] - a[2])
function cumlen(wp)
    c = [0.0]
    for i in 2:length(wp)
        push!(c, c[end] + seglen(wp[i-1], wp[i]))
    end
    return c
end
function point_at(wp, s)            # s = arclength fraction 0..1
    c = cumlen(wp); L = c[end]; d = s * L
    i = clamp(searchsortedlast(c, d), 1, length(wp) - 1)
    f = (d - c[i]) / (c[i+1] - c[i])
    return (wp[i][1] + f * (wp[i+1][1] - wp[i][1]), wp[i][2] + f * (wp[i+1][2] - wp[i][2]))
end
function cell_polyline(wp, s0, s1)  # points from s0 to s1 including the corners in between
    c = cumlen(wp); L = c[end]
    pts = [point_at(wp, s0)]
    for i in 2:length(wp)-1
        s0 < c[i] / L < s1 && push!(pts, wp[i])
    end
    push!(pts, point_at(wp, s1))
    return pts
end
function draw_tube!(p, wp, xs; lw = it_lw)
    N = length(xs)
    for i in 1:N
        pts = cell_polyline(wp, (i - 1) / N, i / N)
        plot!(p, first.(pts), last.(pts); color = :black, lw = lw + 2, label = "", linecap = :butt)
    end
    for i in 1:N
        pts = cell_polyline(wp, (i - 1) / N, i / N)
        plot!(p, first.(pts), last.(pts); color = xcol(xs[i]), lw = lw, label = "", linecap = :butt)
    end
end

# ---- plot ----
p = plot(; aspect_ratio = :equal, framestyle = :none, legend = false, size = it_size, dpi = it_dpi, xlims = (-0.6, 12.4), ylims = (-2.9, 10.1),
         title = "Initial temperature distribution (t = " * string(round(t0, digits = 1)) * " s)", titlefontsize = it_font + 2)
draw_tube!(p, wp_vl, x_vl)
draw_tube!(p, wp_cond, x_cond)
draw_tube!(p, wp_ll, x_ll)

# heater block and evaporator (primary wick + groove)
hb = Shape([1.0, 3.0, 3.0, 1.0], [1.4, 1.4, 2.0, 2.0])
plot!(p, hb; fillcolor = xcol(T_hb), linecolor = :black, label = "")
annotate!(p, 2.0, 1.7, text("HB", it_font - 4, :black))
groove = Shape([1.0, 3.0, 3.0, 1.0], [2.0, 2.0, 2.5, 2.5])
plot!(p, groove; fillcolor = xcol(T_grv), linecolor = :black, label = "")
annotate!(p, 2.0, 2.25, text("groove", it_font - 6, :black))
wick = Shape([1.0, 3.0, 3.0, 1.0], [2.5, 2.5, 3.3, 3.3])
plot!(p, wick; fillcolor = xcol(T_wick), linecolor = :black, label = "")
annotate!(p, 2.0, 2.9, text("wick", it_font - 6, :black))

# CC: half cylinder (dome up, flat side down on the wick). Vapor cap at the top, liquid below.
cx, cy, R = 2.0, 3.3, 1.0
function vapor_theta(f)             # (θ − sin θ)/π = f  ->  θ [rad]
    lo, hi = 0.0, π
    for _ in 1:80
        mid = 0.5 * (lo + hi)
        (mid - sin(mid)) / π < f ? (lo = mid) : (hi = mid)
    end
    return 0.5 * (lo + hi)
end
θv = vapor_theta(clamp(f_v_cc, 0.0, 1.0))
φ = range(0, π, length = 121)
dome = Shape(cx .+ R .* cos.(φ), cy .+ R .* sin.(φ))
plot!(p, dome; fillcolor = xcol(T_cc), linecolor = :black, lw = 2, label = "")
if θv > 1e-3                         # vapor-liquid interface (outline only)
    φc = range(π / 2 - θv / 2, π / 2 + θv / 2, length = 61)
    plot!(p, cx .+ R .* cos.(φc), cy .+ R .* sin.(φc); color = :black, lw = 1, ls = :dash, label = "")
    plot!(p, [cx + R * cos(φc[1]), cx + R * cos(φc[end])], [cy + R * sin(φc[1]), cy + R * sin(φc[end])]; color = :black, lw = 1, ls = :dash, label = "")
end
isfinite(T_cc) && annotate!(p, cx, cy + 0.45, text(string(round(T_cc, digits = 1)) * " °C", it_font - 2, :black))
annotate!(p, cx, cy + R + 0.3, text("CC  (f_v = " * string(round(f_v_cc, digits = 3)) * ")", it_font, :black))

# tube labels
annotate!(p, 7.7, 6.85, text("Vapor line (20 cells)", it_font, :black))
annotate!(p, 9.75, 1.5, text("Condenser (50 cells)", it_font, :black))
annotate!(p, 4.6, 0.35, text("Liquid line (20 cells)", it_font, :black))
annotate!(p, 3.7, 3.4, text("→ VL1", it_font - 2, :gray30))
annotate!(p, 0.0, 2.0, text("LL → CC", it_font - 2, :gray30))

# thermocouple positions
if it_sensors
    sens = [("VL1", wp_vl, 1 / 6), ("VL2", wp_vl, 1 / 2), ("VL3", wp_vl, 5 / 6),
            ("CON1", wp_cond, 0.0), ("CON2", wp_cond, 0.25), ("CON3", wp_cond, 0.5), ("CON4", wp_cond, 0.75), ("CON5", wp_cond, 1.0),
            ("LL1", wp_ll, 0.10 / 1.15), ("LL2", wp_ll, 0.65 / 1.15), ("LL3", wp_ll, 1.11 / 1.15)]
    Tarr = Dict("VL" => x_vl, "CON" => x_cond, "LL" => x_ll)
    function sensT(nm, s)            # temperature at the sensor position (linear interpolation between cell centres)
        T = Tarr[startswith(nm, "VL") ? "VL" : (startswith(nm, "CON") ? "CON" : "LL")]
        N = length(T); k = [i for i in 1:N if isfinite(T[i])]
        isempty(k) && return NaN
        c = [(i - 0.5) / N for i in k]
        return s <= c[1] ? T[k[1]] : (s >= c[end] ? T[k[end]] : (m = searchsortedlast(c, s); T[k[m]] + (T[k[m+1]] - T[k[m]]) * (s - c[m]) / (c[m+1] - c[m])))
    end
    for (nm, wp, s) in sens
        (px, py) = point_at(wp, s)
        scatter!(p, [px], [py]; marker = :circle, ms = 6, color = :white, markerstrokecolor = :black, markerstrokewidth = 2, label = "")
        annotate!(p, px, py + (startswith(nm, "VL") || nm == "CON1" ? 0.35 : (startswith(nm, "CON") ? 0.0 : -0.3)) + 0.0, text(it_values ? nm * " " * string(round(sensT(nm, s), digits = 1)) : nm, it_font - 2, :black, startswith(nm, "CON") && nm != "CON1" ? :left : :center))
    end
end

# colour bar (temperature)
for k in 0:40
    xk = k / 40
    Tk = Tmin + xk * (Tmax - Tmin)
    plot!(p, Shape([9.0 + 3.0 * xk, 9.0 + 3.0 * (k + 1) / 40, 9.0 + 3.0 * (k + 1) / 40, 9.0 + 3.0 * xk], [-1.6, -1.6, -1.25, -1.25]); fillcolor = xcol(Tk), linecolor = xcol(Tk), label = "")
end
annotate!(p, 9.0, -1.95, text(string(round(Tmin, digits = 2)) * " °C", it_font - 2, :black, :left))
annotate!(p, 12.2, -1.95, text(string(round(Tmax, digits = 2)) * " °C", it_font - 2, :black, :right))
annotate!(p, 10.5, -1.0, text((it_field === :T ? "fluid" : "wall") * " temperature [°C]", it_font - 1, :black))
annotate!(p, -0.4, -1.2, text("grey = temperature not available in profiles.csv", it_font - 4, :gray30, :left))

savefig(p, joinpath(resdir, it_out * ".png"))
if it_vector
    savefig(p, joinpath(resdir, it_out * ".pdf"))
    savefig(p, joinpath(resdir, it_out * ".svg"))
end
display(p)
println("Saved: ", joinpath(resdir, it_out * ".png"), it_vector ? " (+ .pdf, .svg)" : "")