# plot_initial_temperature_exp.jl — heat map of the EXPERIMENTAL temperature distribution at t = 0 (start of the data, before heating) on a schematic of the LHP
#
# Run from the project root:  julia scripts/plot_initial_temperature_exp.jl   (edit the settings below; no arguments)
# Input : data/<ie_file>  (same format as plot_exp_onset.jl: column "t" [s], logger column names accepted with or without the leading number)
# Method: the temperature of each sensor at t = 0 is read from the data (average over the ie_avg_s seconds before t = 0, to suppress noise),
#         then LINEARLY INTERPOLATED along each tube between the sensors (constant beyond the end sensors). Nothing is simulated.
#         Sensor positions (fraction of the tube length, measured from the evaporator side for VL, from the vapor-line side for CON,
#         from the condenser side for LL) are the same as in the model:  VL1-3 = 1/6, 1/2, 5/6;  CON1-5 = 0, 1/4, 1/2, 3/4, 1;
#         LL1-3 = 0.10, 0.65, 1.11 m of the 1.15 m liquid line.  HB, RES (CC dome) are single values. groove and wick are not measured (grey).
# Output: res/exp_initial/<ie_out>.png / .pdf / .svg
# Drawing: same schematic as plot_initial_temperature.jl; lengths are not to scale.

# ---------------------------- settings ----------------------------
ie_file    = "exp_5W_startup.csv"                # under data/
ie_out     = "initial_temperature_map_exp"     # output file name (without extension)
ie_t0      = :start                            # t = 0 [s in the file]. :start = the first sample of the file (true t = 0, before heating); :auto = heating onset (first time the heat load exceeds ie_onset_W); or a number
ie_onset_W = 1.0
ie_avg_s   = 10.0                              # average the sensor readings over this many seconds (:start = the first ie_avg_s seconds of the file; otherwise the ie_avg_s seconds before t = 0)
it_Tlim    = :auto                             # colour range [°C], e.g. (20.0, 23.0); :auto = min to max of the sensors
it_cmap    = :turbo                            # colour map (e.g. :turbo, :thermal, :viridis, :coolwarm)
it_sensors = true                              # true: mark the thermocouple positions
it_values  = true                              # true: write the measured temperature [°C] next to each sensor label (1 decimal)
it_size    = (1200, 1200)                      # square figure
it_dpi     = 300
it_font    = 20                                # base font size [pt]
it_lw      = 14                                # drawn tube width [pt]
it_vector  = true                              # also save PDF and SVG

# ---------------------------- code ----------------------------
using Pkg
Pkg.activate(".")
using CSV, DataFrames, Plots

normh(x) = lowercase(replace(replace(strip(replace(string(x), "\ufeff" => "")), r"^\d+_" => ""), r"[ _\"]" => ""))
tofloat(v) = [x isa Number ? Float64(x) : (y = tryparse(Float64, strip(string(x))); y === nothing ? NaN : y) for x in v]
epath = joinpath(@__DIR__, "..", "data", ie_file)
isfile(epath) || error("File not found: " * epath)
dfe = CSV.read(epath, DataFrame)
hdr = names(dfe)
function ecol(name)
    j = findfirst(x -> normh(x) == normh(name), hdr)
    return j === nothing ? nothing : tofloat(dfe[!, j])
end
te = ecol("t"); te === nothing && error("Column \"t\" not found. Columns: " * string(hdr))
qe = ecol("Heat Load"); qe === nothing && (qe = ecol("Q"))
t0 = if ie_t0 === :start
    minimum(filter(isfinite, te))
elseif ie_t0 === :auto
    qe === nothing && error("Heat load column not found; set ie_t0 by hand")
    k = findfirst(i -> isfinite(qe[i]) && qe[i] > ie_onset_W, eachindex(qe))
    k === nothing && error("The heat load never exceeds ie_onset_W")
    te[k]
else
    Float64(ie_t0)
end
println("t = 0 at t = ", round(t0, digits = 2), " s of the file (averaging window ", ie_avg_s, " s before)")
win = ie_t0 === :start ? findall(i -> isfinite(te[i]) && t0 <= te[i] <= t0 + ie_avg_s, eachindex(te)) : findall(i -> isfinite(te[i]) && t0 - ie_avg_s <= te[i] <= t0, eachindex(te))
isempty(win) && error("No data in the averaging window")
function sensor(names...)             # mean of the available columns over the window [°C]
    vals = Float64[]
    for n in names
        c = ecol(n)
        c === nothing && continue
        v = filter(isfinite, c[win])
        isempty(v) || push!(vals, sum(v) / length(v))
    end
    isempty(vals) && (@warn "Sensor column(s) $(names) not found"; return NaN)
    return sum(vals) / length(vals)
end
Tvl  = [sensor("VL1"), sensor("VL2"), sensor("VL3")];                       pvl  = [1 / 6, 1 / 2, 5 / 6]
Tcon = [sensor("CON1"), sensor("CON2"), sensor("CON3"), sensor("CON4"), sensor("CON5")]; pcon = [0.0, 0.25, 0.5, 0.75, 1.0]
Tll  = [sensor("LL1"), sensor("LL2"), sensor("LL3")];                       pll  = [0.10, 0.65, 1.11] ./ 1.15
function interp_line(pos, T, N)       # linear between sensors, constant beyond the ends; value at the N cell centres
    k = [i for i in eachindex(T) if isfinite(T[i])]
    isempty(k) && return fill(NaN, N)
    c = pos[k]; y = T[k]
    f(s) = s <= c[1] ? y[1] : (s >= c[end] ? y[end] : (m = searchsortedlast(c, s); y[m] + (y[m+1] - y[m]) * (s - c[m]) / (c[m+1] - c[m])))
    return [f((i - 0.5) / N) for i in 1:N]
end
x_vl, x_cond, x_ll = interp_line(pvl, Tvl, 20), interp_line(pcon, Tcon, 50), interp_line(pll, Tll, 20)
T_hb, T_grv, T_wick, T_cc = sensor("HB"), NaN, NaN, sensor("RES1", "RES2")
f_v_cc = 0.429                         # drawn CC vapor fraction (initial value of the model; not measured)
const K0 = 273.15
resdir = joinpath(@__DIR__, "..", "res", "exp_initial"); mkpath(resdir)
ie_pos = Dict("VL" => (pvl, Tvl), "CON" => (pcon, Tcon), "LL" => (pll, Tll))

allT = filter(isfinite, vcat(x_vl, x_cond, x_ll, [T_hb, T_grv, T_wick, T_cc]))
isempty(allT) && error("No sensor temperatures found at t = " * string(t0))
Tmin, Tmax = it_Tlim === :auto ? (minimum(allT), maximum(allT)) : (Float64(it_Tlim[1]), Float64(it_Tlim[2]))
Tmax == Tmin && (Tmax = Tmin + 1e-3)
cg = cgrad(it_cmap)
xcol(T) = isnan(T) ? RGB(0.7, 0.7, 0.7) : cg[clamp((T - Tmin) / (Tmax - Tmin), 0.0, 1.0)]
println("Sensor temperature range [°C]: ", round(minimum(allT), digits = 3), " to ", round(maximum(allT), digits = 3), "  (colour scale ", round(Tmin, digits = 3), " to ", round(Tmax, digits = 3), ")")

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
         title = "Experiment: temperature distribution at t = 0 (interpolated between sensors)", titlefontsize = it_font - 4)
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
    function sensT(nm, s)            # measured value of this sensor
        grp = startswith(nm, "VL") ? "VL" : (startswith(nm, "CON") ? "CON" : "LL")
        pos, T = ie_pos[grp]
        return T[argmin(abs.(pos .- s))]
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
annotate!(p, 10.5, -1.0, text("measured temperature [°C]", it_font - 1, :black))
annotate!(p, -0.4, -1.2, text("grey = not measured (groove, wick)", it_font - 4, :gray30, :left))

savefig(p, joinpath(resdir, ie_out * ".png"))
if it_vector
    savefig(p, joinpath(resdir, ie_out * ".pdf"))
    savefig(p, joinpath(resdir, ie_out * ".svg"))
end
display(p)
println("Saved: ", joinpath(resdir, ie_out * ".png"), it_vector ? " (+ .pdf, .svg)" : "")