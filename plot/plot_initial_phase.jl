# plot_initial_phase.jl — plot the initial liquid/vapor distribution (t = 0) on a schematic of the LHP
#
# Run from the project root:  julia scripts/plot_initial_phase.jl   (edit the settings below; no arguments)
# Input : res/<ip_result>/profiles.csv (vapor quality x of every cell at the first saved time) and summary.csv (CC vapor fraction f_v_cc)
#         -> run transient_flat.jl first (any stage); the initial state is the same for the same initial-condition settings.
# Output: res/<ip_result>/initial_phase_map.png / .pdf / .svg
# Drawing: vapor line (VL, 20 cells), condenser (CON, 50 cells, serpentine), liquid line (LL, 20 cells), groove, primary wick, CC (half cylinder).
#          Each tube cell is coloured by its vapor quality x (blue = liquid, red = vapor, in between = two-phase).
#          Cell numbers run from the evaporator side (VL1 ... ), from the vapor-line side (CON1 ...) and from the condenser side (LL1 ... LL3 near the CC).
#          This is a schematic: lengths are not to scale (every cell of a tube has the same drawn length).

# ---------------------------- settings ----------------------------
ip_result  = "transient_flat_ammonia_v9a"      # result folder under res/ (the initial state is read from its first saved time)
ip_out     = "initial_phase_map"               # output file name (without extension)
ip_sensors = true                              # true: mark the thermocouple positions (VL1-3, CON1-5, LL1-3)
ip_size    = (1200, 1200)                      # square figure
ip_dpi     = 300
ip_font    = 20                                # base font size [pt] (all labels scale with it; sensor labels = base-2, small labels = base-1)
ip_lw      = 14                                # drawn tube width [pt]
ip_vector  = true                              # also save PDF and SVG

# ---------------------------- code ----------------------------
using Pkg
Pkg.activate(".")
using CSV, DataFrames, Plots

resdir = isabspath(ip_result) ? ip_result : joinpath(@__DIR__, "..", "res", ip_result)
isfile(joinpath(resdir, "profiles.csv")) || error("profiles.csv not found in " * resdir * "  (run transient_flat.jl first)")
prof = CSV.read(joinpath(resdir, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :x, :T])
summ = CSV.read(joinpath(resdir, "summary.csv"), DataFrame)
t0 = minimum(prof.t)
println("Initial state read at t = ", t0, " s from ", resdir)

function cell_x(comp, N)
    d = prof[(prof.t .== t0) .& (String.(prof.comp) .== comp), :]
    x = fill(NaN, N)
    for r in eachrow(d)
        1 <= r.cell <= N && (x[r.cell] = r.x)
    end
    return x
end
x_vl, x_cond, x_ll = cell_x("vl", 20), cell_x("cond", 50), cell_x("ll", 20)
x_grv = let v = cell_x("grv", 1); isnan(v[1]) ? 1.0 : v[1] end
f_v_cc = "f_v_cc" in names(summ) && isfinite(summ.f_v_cc[1]) ? summ.f_v_cc[1] : NaN
isnan(f_v_cc) && @warn "f_v_cc not found in summary.csv: the CC is drawn with f_v = 0.43 (default)"
isnan(f_v_cc) && (f_v_cc = 0.43)

# colour of vapor quality x: 0 = liquid (blue) ... 1 = vapor (red)
cg = cgrad([:royalblue, :mediumpurple, :crimson])
xcol(x) = isnan(x) ? RGB(0.7, 0.7, 0.7) : cg[clamp(x, 0.0, 1.0)]

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
function draw_tube!(p, wp, xs; lw = ip_lw)
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
p = plot(; aspect_ratio = :equal, framestyle = :none, legend = false, size = ip_size, dpi = ip_dpi, xlims = (-0.6, 12.4), ylims = (-2.9, 10.1),
         title = "Initial liquid/vapor distribution (t = " * string(round(t0, digits = 1)) * " s)", titlefontsize = ip_font + 2)
draw_tube!(p, wp_vl, x_vl)
draw_tube!(p, wp_cond, x_cond)
draw_tube!(p, wp_ll, x_ll)

# heater block and evaporator (primary wick + groove)
hb = Shape([1.0, 3.0, 3.0, 1.0], [1.4, 1.4, 2.0, 2.0])
plot!(p, hb; fillcolor = RGB(0.85, 0.85, 0.85), linecolor = :black, label = "")
annotate!(p, 2.0, 1.7, text("HB", ip_font - 4, :black))
groove = Shape([1.0, 3.0, 3.0, 1.0], [2.0, 2.0, 2.5, 2.5])
plot!(p, groove; fillcolor = xcol(x_grv), linecolor = :black, label = "")
annotate!(p, 2.0, 2.25, text("groove", ip_font - 6, :white))
wick = Shape([1.0, 3.0, 3.0, 1.0], [2.5, 2.5, 3.3, 3.3])
plot!(p, wick; fillcolor = xcol(0.0), fillalpha = 0.55, linecolor = :black, label = "")
annotate!(p, 2.0, 2.9, text("wick", ip_font - 6, :black))

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
plot!(p, dome; fillcolor = xcol(0.0), fillalpha = 0.55, linecolor = :black, lw = 2, label = "")
φc = range(π / 2 - θv / 2, π / 2 + θv / 2, length = 61)
if θv > 1e-3
    cap = Shape(cx .+ R .* cos.(φc), cy .+ R .* sin.(φc))
    plot!(p, cap; fillcolor = xcol(1.0), fillalpha = 0.85, linecolor = :black, lw = 1, label = "")
end
annotate!(p, cx, cy + 0.18, text("vapor", ip_font - 6, :white))
annotate!(p, cx, cy + 0.65 * R * cos(θv / 2) * 0.55 - 0.0, text("liquid", ip_font - 6, :black))
annotate!(p, cx, cy + R + 0.3, text("CC  (f_v = " * string(round(f_v_cc, digits = 3)) * ")", ip_font, :black))

# tube labels
annotate!(p, 7.7, 6.85, text("Vapor line (20 cells)", ip_font, :black))
annotate!(p, 9.75, 1.5, text("Condenser (50 cells)", ip_font, :black))
annotate!(p, 4.6, 0.35, text("Liquid line (20 cells)", ip_font, :black))
annotate!(p, 3.7, 3.4, text("→ VL1", ip_font - 2, :gray30))
annotate!(p, 0.0, 2.0, text("LL → CC", ip_font - 2, :gray30))

# thermocouple positions
if ip_sensors
    sens = [("VL1", wp_vl, 1 / 6), ("VL2", wp_vl, 1 / 2), ("VL3", wp_vl, 5 / 6),
            ("CON1", wp_cond, 0.0), ("CON2", wp_cond, 0.25), ("CON3", wp_cond, 0.5), ("CON4", wp_cond, 0.75), ("CON5", wp_cond, 1.0),
            ("LL1", wp_ll, 0.10 / 1.15), ("LL2", wp_ll, 0.65 / 1.15), ("LL3", wp_ll, 1.11 / 1.15)]
    for (nm, wp, s) in sens
        (px, py) = point_at(wp, s)
        scatter!(p, [px], [py]; marker = :circle, ms = 6, color = :white, markerstrokecolor = :black, markerstrokewidth = 2, label = "")
        annotate!(p, px, py + (startswith(nm, "VL") ? 0.35 : (startswith(nm, "CON") ? 0.0 : -0.3)) + 0.0, text(nm, ip_font - 2, :black, startswith(nm, "CON") ? :left : :center))
    end
end

# colour bar (vapor quality)
for k in 0:40
    xk = k / 40
    plot!(p, Shape([9.0 + 3.0 * xk, 9.0 + 3.0 * (k + 1) / 40, 9.0 + 3.0 * (k + 1) / 40, 9.0 + 3.0 * xk], [-1.6, -1.6, -1.25, -1.25]); fillcolor = xcol(xk), linecolor = xcol(xk), label = "")
end
annotate!(p, 9.0, -1.95, text("liquid (x = 0)", ip_font - 2, :black, :left))
annotate!(p, 12.2, -1.95, text("vapor (x = 1)", ip_font - 2, :black, :right))
annotate!(p, 10.5, -1.0, text("vapor quality x", ip_font - 1, :black))

# counts
nv(x) = count(v -> isfinite(v) && v > 0.5, x)
n2(x) = count(v -> isfinite(v) && 0.001 < v < 0.999, x)
info = "vapor cells (x > 0.5):  VL " * string(nv(x_vl)) * "/20,  CON " * string(nv(x_cond)) * "/50,  LL " * string(nv(x_ll)) * "/20" *
       "     two-phase cells (0.001 < x < 0.999):  VL " * string(n2(x_vl)) * ",  CON " * string(n2(x_cond)) * ",  LL " * string(n2(x_ll))
info1 = "vapor cells (x > 0.5):  VL " * string(nv(x_vl)) * "/20,  CON " * string(nv(x_cond)) * "/50,  LL " * string(nv(x_ll)) * "/20"
info2 = "two-phase cells (0.001 < x < 0.999):  VL " * string(n2(x_vl)) * ",  CON " * string(n2(x_cond)) * ",  LL " * string(n2(x_ll))
annotate!(p, -0.4, -1.2, text(info1, ip_font - 4, :black, :left))
annotate!(p, -0.4, -1.8, text(info2, ip_font - 4, :black, :left))
println(info)

savefig(p, joinpath(resdir, ip_out * ".png"))
if ip_vector
    savefig(p, joinpath(resdir, ip_out * ".pdf"))
    savefig(p, joinpath(resdir, ip_out * ".svg"))
end
display(p)
println("Saved: ", joinpath(resdir, ip_out * ".png"), ip_vector ? " (+ .pdf, .svg)" : "")