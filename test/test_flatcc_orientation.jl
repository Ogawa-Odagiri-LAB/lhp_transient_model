# =============================================================================
# test_flatcc_orientation.jl
#
# 平板型CCの姿勢 (FlatCCGeometry.orientation) の切り替えの確認。
#   テスト1: 上面加熱の液の円弓形の面積分率が 1 − f_v に一致するか (幾何の整合)
#   テスト2: 上面加熱の接触面積の合計 = 曲面 + 両端面 + 平面 (S_total) か
#   テスト3: 下面加熱の自然対流が、改修前の式と完全に一致するか (回帰)
#   テスト4: 同じ状態での下面加熱 / 上面加熱の熱 Q_fluid の比較 (目安)
# 実行: プロジェクトのルートで julia --project=. test/test_flatcc_orientation.jl
# =============================================================================

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using LHPTransientModel
using ComponentArrays: ComponentArray
const M = LHPTransientModel

tables = init_fluid_lookups("ammonia", T_min = 230.0, T_max = 360.0, P_max = 4.5e6)
g = 9.8

cc_bot = M.FlatCCGeometry(N = 1)
cc_top = M.FlatCCGeometry(N = 1, orientation = :heat_top)
p_bot = M.flatcc_params(cc_bot)
p_top = M.flatcc_params(cc_top)
R, dl = p_top.R, p_top.dl

println("="^60)
println("テスト1: 上面加熱の円弓形の面積分率 = 1 − f_v")
println("="^60)
maxerr = 0.0
for f_v in range(0.01, 0.99, length = 200)
    geom = M.cc_bottom_liquid_geometry(f_v, R, dl)
    θ = geom.θ_liq
    A_seg = 0.5 * R^2 * (θ - sin(θ))            # 円弓形の面積
    f_l = A_seg / (0.5 * π * R^2)               # 半円に占める分率
    global maxerr = max(maxerr, abs(f_l - (1 - geom.f_eff)))
    # 液深さ: 円弓形の弦の高さ R(1 − cos(θ/2))
end
println("  最大誤差 |f_l − (1 − f_v)| = ", maxerr, "  (1e-4 以下であれば良好)")

println()
println("="^60)
println("テスト2: 上面加熱の接触面積の合計")
println("="^60)
for f_v in (0.05, 0.3, 0.5, 0.8, 0.95)
    geom = M.cc_bottom_liquid_geometry(f_v, R, dl)
    θ = geom.θ_liq
    A_liq = R * θ * dl + R^2 * (θ - sin(θ))
    A_vap = (p_top.S_inner - A_liq) + p_top.A_flat
    println("  f_v=$(f_v): h_liq=$(round(geom.h_liq*1e3, digits=2)) mm, h_vap=$(round(geom.h_vap*1e3, digits=2)) mm, ",
            "A_liq+A_vap−S_total = ", A_liq + A_vap - p_top.S_total)
end

# CC の状態 (飽和、乾き度 x から ρm を作る)
function cc_state(T, x)
    sp = M.sat_props(T, tables)
    v = x / sp.ρ_v + (1 - x) / sp.ρ_l
    return ComponentArray(ρm = [1 / v], u = [0.0], Em = [0.0], P = [M.Psat(T, tables)], T = [T], x = [x])
end

println()
println("="^60)
println("テスト3: 下面加熱が改修前の式と一致するか (回帰)")
println("="^60)
let T = 300.0, T_sw = 301.0, x = 0.01
    u = cc_state(T, x)
    nc = M.flatcc_natural_convection(u, T_sw, p_bot, g, tables)
    # 改修前の式 (residuals_flatcc.jl の旧版をそのまま書き下したもの)
    sp = M.sat_props(T, tables)
    f_v = M.vapor_volume_fraction(u.ρm[1], u.x[1], sp.ρ_v)
    geom = M.cc_vapor_geometry(f_v, p_bot.R, p_bot.dl)
    θ = geom.θ_vap
    A_vap = p_bot.R * θ * p_bot.dl + p_bot.R^2 * (θ - sin(θ))
    A_liq = p_bot.S_total - A_vap
    Q_old = nc.h_l * A_liq * (T_sw - T) + nc.h_v * A_vap * (T_sw - T)
    println("  Q_fluid (改修後) = ", nc.Q_fluid, " W, 改修前の式 = ", Q_old, " W, 一致: ", nc.Q_fluid === Q_old)
end

println()
println("="^60)
println("テスト4: 同じ状態での Q_fluid の比較 (T_swick − T_cc = 1 K)")
println("="^60)
for x in (0.002, 0.01, 0.05)
    u = cc_state(300.0, x)
    nb = M.flatcc_natural_convection(u, 301.0, p_bot, g, tables)
    nt = M.flatcc_natural_convection(u, 301.0, p_top, g, tables)
    println("  x=$(x): f_v=$(round(nt.geom.f_eff, digits=3)), 下面加熱 Q=$(round(nb.Q_fluid, digits=3)) W, ",
            "上面加熱 Q=$(round(nt.Q_fluid, digits=3)) W (Q_l=$(round(nt.Q_l, digits=3)), Q_v=$(round(nt.Q_v, digits=3)), ",
            "平面側 Ra=$(round(nt.Ra_flat, sigdigits=3)), h=$(round(nt.h_flat, digits=1)) W/m²K)")
end