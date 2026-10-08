# =============================================================================
# test_geometry_cc.jl
#
# smoothing.jl (soft_clamp の修正版) と geometry_cc.jl の動作確認。
# =============================================================================

const SRC = joinpath(@__DIR__, "..", "src")   # テストは ./test/、ソースは ./src/
include(joinpath(SRC, "smoothing.jl"))
include(joinpath(SRC, "geometry_cc.jl"))

println("="^60)
println("テスト1: soft_clamp (修正版) — 中間域が恒等的に通るか")
println("="^60)
for x in (-0.5, 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.5)
    println("  soft_clamp($x, 0, 1) = ", soft_clamp(x, 0.0, 1.0))
end
println("期待: 0.1〜0.9 はほぼ入力と同じ値、-0.5→≈0、1.5→≈1")

println()
println("="^60)
println("テスト2: f_l → θ_liq → f_l の往復整合 (テーブルの精度)")
println("="^60)
maxerr = 0.0
for f in range(1e-3, 1 - 1e-3, length=1000)
    θ_liq = theta_of_liquid_fraction(f)
    maxerr = max(maxerr, abs(f_liquid_of_theta(θ_liq) - f))
end
println("  往復の最大誤差 |f_l(θ_liq(f))-f| = ", maxerr, " (1e-4 以下であれば良好)")

println()
println("="^60)
println("テスト3: cc_liquid_geometry (直径35mm, セル長70mm の例)")
println("="^60)
R_cc = 0.035 / 2   # (physics.jl の気体定数 R との衝突を避ける)
dl = 0.07
for f in (-0.2, 0.0, 0.1, 0.5, 0.9, 1.0, 1.05)
    g = cc_liquid_geometry(f, R_cc, dl)
    println("  f_l=$f → f_eff=$(round(g.f_eff, digits=5)), θ_liq=$(round(g.θ_liq, digits=4)) rad, ",
            "h_liq=$(round(g.h_liq*1e3, digits=3)) mm, ",
            "A_liq+A_vap=$(g.A_liq + g.A_vap) (期待: $(2π * R_cc * dl))")
end
g = cc_liquid_geometry(0.5, R_cc, dl)
println("  f_l=0.5 で h_liq/R_cc = ", g.h_liq / R_cc, " (期待: 1.0)")