# =============================================================================
# test_types.jl
#
# types.jl の initialize_state / scale_state / restore_state の動作確認。
# T_wall.evap (蒸発器ケース) は長さ1の配列で、他の壁と同じ ./= でスケーリングされる。
# =============================================================================

using ComponentArrays, Parameters, StaticArrays

const SRC = joinpath(@__DIR__, "..", "src")   # テストは ./test/、ソースは ./src/
include(joinpath(SRC, "types.jl"))

wick = WickGeometry(N=3, φ=0.65, ρ_wick=8900.0, Cp_wick=444.0)
geo = LHPGeometry(wick,
                  GrvGeometry(N=1),
                  TLGeometry(N=5, L=2.0, D_i=0.00435, D_o=0.00635),
                  CondGeometry(N=5, L=1.0, D_i=0.00435, D_o=0.00635),
                  TLGeometry(N=5, L=2.0, D_i=0.00218, D_o=0.00318),
                  CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037))
cond = LHPTransientCondition(Q_apply=20.0, t_init=0.0, t_end=1.0)

u = initialize_state(geo)
u.vl.ρm .= 20.0;  u.vl.u .= 1.0;  u.vl.Em .= 3.0e5;  u.vl.P .= 1.0e6;  u.vl.T .= 300.0;  u.vl.x .= 0.5
u.wick.E .= 100.0;  u.wick.T .= 295.0
u.T_wall.evap .= 310.0
u.T_wall.vl .= 305.0;  u.T_wall.cond .= 290.0;  u.T_wall.ll .= 295.0;  u.T_wall.cc .= 298.0

println("="^64)
println("テスト1: scale_state で T_wall.evap (長さ1) と T_wall.vl の両方が更新される")
println("="^64)
s = scale_state(u, cond)
println("  T_wall.evap[1] : ", u.T_wall.evap[1], " → ", s.T_wall.evap[1],
        "  (期待: ", 310.0 / cond.T_scale, ")")
println("  T_wall.vl[1]: ", u.T_wall.vl[1], " → ", s.T_wall.vl[1],
        "  (期待: ", 305.0 / cond.T_scale, ")")
println("  元の u は変更されていないか: u.T_wall.evap[1] = ", u.T_wall.evap[1], " (310.0 のはず)")

println()
println("="^64)
println("テスト2: restore_state(scale_state(u)) が元に戻る")
println("="^64)
r = restore_state(s, cond)
maxerr = maximum(abs, getdata(r) .- getdata(u))
println("  最大誤差 = ", maxerr, " (丸め誤差程度のはず)")

println()
println("="^64)
println("テスト3: differential_vars のマスク (CCの u は代数変数、T_wall は全て微分変数)")
println("="^64)
mask = initialize_differential_vars(geo)
println("  cc.u (代数 = false のはず): ", mask.cc.u)
println("  cc.ρm (微分 = true のはず): ", mask.cc.ρm)
println("  vl.u (微分 = true のはず): ", mask.vl.u)
println("  T_wall.evap (微分 = true のはず): ", mask.T_wall.evap[1])
println("  状態ベクトルとマスクの長さが一致: ", length(getdata(u)) == length(getdata(mask)))