# =============================================================================
# test_residuals_wall.jl
#
# residuals_wall.jl の動作確認。流体は使わず、壁の熱収支の符号・保存性だけを
# 確認する (流体側 q と壁側 Q_fluid が同じ h_conv を共有することの確認を含む)。
#
# 実行前に、既存の physics.jl から k_sus304, c_sus304, ρ_sus304 を使える
# ようにしておくこと (下の include は、お手元の構成に合わせて調整)。
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays

const SRC = joinpath(@__DIR__, "..", "src")   # テストは ./test/、ソースは ./src/
include(joinpath(SRC, "physics.jl"))        # k_sus304, c_sus304, ρ_sus304 (既存)
include(joinpath(SRC, "smoothing.jl"))
include(joinpath(SRC, "types.jl"))
# heattransfer.jl の関数は引数の型注釈に FluidTables を使うため、
# proptables.jl / propaccess.jl が先に必要 (流体テーブル自体は構築しない)
include(joinpath(SRC, "proptables.jl"))
include(joinpath(SRC, "propaccess.jl"))
include(joinpath(SRC, "correlations.jl"))   # h_nc_holman
include(joinpath(SRC, "heattransfer.jl"))   # h_ambient_smooth
include(joinpath(SRC, "residuals_wall.jl"))

geo = TLGeometry(N=5, L=0.5, D_i=0.00435, D_o=0.00635)
wp = wall_params(geo)
N = geo.N

println("="^64)
println("テスト1: 熱平衡 (Tw=T流体=T周囲, 勾配なし) → 残差 ≈ 0")
println("="^64)
Tw = fill(300.0, N); Tf = fill(300.0, N); h = fill(1000.0, N)
du = zeros(N); res = zeros(N)
up = wall_ghost(300.0, wp); dn = wall_ghost(300.0, wp)
residual_wall_line!(res, du, Tw, Tf, h, wp, up, dn, 300.0, :pipe)
println("  max|res| = ", maximum(abs, res))

println()
println("="^64)
println("テスト2: 流体が壁より熱い → 壁が加熱される (resid<0 : dTw/dt>0 で釣り合う)")
println("="^64)
Tf = fill(320.0, N)
residual_wall_line!(res, du, Tw, Tf, h, wp, up, dn, 300.0, :pipe)
println("  resid の符号 (全て負なら壁が加熱される側): ", sign.(res))

println()
println("="^64)
println("テスト3: 凝縮器 (T_cond=285K < 壁300K) → 壁が冷やされる (resid>0)")
println("="^64)
Tf = fill(300.0, N)
residual_wall_line!(res, du, Tw, Tf, h, wp, up, dn, 285.0, :condenser)
println("  resid の符号 (全て正なら壁が冷却される側): ", sign.(res))

println()
println("="^64)
println("テスト4: 蒸気管 (T_amb=290K < 壁300K) → 放熱 (resid>0)。式30 更新後の符号確認")
println("="^64)
residual_wall_line!(res, du, Tw, Tf, h, wp, up, dn, 290.0, :pipe)
println("  resid の符号 (全て正なら放熱側): ", sign.(res))

println()
println("="^64)
println("テスト5: 軸方向熱伝導のコンダクタンスが素直な kA/dl に一致するか")
println("="^64)
G = wall_axial_conductance(300.0, wp, 300.0, wp)
println("  G = ", G, " W/K")
println("  k A / dl = ", k_sus304(300.0) * wp.A_cross / wp.dl, " W/K (一致するはず)")