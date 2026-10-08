# =============================================================================
# test_proptables.jl
#
# flash_ρe_once の動作確認用スクリプト。
# proptables.jl を include した状態で、単相液・単相蒸気・二相の3パターンを
# テストする。各テストでは、既知の (T,P) または (P,Q) から逆算した (ρ,e) を
# flash_ρe_once に渡し、元の (T,P,x) が復元されるかを確認する。
# =============================================================================

using CoolProp
include("../src/proptables.jl")   # with_abstract_state, flash_ρe_once を読み込む

fluid = "propylene"

# インデックスの事前取得 (build_ρe_table と同じ手順)
idx_DU = CoolProp.get_input_pair_index("DmassUmass_INPUTS")
idx_P  = CoolProp.get_param_index("P")
idx_T  = CoolProp.get_param_index("T")
idx_Q  = CoolProp.get_param_index("Q")

println("="^60)
println("テスト1: 過冷液 (T=300K, P=101325Pa を想定)")
println("="^60)
ρ1 = CoolProp.PropsSI("Dmass", "T", 268.0, "Q", 0, fluid)
e1 = CoolProp.PropsSI("Umass", "T", 268.0, "Q", 0, fluid)
ρ2 = CoolProp.PropsSI("Dmass", "T", 268.0, "Q", 1, fluid)
e2 = CoolProp.PropsSI("Umass", "T", 268.0, "Q", 1, fluid)
ρ3 = CoolProp.PropsSI("Dmass", "T", 340.0, "Q", 0, fluid)
e3 = CoolProp.PropsSI("Umass", "T", 340.0, "Q", 0, fluid)
ρ4 = CoolProp.PropsSI("Dmass", "T", 340.0, "Q", 1, fluid)
e4 = CoolProp.PropsSI("Umass", "T", 340.0, "Q", 1, fluid)
println("入力: ρ=$(ρ1) kg/m3, e=$(e1) J/kg")

println()
println("="^60)
println("テスト4 (任意): 小さな(ρ,e)格子を実際に構築して中身を見る")
println("="^60)
# 5×5程度の小さい格子ならすぐ終わるので、テーブル構築の挙動確認に使える
ρ_test_range = range(min(ρ1, ρ2, ρ3, ρ4) * 0.9, max(ρ1, ρ2, ρ3, ρ4) * 1.1, length=10)
e_test_range = range(min(e1, e2, e3, e4) * 0.9, max(e1, e2, e3, e4) * 1.1, length=10)

P_grid, T_grid, x_grid = build_ρe_table(fluid, ρ_test_range, e_test_range; verbose=true)

println("P_grid =")
display(P_grid)
println("\nT_grid =")
display(T_grid)
println("\nx_grid =")
display(x_grid)

T_sat_range = range(263.0, 343.0, 10)
@show build_saturation_tables(fluid, T_sat_range)