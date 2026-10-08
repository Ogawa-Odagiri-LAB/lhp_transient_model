# =============================================================================
# test_proptables.jl
#
# flash_ρe_once の動作確認用スクリプト。
# proptables.jl を include した状態で、単相液・単相蒸気・二相の3パターンを
# テストする。各テストでは、既知の (T,P) または (P,Q) から逆算した (ρ,e) を
# flash_ρe_once に渡し、元の (T,P,x) が復元されるかを確認する。
# =============================================================================

using CoolProp
const SRC = joinpath(@__DIR__, "..", "src")   # テストは ./test/、ソースは ./src/
include(joinpath(SRC, "proptables.jl"))   # with_abstract_state, flash_ρe_once を読み込む

fluid = "Water"

# インデックスの事前取得 (build_ρe_table と同じ手順)
idx_DU = CoolProp.get_input_pair_index("DmassUmass_INPUTS")
idx_P  = CoolProp.get_param_index("P")
idx_T  = CoolProp.get_param_index("T")
idx_Q  = CoolProp.get_param_index("Q")
idx_V  = CoolProp.get_param_index("V")
idx_L  = CoolProp.get_param_index("L")
idx_C  = CoolProp.get_param_index("C")

println("="^60)
println("テスト1: 過冷液 (T=300K, P=101325Pa を想定)")
println("="^60)
ρ1 = CoolProp.PropsSI("Dmass", "T", 300.0, "P", 101325.0, fluid)
e1 = CoolProp.PropsSI("Umass", "T", 300.0, "P", 101325.0, fluid)
println("入力: ρ=$(ρ1) kg/m3, e=$(e1) J/kg")

with_abstract_state("HEOS", fluid) do handle
    P, T, x, μ, k, cp = flash_ρe_once(handle, ρ1, e1, fluid, idx_DU, idx_P, idx_T, idx_Q, idx_V, idx_L, idx_C)
    println("出力: P=$(P) Pa, T=$(T) K, x=$(x)")
    println("期待値: P≈101325, T≈300, x=0.0 (過冷液)")
end

println()
println("="^60)
println("テスト2: 過熱蒸気 (T=500K, P=101325Pa を想定)")
println("="^60)
ρ2 = CoolProp.PropsSI("Dmass", "T", 500.0, "P", 101325.0, fluid)
e2 = CoolProp.PropsSI("Umass", "T", 500.0, "P", 101325.0, fluid)
println("入力: ρ=$(ρ2) kg/m3, e=$(e2) J/kg")

with_abstract_state("HEOS", fluid) do handle
    P, T, x, μ, k, cp = flash_ρe_once(handle, ρ2, e2, fluid, idx_DU, idx_P, idx_T, idx_Q, idx_V, idx_L, idx_C)
    println("出力: P=$(P) Pa, T=$(T) K, x=$(x)")
    println("期待値: P≈101325, T≈500, x=1.0 (過熱蒸気)")
end

println()
println("="^60)
println("テスト3: 二相 (P=101325Pa, Q=0.5 を想定)")
println("="^60)
ρ3 = CoolProp.PropsSI("Dmass", "P", 101325.0, "Q", 0.5, fluid)
e3 = CoolProp.PropsSI("Umass", "P", 101325.0, "Q", 0.5, fluid)
println("入力: ρ=$(ρ3) kg/m3, e=$(e3) J/kg")

with_abstract_state("HEOS", fluid) do handle
    P, T, x, μ, k, cp = flash_ρe_once(handle, ρ3, e3, fluid, idx_DU, idx_P, idx_T, idx_Q, idx_V, idx_L, idx_C)
    println("出力: P=$(P) Pa, T=$(T) K, x=$(x)")
    println("期待値: P≈101325, T≈373.1 (1atmでの沸点), x≈0.5")
end

println()
println("="^60)
println("テスト4 (任意): 小さな(ρ,e)格子を実際に構築して中身を見る")
println("="^60)
# 5×5程度の小さい格子ならすぐ終わるので、テーブル構築の挙動確認に使える
ρ_test_range = range(ρ1 * 0.9, ρ2 * 1.1, length=5)
e_test_range = range(min(e1, e2, e3) * 0.9, max(e1, e2, e3) * 1.1, length=5)

P_grid, T_grid, x_grid = build_ρe_table(fluid, ρ_test_range, e_test_range; verbose=true)

println("P_grid =")
display(P_grid)
println("\nT_grid =")
display(T_grid)
println("\nx_grid =")
display(x_grid)