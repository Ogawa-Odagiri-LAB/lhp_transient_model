# =============================================================================
# diagnose_grid.jl
#
# (ρ,e)格子の各点でflashが成功するかどうかを可視化する。
# fill_invalid_cells!を呼ぶ前の、生の成否パターンを見ることで、
# 失敗が「四隅に偏っている(矩形はみ出し)」のか「ランダムに散らばっている
# (別の問題)」のかを切り分ける。
# =============================================================================

using CoolProp
include("../src/proptables.jl")

fluid = "Water"   # ご自身の作動流体に置き換えてください

idx_DU = CoolProp.get_input_pair_index("DmassUmass_INPUTS")
idx_P  = CoolProp.get_param_index("P")
idx_T  = CoolProp.get_param_index("T")
idx_Q  = CoolProp.get_param_index("Q")

# 運転範囲は、お使いの init_fluid_lookups に渡しているものと合わせてください
T_min, T_max = 280.0, 400.0
P_min, P_max = 1.0e5, 2.0e5

ρ_lo, ρ_hi, e_lo, e_hi = estimate_ρe_bounds(fluid, T_min, T_max, P_min, P_max)
println("見積もられた範囲: ρ∈[$ρ_lo, $ρ_hi], e∈[$e_lo, $e_hi]")

Nρ, Ne = 25, 25
ρ_range = range(ρ_lo, ρ_hi, length=Nρ)
e_range = range(e_lo, e_hi, length=Ne)

valid = falses(Nρ, Ne)
with_abstract_state("HEOS", fluid) do handle
    for (i, ρ) in enumerate(ρ_range), (j, e) in enumerate(e_range)
        try
            flash_ρe_once(handle, ρ, e, fluid, idx_DU, idx_P, idx_T, idx_Q)
            valid[i, j] = true
        catch
            valid[i, j] = false
        end
    end
end

n_valid = count(valid)
println("成功率: $(n_valid) / $(Nρ*Ne) = $(round(100*n_valid/(Nρ*Ne), digits=1))%")
println()
println("成否マップ (■=成功, ・=失敗)")
println("  縦軸: ρ (上に行くほど大きい)   横軸: e (右に行くほど大きい)")
println()
for i in Nρ:-1:1
    println(join([valid[i, j] ? "■" : "・" for j in 1:Ne]))
end