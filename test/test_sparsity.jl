# =============================================================================
# test_sparsity.jl
#
# sparsity.jl の動作確認。
#   1. パターンの数値的な検出 (非ゼロの割合、色の数)
#   2. 彩色の妥当性 (同じ色の列が、同じ行に非ゼロを持たない)
#   3. パターンの網羅性: 別の状態で作った密な有限差分ヤコビアンの有意な非ゼロ要素が、
#      すべてパターンに含まれること (含まれない要素があれば、検出の見落とし)
#   4. 彩色した有限差分ヤコビアンが、密な有限差分ヤコビアンと一致すること
#   5. ブロックごとの結合の要約 (人が見て、物理的な結合と合っているか確認する)
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays,
      SparseArrays, LinearAlgebra, Random

const SRC = joinpath(@__DIR__, "..", "src")   # テストは ./test/、ソースは ./src/
for f in ("physics", "smoothing", "types", "proptables", "propaccess", "correlations",
          "heattransfer", "geometry_cc", "residuals_transport", "residuals_wall",
          "residuals_wick", "residuals_cc", "residuals_evapcase", "residuals_grv",
          "residuals_main", "initialcondition", "sparsity")
    include(joinpath(SRC, f * ".jl"))
end

fluid = "ammonia"
println("テーブル構築中 (初回は数分かかる場合があります。キャッシュされます)...")
tables = init_fluid_lookups(fluid; T_min=230.0, T_max=360.0, P_min=1.0e5, P_max=4.5e6,
                            Ne=400, Δδ=1.0, N_T_sat=200)

# 小さめの分割 (検証の評価回数を抑える)
wick = WickGeometry(N=3, φ=0.65, ρ_wick=8900.0, Cp_wick=444.0)
geo = LHPGeometry(wick,
                  GrvGeometry(N=1),
                  TLGeometry(N=4, L=1.0, D_i=0.00435, D_o=0.00635),
                  CondGeometry(N=4, L=1.0, D_i=0.00435, D_o=0.00635),
                  TLGeometry(N=4, L=1.0, D_i=0.00218, D_o=0.00318),
                  CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037))
condition = LHPTransientCondition(Q_apply=20.0, t_init=0.0, t_end=100.0, fluid=fluid)

u0, du0, p = initial_condition(geo, condition, tables; verbose=false)
u0vec = copy(getdata(u0))
n = length(u0vec)
diffmask = Bool.(getdata(initialize_differential_vars(geo)))
scales = variable_scales(geo, condition)
t0 = condition.t_init
println("変数の数 n = ", n, " (微分変数 ", count(diffmask), ", 代数変数 ", n - count(diffmask), ")")

println("="^72)
println("テスト1: パターンの検出")
println("="^72)
t_detect = @elapsed pattern = detect_jacobian_pattern(p, u0vec, t0, scales, diffmask)
println("  検出時間: ", round(t_detect, digits=2), " s")
println("  非ゼロ: ", count(pattern), " / ", n^2, " (", round(100 * count(pattern) / n^2, digits=2), " %)")
println("  行ごとの非ゼロ数の最大: ", maximum(sum(pattern, dims=2)),
        ",  列ごとの最大: ", maximum(sum(pattern, dims=1)))
println("  対角がすべて含まれる: ", all(pattern[i, i] for i in 1:n))

println()
println("="^72)
println("テスト2: 彩色")
println("="^72)
color = greedy_coloring(pattern)
println("  色の数: ", maximum(color), " (密行列の差分なら ", n, " 回の残差評価 → ", maximum(color), " 回)")
println("  彩色の妥当性 (同色の列が行を共有しない): ", check_coloring(pattern, color), " (true のはず)")

# --- 密な有限差分ヤコビアン (検証用) ---
function dense_fd_jacobian(u, du, γ)
    F0 = zeros(n); F1 = zeros(n); J = zeros(n, n)
    f_dae_residual!(F0, du, u, p, t0)
    for j in 1:n
        d = 1.0e-7 * max(abs(u[j]), 1.0e-3 * scales[j])
        u1 = copy(u); du1 = copy(du)
        u1[j] += d; du1[j] += γ * d
        f_dae_residual!(F1, du1, u1, p, t0)
        J[:, j] = (F1 .- F0) ./ d
    end
    return J
end

println()
println("="^72)
println("テスト3: パターンの網羅性 (別の状態で、密な差分ヤコビアンの有意な要素がパターン内にあるか)")
println("="^72)
rng = MersenneTwister(123)
states = []
for trial in 1:3
    u, du = random_probe_state(u0vec, p, scales, diffmask, MersenneTwister(1000 + trial))
    du .*= 1.0e-3
    push!(states, (u, du))
end
for (k, (u, du)) in enumerate(states)
    J = dense_fd_jacobian(u, du, 10.0)
    rowmax = [maximum(abs, J[i, :]) + 1.0e-300 for i in 1:n]
    miss = [(i, j) for i in 1:n, j in 1:n if abs(J[i, j]) > 1.0e-7 * rowmax[i] && !pattern[i, j]]
    println("  状態$(k): パターン外の有意な要素 = ", length(miss), " 個")
    if !isempty(miss)
        labels = variable_labels(p.axes, n)
        for (i, j) in first(miss, 10)
            println("     行 ", labels[i], " ← 列 ", labels[j], "  |J|=", abs(J[i, j]))
        end
    end
end
println("  (0個なら、検出は、この状態で見落としなし。1個以上あれば、nprobe を増やすか、")
println("   random_probe_state の摂動を見直す)")

println()
println("="^72)
println("テスト4: 彩色した差分ヤコビアン = 密な差分ヤコビアン")
println("="^72)
jf = ColoredFDJacobian(p, pattern, color, scales)
J_sparse = copy(jf.prototype)
for (k, (u, du)) in enumerate(states)
    γ = 10.0 * k
    Jd = dense_fd_jacobian(u, du, γ)
    jf(J_sparse, du, u, p, γ, t0)
    Jc = Matrix(J_sparse)
    err = abs.(Jd .- Jc)
    rowmax = [maximum(abs, Jd[i, :]) + 1.0e-300 for i in 1:n]
    rel = maximum(err[i, j] / rowmax[i] for i in 1:n, j in 1:n)
    println("  状態$(k) (γ=$(γ)): 差の最大 = ", maximum(err), ",  行最大に対する相対 = ", rel)
end
println("  (相対 1e-6 程度以下なら良好。パターンが覆っていれば、ほぼ一致する)")

println()
println("="^72)
println("テスト5: ブロックごとの結合の要約")
println("="^72)
pattern_summary(pattern, p.axes)