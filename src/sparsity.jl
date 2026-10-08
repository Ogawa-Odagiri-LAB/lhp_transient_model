# =============================================================================
# sparsity.jl (ヤコビアンの疎構造と、彩色した有限差分ヤコビアン)
#
# IDA は、各ステップで、ヤコビアン J = ∂F/∂u + cj ∂F/∂du を使って Newton 反復を
# 行う。密行列で、Sundials の内部の差分に任せると、n 個の変数について n 回の
# 残差評価が要る (n は数百になりうる)。J は疎 (各行が数個〜十数個の変数にしか
# 依存しない) なので、同時に摂動しても干渉しない列をまとめて (彩色)、評価回数を
# 色の数まで減らせる。疎行列の線形ソルバー (KLU) も使えるようになる。
#
# 【パターンの作り方: 数値的な検出】
# 当初は手動で依存関係を列挙する計画だったが、残差は共有する量
# (iface, Q_ec_grv, h_conv など) を経由して、遠く離れた変数同士が結合しており
# (蒸発器ケースの行がグルーブの全セルに依存する、など)、列挙の漏れが
# 静かに Newton 反復の収束を悪化させるので、数値的に検出する:
#   各変数 (u と du のそれぞれ) を1つずつ動かし、残差の各行が「ビット単位で」
#   変わるかを調べる。依存しない変数は、残差の計算が完全に同一なので、
#   浮動小数点でも厳密に同じ値になる。Dual数 (自動微分) は不要。
# 検出の見落としを避けるため:
#   1. 複数のランダムな状態の和集合を取る
#   2. 速度は、風上ブレンドの幅 ε_u 程度の小さな値にする (速度が大きいと
#      tanh が飽和し、上流側のセルへの依存が数値的に消えるが、流れが反転
#      したときには必要になる)
#   3. 対角成分は必ず含める (KLU の構造的特異性を避ける。代数変数の行は
#      ∂F/∂u の対角、微分変数の行は cj ∂F/∂du の対角が非ゼロ)
# 検証は test_sparsity.jl で、別の状態の密な有限差分ヤコビアンと比較する。
# =============================================================================

using Random
using SparseArrays: SparseMatrixCSC, sparse, nonzeros, rowvals, nzrange, findnz
using ComponentArrays: ComponentArray, getdata, getaxes

# -----------------------------------------------------------------------------
# 1. 変数の大きさ (摂動幅の下限に使う)
# -----------------------------------------------------------------------------

"""
    variable_scales(geo, condition) -> Vector{Float64}

状態ベクトルと同じ並びの、各変数の代表的な大きさ (types.jl の *_scale)。
有限差分の摂動幅の下限と、IDA の絶対許容誤差に使う。
"""
function variable_scales(geo::LHPGeometry, condition::LHPTransientCondition)
    ones_state = initialize_state(geo)
    getdata(ones_state) .= 1.0
    return copy(getdata(restore_state(ones_state, condition)))
end

# -----------------------------------------------------------------------------
# 2. パターンの検出
# -----------------------------------------------------------------------------

"""
    random_probe_state(u0, scales, condition, rng) -> (u, du)

パターン検出用の状態。IC の近傍で、各変数を相対 1e-3 だけランダムに動かし、
速度は風上ブレンドの幅 ε_u 程度の小さな値 (±2ε_u) にする。du は、微分変数に
ランダムな値を入れる (F は du について線形なので、値の大きさは関係ない)。
物理的に妥当な範囲 (物性テーブルの範囲内) に留めるため、摂動は小さくする。
"""
function random_probe_state(u0::AbstractVector, p::LHPParams, scales, diffmask, rng)
    n = length(u0)
    u = u0 .* (1 .+ 1.0e-3 .* randn(rng, n))
    # 速度: 風上ブレンドが飽和しない大きさ。ComponentArray で u の位置を取り出す
    cu = ComponentArray(u, p.axes)
    ε = p.cond.ε_u
    for comp in (:grv, :vl, :cond, :ll, :cc)
        v = getproperty(getproperty(cu, comp), :u)
        v .= (2 .* rand(rng, length(v)) .- 1) .* 2ε
    end
    du = zeros(n)
    du[diffmask] .= randn(rng, count(diffmask))
    return u, du
end

"""
    detect_jacobian_pattern(p, u0, t0, scales, diffmask; nprobe=3, rel=1e-4, seed=0)
        -> BitMatrix  (n×n, 行 = 残差、列 = 変数)

残差 F(t, u, du) のヤコビアン (∂F/∂u と ∂F/∂du の和集合) の非ゼロ構造を、
各変数を1つずつ動かして、残差がビット単位で変わる行を調べて検出する。
nprobe 個のランダムな状態の和集合。対角は必ず含める。

du の摂動は、F が du について線形なので、大きさ1で十分。u の摂動は、
rel × max(|u_j|, 1e-3 scale_j)。
"""
function detect_jacobian_pattern(p::LHPParams, u0::AbstractVector, t0, scales, diffmask;
    nprobe::Int=3, rel::Float64=1.0e-4, seed::Int=0)

    n = length(u0)
    rng = Random.MersenneTwister(seed)
    pattern = falses(n, n)
    F0 = zeros(n)
    F1 = zeros(n)
    for _ in 1:nprobe
        u, du = random_probe_state(u0, p, scales, diffmask, rng)
        f_dae_residual!(F0, du, u, p, t0)
        for j in 1:n
            # u を動かす
            δ = rel * max(abs(u[j]), 1.0e-3 * scales[j])
            uj = u[j]
            u[j] = uj + δ
            f_dae_residual!(F1, du, u, p, t0)
            u[j] = uj
            for i in 1:n
                F1[i] != F0[i] && (pattern[i, j] = true)
            end
            # du を動かす
            duj = du[j]
            du[j] = duj + 1.0
            f_dae_residual!(F1, du, u, p, t0)
            du[j] = duj
            for i in 1:n
                F1[i] != F0[i] && (pattern[i, j] = true)
            end
        end
    end
    for i in 1:n
        pattern[i, i] = true
    end
    return pattern
end

# -----------------------------------------------------------------------------
# 3. 彩色
# -----------------------------------------------------------------------------

"""
    greedy_coloring(pattern) -> Vector{Int}   (各列の色 1..ncolors)

列の貪欲な彩色。同じ行に非ゼロを持つ列は、同時に摂動すると干渉するので、
同じ色にできない。非ゼロの多い列から順に、使えない色を避けて最小の色を割り当てる。
"""
function greedy_coloring(pattern::AbstractMatrix{Bool})
    nr, nc = size(pattern)
    rows_of = [findall(@view pattern[:, j]) for j in 1:nc]
    cols_of = [findall(@view pattern[i, :]) for i in 1:nr]
    color_order = sortperm(length.(rows_of), rev=true)
    color = zeros(Int, nc)
    for j in color_order
        forbidden = Set{Int}()
        for i in rows_of[j]
            for k in cols_of[i]
                color[k] > 0 && push!(forbidden, color[k])
            end
        end
        c = 1
        while c in forbidden
            c += 1
        end
        color[j] = c
    end
    return color
end

"""
    check_coloring(pattern, color) -> Bool

彩色の妥当性: 同じ色の列が、同じ行に非ゼロを持たないこと。
"""
function check_coloring(pattern::AbstractMatrix{Bool}, color::Vector{Int})
    for c in 1:maximum(color)
        cols = findall(==(c), color)
        for i in 1:size(pattern, 1)
            count(j -> pattern[i, j], cols) > 1 && return false
        end
    end
    return true
end

# -----------------------------------------------------------------------------
# 4. 彩色した有限差分ヤコビアン (IDA の jac として使う)
# -----------------------------------------------------------------------------

"""
    ColoredFDJacobian(p, pattern, color, scales)

IDA の DAEFunction の jac に渡す、彩色した有限差分のヤコビアン。
  J = ∂F/∂u + γ ∂F/∂du       (γ = IDA の cj)
色ごとに、その色の全列を同時に動かす (u に δ、du に γδ)。F は du について
線形なので、1回の評価で、その色の各列の J の列 (u と du の寄与の和) が、
δ の1次の誤差で得られる (IDA 自身の差分と同じ方法)。
J は prototype (パターンの疎行列) と同じ構造の SparseMatrixCSC。
"""
struct ColoredFDJacobian{P}
    p::P
    prototype::SparseMatrixCSC{Float64,Int}
    color_cols::Vector{Vector{Int}}
    scales::Vector{Float64}
    F0::Vector{Float64}
    F1::Vector{Float64}
    u1::Vector{Float64}
    du1::Vector{Float64}
    δ::Vector{Float64}
end

function ColoredFDJacobian(p, pattern::AbstractMatrix{Bool}, color::Vector{Int}, scales)
    n = size(pattern, 1)
    I, J = findnz(sparse(pattern))
    prototype = sparse(I, J, ones(length(I)), n, n)
    color_cols = [findall(==(c), color) for c in 1:maximum(color)]
    return ColoredFDJacobian(p, prototype, color_cols, collect(Float64, scales),
                             zeros(n), zeros(n), zeros(n), zeros(n), zeros(n))
end

function (jf::ColoredFDJacobian)(J, du, u, p, γ, t)
    f_dae_residual!(jf.F0, du, u, jf.p, t)
    nz = nonzeros(J)
    rv = rowvals(J)
    fill!(nz, 0.0)
    for cols in jf.color_cols
        copyto!(jf.u1, u)
        copyto!(jf.du1, du)
        for j in cols
            d = sqrt(eps(Float64)) * max(abs(u[j]), 1.0e-3 * jf.scales[j])
            jf.δ[j] = d
            jf.u1[j] += d
            jf.du1[j] += γ * d
        end
        f_dae_residual!(jf.F1, jf.du1, jf.u1, jf.p, t)
        for j in cols
            for k in nzrange(J, j)
                i = rv[k]
                nz[k] = (jf.F1[i] - jf.F0[i]) / jf.δ[j]
            end
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------
# 5. パターンの要約 (人が読んで確認する用)
# -----------------------------------------------------------------------------

"""
    variable_labels(axes, n) -> Vector{String}

状態ベクトルの各要素の名前 ("comp.var[k]")。
"""
function variable_labels(axes, n::Int)
    idx = ComponentArray(collect(1:n), axes)
    labels = fill("", n)
    for comp in keys(idx)
        sub = getproperty(idx, comp)
        for var in keys(sub)
            v = getproperty(sub, var)
            for (k, i) in enumerate(v)
                labels[i] = "$(comp).$(var)[$(k)]"
            end
        end
    end
    return labels
end

"""
    pattern_summary(pattern, axes)

ブロック (comp.var) ごとに、各行が依存するブロックと、その非ゼロ数を表示する。
手動で期待した結合と、検出された結合を、人が見比べるための要約。
"""
function pattern_summary(pattern::AbstractMatrix{Bool}, axes)
    n = size(pattern, 1)
    labels = variable_labels(axes, n)
    block(l) = replace(l, r"\[\d+\]$" => "")
    blocks = unique(block.(labels))
    bidx = Dict(b => findall(l -> block(l) == b, labels) for b in blocks)
    println("パターン: $(count(pattern)) / $(n * n) 非ゼロ ($(round(100 * count(pattern) / n^2, digits=2)) %)")
    println("行ごとの非ゼロ数の最大: ", maximum(sum(pattern, dims=2)),
            ",  列ごとの非ゼロ数の最大: ", maximum(sum(pattern, dims=1)))
    println("各ブロックの行が、どのブロックの変数に依存するか (括弧内は、その行のブロック内での、依存する変数の数の最大):")
    for b in blocks
        rows = bidx[b]
        deps = Dict{String,Int}()
        for i in rows, j in 1:n
            if pattern[i, j]
                bj = block(labels[j])
                deps[bj] = max(get(deps, bj, 0), count(jj -> pattern[i, jj] && block(labels[jj]) == bj, 1:n))
            end
        end
        println("  ", rpad(b, 12), " ← ", join(["$(k)($(deps[k]))" for k in blocks if haskey(deps, k)], ", "))
    end
end