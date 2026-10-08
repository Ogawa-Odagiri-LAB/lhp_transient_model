# =============================================================================
# analyze_oscillation.jl (液管の速度振動の包絡線の解析。保存済みの CSV だけを使う)
#
# 目的: 液管の速度の振動が、時間とともに、減衰するか、減衰した後に何が残るか、を調べる。
#
# 診断 (test_diagnose_oscillation.jl) の結果:
#   ・ステップ後 (t=60〜64 s) の振動は、rtol を 100 倍変えても、σ_u ≈ 0.027 m/s、周波数 0.75〜1.25 Hz
#     で、変わらない → 許容誤差の雑音ではない
#   ・液管の両端と平均の速度が、1% 以内で一致 → 液柱が、剛体として往復
# このスクリプトは、保存済みの結果 (saveat=1 s) から、
#   1. 窓ごとの σ(u_ll) と、減衰する局面の e-fold 時間 (floor の 2 倍を超える窓だけを近似する)
#   2. 減衰した後に残る、持続するゆらぎ (floor) の大きさ (平均流速に対する比)
#   3. 圧力の突発的な増加の時刻と、凝縮器の前線のセルの状態の変化 (液セル数、乾きセル数) の時刻の照合
# を行う。
#
# 【経緯】 最初の版は、減衰の直線近似に、floor の窓まで含めてしまい (R²=0.27)、意味のない
# e-fold 時間 (170 s) を出した。また、4 s の窓の前半と後半の σ の比較は、0.75 Hz と 1.25 Hz の
# 2 つのモードの、うなり (周期 約 2 s) の影響を受け、減衰の判定には、使えなかった。
#
# 注意: saveat=1 s では、1 Hz 付近の振動は、折り返して見える (Nyquist 周波数 0.5 Hz)。
# 窓の標準偏差は、振動の振幅の目安としては使える (窓の中で位相がばらつくため) が、
# 周波数は分からない。周波数は、test_diagnose_oscillation.jl で調べる。
#
# 使い方:  julia --project analyze_oscillation.jl   (または REPL で include)
# =============================================================================

using CSV, DataFrames, Statistics

# ---------------------------- 設定 ----------------------------
RES_DIR = joinpath(@__DIR__, "..", "res", "transient_ammonia")
T_STEP = 200.0          # Q のステップアップの開始時刻 [s] (transient.jl の t_step)
TAU_STEP = 1.0          # 立ち上がり時間 [s]
W = 10.0                # 窓の幅 [s]
FRICTION_TAU = 1.3      # 層流摩擦だけの場合の、振幅の e-fold 時間 [s] (test_diagnose_oscillation の見積もり)
T_QUIET = 40.0          # 起動の過渡の後、準定常とみなす時刻 [s] (floor の推定に使う)
DECAY_FLOOR_FACTOR = 1.5   # 減衰の近似に使う窓: σ が、floor のこの倍数を超えている窓 (大きいほど、floor の揺らぎを避けるが、窓が減る)
SPIKE_FACTOR = 6.0      # 圧力の突発的な増加とみなす閾値: 中央値からのずれが、(MAD に基づく標準偏差) の何倍か
SPIKE_MIN = 80.0        # 突発的な増加とみなす、最小の圧力差 [Pa]
MATCH_WINDOW = 5.0      # 事象の照合: 前線のセルの状態の変化と、圧力の増加が、何秒以内なら、対応とみなすか [s]
X_LIQ = 1.0e-3          # 凝縮器のセルが「液」とみなす乾き度の上限
X_DRY = 0.999           # 凝縮器のセルが「乾き蒸気」とみなす乾き度の下限
# --------------------------------------------------------------

prof = CSV.read(joinpath(RES_DIR, "profiles.csv"), DataFrame)

# 成分 comp の、量 var の、セル平均の時系列
function series(comp::AbstractString, var::Symbol)
    sub = prof[prof.comp .== comp, :]
    g = combine(groupby(sub, :t), var => mean => :v)
    sort!(g, :t)
    return collect(g.t), collect(g.v)
end

t, u_ll = series("ll", :u)
_, P_ll = series("ll", :P)
_, P_cc = series("cc", :P)
_, u_cond = series("cond", :u)
dP_ll = P_ll .- P_cc

println("結果: ", RES_DIR, "  (t = ", t[1], " 〜 ", t[end], " s, ", length(t), " 点)")
println()
println("窓 ", W, " s ごとの、液管の平均速度と、液管の P−P_cc (セル平均)、凝縮器の平均速度:")
println("   窓 [s]            u_ll 平均    σ(u_ll)   最大−最小    ΔP_ll 平均   σ(ΔP_ll)    u_cond 平均  σ(u_cond)")
stds = Float64[]; mids = Float64[]
t0 = t[1]
for k in 0:Int(floor((t[end] - t0) / W)) - 1
    win_a = t0 + k * W; win_b = win_a + W
    idx = findall(τ -> win_a <= τ < win_b, t)
    length(idx) < 3 && continue
    s_u = std(u_ll[idx])
    push!(stds, s_u); push!(mids, 0.5 * (win_a + win_b))
    mark = (win_a <= T_STEP < win_b) ? "  ← ステップ" : ""
    println("   ", lpad(round(Int, win_a), 4), "〜", lpad(round(Int, win_b), 4), "   ",
            lpad(round(mean(u_ll[idx]), digits = 5), 10), lpad(round(s_u, sigdigits = 3), 11),
            lpad(round(maximum(u_ll[idx]) - minimum(u_ll[idx]), sigdigits = 3), 11), "   ",
            lpad(round(mean(dP_ll[idx]), digits = 1), 10), lpad(round(std(dP_ll[idx]), sigdigits = 3), 11), "    ",
            lpad(round(mean(u_cond[idx]), digits = 5), 10), lpad(round(std(u_cond[idx]), sigdigits = 3), 11), mark)
end

# 減衰する局面の、e-fold 時間: ln σ の直線近似。floor (減衰後に残るゆらぎ) の 2 倍を超える窓だけを使う。
# floor は、準定常 (T_QUIET 以降) の窓の σ の中央値で推定する (ステップの過渡の窓は除く)。
function decay_fit(mids, stds, first_idx, floor_σ; label = "")
    # first_idx から始めて、σ が floor の DECAY_FLOOR_FACTOR 倍以下になる手前までの、連続した窓
    last_idx = first_idx
    while last_idx < length(stds) && stds[last_idx+1] > DECAY_FLOOR_FACTOR * floor_σ && stds[last_idx+1] > 0
        last_idx += 1
    end
    if last_idx - first_idx + 1 < 3
        println("  ", label, ": 減衰する局面の窓が足りません (", last_idx - first_idx + 1, " 個)")
        return nothing
    end
    fx = mids[first_idx:last_idx]; fy = log.(stds[first_idx:last_idx])
    fxm = mean(fx); fym = mean(fy)
    fit_b = sum((fx .- fxm) .* (fy .- fym)) / sum((fx .- fxm) .^ 2)
    r2 = 1 - sum((fy .- (fym .+ fit_b .* (fx .- fxm))) .^ 2) / sum((fy .- fym) .^ 2)
    println("  ", label, ": t_mid = ", round(fx[1], digits = 1), " 〜 ", round(fx[end], digits = 1), " s (", length(fx),
            " 窓), log σ の傾き ", round(fit_b, sigdigits = 3), " 1/s, R² = ", round(r2, digits = 3))
    length(fx) < 5 && println("    (窓が ", length(fx), " 個と少ないので、e-fold 時間の誤差は大きい。目安として扱うこと)")
    if fit_b < 0
        println("    → e-fold 時間 τ = ", round(-1 / fit_b, digits = 1), " s (層流摩擦だけの理論値 ", FRICTION_TAU, " s の ",
                round(-1 / fit_b / FRICTION_TAU, digits = 1), " 倍)")
    else
        println("    → 減衰していない")
    end
    return -1 / fit_b
end

println()
quiet = findall(i -> mids[i] > T_QUIET && !(T_STEP - W <= mids[i] <= T_STEP + 5 * W), eachindex(mids))
floor_σ = isempty(quiet) ? NaN : median(stds[quiet])
println("floor (準定常の σ(u_ll) の中央値; t>", T_QUIET, " s、ステップの過渡を除く) = ", round(floor_σ, sigdigits = 3), " m/s")
println("減衰する局面の e-fold 時間 (floor の ", DECAY_FLOOR_FACTOR, " 倍を超える窓だけを、近似):")
i_start = findfirst(i -> mids[i] > 0, eachindex(mids))
decay_fit(mids, stds, i_start, floor_σ; label = "起動直後")
i_step = findfirst(i -> mids[i] > T_STEP + TAU_STEP + W / 2, eachindex(mids))    # ステップの立ち上がりの直後の窓から
i_step === nothing || decay_fit(mids, stds, i_step, floor_σ; label = "ステップ後")

# ステップ前の準定常 (起動の過渡の後) の、σ(u_ll)
pre = findall(i -> 40.0 < mids[i] < T_STEP - W, eachindex(mids))
if !isempty(pre)
    println()
    println("ステップ前の準定常 (t=40 〜 ", T_STEP - W, " s) の σ(u_ll): 平均 ", round(mean(stds[pre]), sigdigits = 3),
            " m/s (最小 ", round(minimum(stds[pre]), sigdigits = 3), ", 最大 ", round(maximum(stds[pre]), sigdigits = 3), ")")
    println("  平均の流速は、", round(mean(u_ll[findall(τ -> 40.0 < τ < T_STEP, t)]), digits = 5),
            " m/s (期待: 20 W で、約 0.0076 m/s)")
end
# ---------------------------------------------------------------------------
# 圧力の突発的な増加と、凝縮器の前線のセルの状態の変化の照合
# ---------------------------------------------------------------------------
println()
println("="^76)
println("事象の照合: 液管の圧力の突発的な増加 と 凝縮器の前線のセルの状態の変化")
println("="^76)
cond = prof[prof.comp .== "cond", :]
gc = combine(groupby(cond, :t),
             :x => (x -> count(<(X_LIQ), x)) => :n_liq,       # 液 (飽和液線の下) のセル数
             :x => (x -> count(>(X_DRY), x)) => :n_dry,       # 乾き蒸気のセル数
             :x => (x -> count(xx -> X_LIQ <= xx <= X_DRY, x)) => :n_2p)    # 二相のセル数
sort!(gc, :t)
println("凝縮器 (", maximum(cond.cell), " セル) の状態: 液 x<", X_LIQ, ", 二相, 乾き x>", X_DRY)
println("  t=", round(gc.t[1], digits = 0), " s: 液 ", gc.n_liq[1], ", 二相 ", gc.n_2p[1], ", 乾き ", gc.n_dry[1])
println("  t=", round(gc.t[end], digits = 0), " s: 液 ", gc.n_liq[end], ", 二相 ", gc.n_2p[end], ", 乾き ", gc.n_dry[end])

# 状態の変化 (液セル数または乾きセル数が変わった時刻)。起動直後 (T_QUIET まで) と、ステップの過渡は除く
changes = Float64[]
for i in 2:nrow(gc)
    (gc.n_liq[i] != gc.n_liq[i-1] || gc.n_dry[i] != gc.n_dry[i-1]) || continue
    (gc.t[i] < T_QUIET || T_STEP - 2 <= gc.t[i] <= T_STEP + 30) && continue
    push!(changes, gc.t[i])
end
println("前線のセルの状態が変わった時刻 (液セル数または乾きセル数が変化。t>", T_QUIET, " s、ステップの過渡を除く): ", length(changes), " 回")
println("  ", join([string(round(c, digits = 0)) for c in changes[1:min(end, 40)]], ", "), length(changes) > 40 ? ", ..." : "")

# 圧力の突発的な増加: ΔP_ll の、中央値からのずれが、閾値を超えたサンプル
sel_t = findall(τ -> τ > T_QUIET && !(T_STEP - 2 <= τ <= T_STEP + 30), t)
dp_med = median(dP_ll[sel_t])
dp_mad = 1.4826 * median(abs.(dP_ll[sel_t] .- dp_med))
thr = max(SPIKE_MIN, SPIKE_FACTOR * dp_mad)
spikes = [τ for (τ, d) in zip(t[sel_t], dP_ll[sel_t]) if abs(d - dp_med) > thr]
println("液管の圧力 ΔP_ll の、突発的な増加: 中央値 ", round(dp_med, digits = 1), " Pa, MAD に基づく σ ", round(dp_mad, digits = 1),
        " Pa → 閾値 ", round(thr, digits = 1), " Pa; 閾値を超えたサンプル ", length(spikes), " 個")
println("  (saveat=1 s なので、継続時間が 1 s より短いスパイクは、サンプルに当たらないことがある)")
if !isempty(spikes)
    println("  圧力の増加の時刻, ΔP_ll [Pa], 最も近い前線の状態変化までの時間 [s]:")
    n_match = 0
    for τ in spikes
        d = dP_ll[findfirst(==(τ), t)]
        dt_near = isempty(changes) ? NaN : minimum(abs.(changes .- τ))
        matched = dt_near <= MATCH_WINDOW
        n_match += matched
        println("    t = ", lpad(round(τ, digits = 0), 6), " s   ΔP = ", lpad(round(d, digits = 0), 6), "   最も近い状態変化まで ",
                lpad(round(dt_near, digits = 1), 6), " s", matched ? "   ← 対応" : "")
    end
    println("  → ", n_match, " / ", length(spikes), " 個の圧力の増加が、前線の状態変化の ±", MATCH_WINDOW, " s 以内")
    # 比較の基準: ランダムな時刻が、状態変化の ±MATCH_WINDOW 以内に入る確率
    if !isempty(changes)
        span = maximum(t[sel_t]) - minimum(t[sel_t])
        p_rand = min(1.0, length(changes) * 2 * MATCH_WINDOW / span)
        println("     (参考: 時刻をランダムに選んだとき、状態変化の ±", MATCH_WINDOW, " s 以内に入る確率は約 ",
                round(100 * p_rand, digits = 0), " %。これより、はるかに高ければ、対応がある)")
    end
end

println()
println()
println("読み方: ・減衰する局面の e-fold 時間が、摩擦だけの理論値 (約 1.3 s) より長いのは、減衰が、液管の摩擦だけで決まらないため")
println("        ・floor (減衰後に残るゆらぎ) が、平均流速の数割あるのは、持続する励起があるため")
println("        ・圧力の増加が、前線のセルの状態変化と対応していれば、前線のセルが、状態方程式の硬さの変わる所を")
println("          またぐ事象が、液柱を叩いている (仮説)。対応が、ランダムな確率と同程度なら、別の原因")