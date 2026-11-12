# cc_mix_run.jl (未実行) : CC の液の混合 (噴流・プルーム) を入れた本体の :short (100 s) を、混合なし / 噴流 / プルーム / 両方 で順に走らせる。
# 使い方 (プロジェクトのルートで):  julia --project=. scripts/cc_mix_run.jl [none jet plume both]
#   引数なし = 4 つ全部。結果は res/transient_flat_ammonia<run_tag>/ (run_tag = _mix_none / _mix_jet / _mix_plume / _mix_both)。上書きしない。
# 比較の注意: SWEEP を定義すると transient_flat.jl は run_tag を SWEEP.tag にし、保存点に 13〜45 s の 0.25 s 刻みが加わる (物理は同じ)。
#   _mix_none は、従来の :short と一致するはずの確認用 (受け入れ条件 1)。
cases = isempty(ARGS) ? ["none", "jet", "plume", "both"] : ARGS
for cs in cases
    global SWEEP = (stage = :short, t_end = 100.0, tag = "_mix_" * cs, cc_mix = Symbol(cs))
    println("\n========== cc_mix = :", cs, " ==========")
    include(joinpath(@__DIR__, "transient_flat.jl"))
end