using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# ic_hev_sweep.jl … C2 の条件で、界面蒸発の熱伝達係数 h_ev だけを振る感度計算 (2026-11-06)
#   使い方: julia ic_hev_sweep.jl                 (全部を順に)
#           julia ic_hev_sweep.jl 15000 100000    (指定した値だけ。2 つの端末で分けて並列に走らせてよい)
#   置き場所: scripts/ (transient_flat.jl・ic_cases.jl と同じフォルダ)。
#             transient_flat.jl は h_ev を SWEEP で受け付ける版 (11/6 版) を使うこと。
#   しくみ: ic_cases.jl と同じ。Main.SWEEP を定義して transient_flat.jl を 1 回ずつ include する (stage = :short)。
#   結果: res/transient_flat_ammonia_ic_C2hev<値>/ (summary.csv, profiles.csv, ...)。
#         比較は python3 scripts/claims_ic_compare.py (フォルダ名から自動で拾う)。
#
#   条件: C2 と同じ (蒸気管 11 セル、グルーブの液 0.05、ΔT_onset 0.55 K、CC 16 層、:short t_end 100 s) で、h_ev だけ変える。
#     h_ev [W/m²K]: 15000、35000 (= C2 と同じ値)、55000、100000。界面積 A_ec_wick = 3.645e-4 m² なので、
#     h_ev·A = 5.5、12.8、20.0、36.5 W/K。
#   35000 は C2 と同じ計算なので、hev_RERUN_BASE = false なら計算せず、res/transient_flat_ammonia_ic_C2/ の結果を使う
#   (C2 は 11/3 の transient_flat.jl で計算。その後の変更は SWEEP のキーの追加だけで、既定の値は同じ)。
#   注意: h_ev は起動直後だけでなく全期間に効く (ウィック界面の蒸発 = ループの定常の熱抵抗にも効く)。
#         h_ev が変わると、ケースから界面への熱が変わるので、沸騰の開始後の EC・HB の温度と、蒸発の立ち上がりの速さが変わる。
#         沸騰の開始の判定 (ΔT_onset 0.55 K) は h_ev によらない (受熱面 P の温度で判定)。
# =============================================================================

hev_values = [15000.0, 35000.0, 55000.0, 61000.0, 100000.0]
hev_t_end = 100.0             # 計算の終わり [s] (加熱開始は 10 s → τ = 90 s まで)
hev_RERUN_BASE = false        # true なら 35000 も新しく計算する (結果は _ic_C2hev35000)
hev_SKIP_EXISTING = true      # summary.csv が保存済みの値は飛ばす (途中で止まったときの再開用)

hev_common = (stage = :short, t_end = hev_t_end, dT_onset = 0.55, f_l_grv_init = 0.05, n_vap_vl_init = 11)
hev_select = isempty(ARGS) ? hev_values : [parse(Float64, a) for a in ARGS]
hev_name(v) = "C2hev" * string(round(Int, v))
hev_dir(name) = joinpath(@__DIR__, "..", "res", "transient_flat_ammonia_ic_" * name)

# CC_LUMPED (src/types.jl の const) は false (CC 16 層) のときだけ実行する (C2 と同じ CC にするため)
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel
if LHPTransientModel.CC_LUMPED
    error("CC_LUMPED = true になっている。C2 と同じ CC (16 層) にするため、src/types.jl を false に戻して Julia を再起動すること")
end

hev_status = Tuple{String,String,Float64}[]
for v in hev_select
    v in hev_values || @warn "h_ev = $(v) は既定の値の一覧 ($(hev_values)) に無いが、指定どおり計算する"
    name = hev_name(v)
    if v == 35000.0 && !hev_RERUN_BASE
        f = joinpath(hev_dir("C2"), "summary.csv")
        println("\n(計算しない) h_ev = 35000 は C2 と同じ: ", isfile(f) ? abspath(f) : "res/transient_flat_ammonia_ic_C2/summary.csv が無い (ic_cases.jl C2 を先に)")
        push!(hev_status, (name, "C2 の結果を使う", NaN))
        continue
    end
    if hev_SKIP_EXISTING && isfile(joinpath(hev_dir(name), "summary.csv"))
        println("\n(飛ばす: 保存済み) ", name)
        push!(hev_status, (name, "保存済み", NaN))
        continue
    end
    println("\n", "="^70, "\n  h_ev = ", v, " W/m²K  (", name, ")\n", "="^70)
    global SWEEP = merge(hev_common, (h_ev = v, tag = "_ic_" * name))
    t0 = time()
    try
        include(joinpath(@__DIR__, "transient_flat.jl"))
        push!(hev_status, (name, string(sol.retcode) * " (最終時刻 " * string(round(sol.t[end], digits = 2)) * " s)", time() - t0))
    catch err
        @warn "h_ev = $(v) の実行が途中で失敗" err
        push!(hev_status, (name, "失敗: " * sprint(showerror, err)[1:min(end, 200)], time() - t0))
    end
end

println("\n", "="^70, "\nまとめ (ケース | 結果 | 計算時間 [s]):")
for (n, s, dt) in hev_status
    println("  ", rpad(n, 12), " ", s, "  ", round(dt, digits = 0))
end
println("結果のフォルダ: ", abspath(joinpath(@__DIR__, "..", "res")), "/transient_flat_ammonia_ic_C2hev<値>/")
println("比較: python3 scripts/claims_ic_compare.py (C2 と並べて表示される)")