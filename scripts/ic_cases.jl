using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# ic_cases.jl … 初期の気液分布を変えた起動計算 (主張 1〜4 の発表用。2026-11-03)
#   使い方: julia ic_cases.jl            (全ケースを順に)
#           julia ic_cases.jl C1 C3      (指定したケースだけ。2 つの端末で分けて並列に走らせてよい)
#   置き場所: scripts/ (transient_flat.jl と同じフォルダ)。transient_flat.jl は SWEEP に
#             n_vap_vl_init, f_l_grv_init, dT_onset を受け付ける版 (11/3 版) を使うこと。
#   しくみ: h_sweep.jl と同じ。Main.SWEEP を定義して transient_flat.jl を 1 回ずつ include する (stage = :short)。
#   結果: res/transient_flat_ammonia_ic_<ケース>/ (summary.csv, profiles.csv, ...)。前の結果は上書きしない (タグが違う)。
#   ログ: 各ケースの標準出力の「封入量の収支」(CC の液体積分率が 0〜1 に収まっているか) を必ず確認する。
#
#   共通の条件:
#     ΔT_onset = 0.55 K   … 界面蒸発のゲートの過熱度。実測 (1-b v2) の降下直前の受熱面 P の過熱度 0.53〜0.58 K の中間。
#                            r_c_wick (フィッティング値 1.85 µm) は使わない (transient_flat.jl で換算)。
#     f_l_grv_init = 0.05 … グルーブの初期液体積分率 (溝の体積の 5 %、隅に液がある程度)。本体の f_full_grv (0.08) 未満。
#     ほかは transient_flat.jl の現行の値のまま (h_ev 35000、h_con 16500、T_amb 21.5 °C、C_ec_P 8.374 など)。
#   ケース (蒸気管 20 セル。VL1 はセル 3〜4、VL2 は 10、VL3 は 17 の位置):
#     C1     蒸気管はほぼ全部蒸気 (20 セル)          … 極端
#     C2     実測の加熱前の壁温から (11 セル。現行)   … 基準
#     C3     蒸気管はほぼ全部液 (蒸発器側 2 セルだけ蒸気) … 極端
#     C2dry  C2 でグルーブを蒸気のみ (f_l_grv = 0)    … グルーブの液の有無の感度
#     C2lo / C2hi  C2 で ΔT_onset 0.53 / 0.58 K       … 実測の過熱度の幅 (時間に余裕があれば)
#     C4     C2 でグルーブを満水 (grv_init = :liquid、過冷却 0.1 K の液)   … 蒸発器側の極端 (11/4 追加)
#     C4h    C2 でグルーブの液 0.5 (飽和の気液二相)                       … 1-c で整合した液量 (溝の 14〜53 %) の上端 (11/4 追加)
#     C2lump C2 で CC を集中定数 (液 1 層、ほぼ気液平衡の CC) にしたもの (11/4 追加。発表の Results の「CC の非平衡」の比較用)
#            【実行の前に】src/types.jl の `const CC_LUMPED = false` を true に書き換え、Julia を再起動してから
#             julia ic_cases.jl C2lump。終わったら false に戻す (解析チャットからは src/ を書き換えないので、手元で一時的に)。
#            CC_LUMPED の状態とケースが合わないと実行しない (C2lump は true のときだけ、ほかは false のときだけ)。
#     注: 本体のグルーブは成層の出口で、液が f_full_grv (0.08) を超えると液も出口から出る。満水・0.5 は、流れが始まると
#         超えた分が蒸気管へ押し出される (= 蒸発器側の液が蒸気管の液に加わる) ことを意味する。計算が不安定になる可能性がある。
# =============================================================================

ic_t_end = 100.0          # 計算の終わり [s] (加熱開始は 10 s → τ = 90 s まで。CON1 の急上昇・P_RES のピークは τ ≈ 47 s)
ic_SKIP_EXISTING = true   # summary.csv が保存済みのケースは飛ばす (途中で止まったときの再開用)

ic_common = (stage = :short, t_end = ic_t_end, dT_onset = 0.55, f_l_grv_init = 0.05)
ic_all = [
    ("C1",    (n_vap_vl_init = 20,)),
    ("C2",    (n_vap_vl_init = 11,)),
    ("C3",    (n_vap_vl_init = 2,)),
    ("C2dry", (n_vap_vl_init = 11, f_l_grv_init = 0.0)),
    ("C2lo",  (n_vap_vl_init = 11, dT_onset = 0.53)),
    ("C2hi",  (n_vap_vl_init = 11, dT_onset = 0.58)),
    ("C4",    (n_vap_vl_init = 11, grv_init = :liquid, f_l_grv_init = 0.0)),
    ("C4h",   (n_vap_vl_init = 11, f_l_grv_init = 0.5)),
    ("C2lump", (n_vap_vl_init = 11,)),
]
ic_select = isempty(ARGS) ? [c[1] for c in ic_all if !(c[1] in ("C4", "C4h", "C2lump"))] : ARGS   # C4・C4h・C2lump は名前を指定したときだけ

# CC_LUMPED (src/types.jl の const) の状態を確かめる。ケース名が "lump" で終わるものは true のときだけ、ほかは false のときだけ実行する
# (CC の層数が違う結果を、同じ名前のフォルダに混ぜないため)
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel
ic_lump = LHPTransientModel.CC_LUMPED
println("CC_LUMPED = ", ic_lump, ic_lump ? "  (CC は液 1 層: ほぼ気液平衡)" : "  (CC は液 16 層: 非平衡)")
ic_dir(name) = joinpath(@__DIR__, "..", "res", "transient_flat_ammonia_ic_" * name)

ic_status = Tuple{String,String,Float64}[]
for (name, pars) in ic_all
    name in ic_select || continue
    if endswith(name, "lump") != ic_lump
        msg = ic_lump ? "CC_LUMPED = true なので、lump でないケースは実行しない (src/types.jl を false に戻して Julia を再起動)" :
                        "CC_LUMPED = false なので実行しない (src/types.jl の CC_LUMPED を true にして Julia を再起動)"
        println("\n(飛ばす) ", name, ": ", msg)
        push!(ic_status, (name, "飛ばした: " * msg, NaN))
        continue
    end
    if ic_SKIP_EXISTING && isfile(joinpath(ic_dir(name), "summary.csv"))
        println("\n(飛ばす: 保存済み) ", name)
        push!(ic_status, (name, "保存済み", NaN))
        continue
    end
    println("\n", "="^70, "\n  ケース ", name, "  ", merge(ic_common, pars), "\n", "="^70)
    global SWEEP = merge(ic_common, pars, (tag = "_ic_" * name,))
    t0 = time()
    try
        include(joinpath(@__DIR__, "transient_flat.jl"))
        push!(ic_status, (name, string(sol.retcode) * " (最終時刻 " * string(round(sol.t[end], digits = 2)) * " s)", time() - t0))
    catch err
        @warn "ケース $(name) の実行が途中で失敗" err
        push!(ic_status, (name, "失敗: " * sprint(showerror, err)[1:min(end, 200)], time() - t0))
    end
end

println("\n", "="^70, "\nまとめ (ケース | 結果 | 計算時間 [s]):")
for (n, s, dt) in ic_status
    println("  ", rpad(n, 7), " ", s, "  ", round(dt, digits = 0))
end
println("結果のフォルダ: ", abspath(joinpath(@__DIR__, "..", "res")), "/transient_flat_ammonia_ic_<ケース>/")
println("比較: python3 scripts/claims_ic_compare.py (Claude 側でも実行できる。summary.csv を渡す)")