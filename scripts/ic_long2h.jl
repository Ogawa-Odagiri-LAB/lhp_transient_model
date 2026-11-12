using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# ic_long2h.jl … C2 と C2lump を :long2h (6910 s = 加熱開始 +6900 s) で実行する (2026-11-04)
#   目的: CC の液を 16 層 (C2) と 1 層 (C2lump、集中定数) にした 2 つの計算が、定常 (長時間) で同じ T_sat・RES に収束するかを調べる
#         (src/residuals_flatcc.jl の構造の確認から、定常の差は RES − T_sat で 0.1 K 以下と見積もった。それを確かめる)。
#   使い方: julia ic_long2h.jl      (scripts/ に置く。transient_flat.jl と同じフォルダ。引数なし)
#   どちらを走らせるかは、src/types.jl の `const CC_LUMPED` で決まる:
#       CC_LUMPED = false → C2      (CC の液 16 層)      結果: res/transient_flat_ammonia_2h_ic_C2/
#       CC_LUMPED = true  → C2lump  (CC の液 1 層)       結果: res/transient_flat_ammonia_2h_ic_C2lump/
#   CC_LUMPED は const なので、切り替えたら Julia を再起動する。src/ を書き換えると再コンパイルが入るので、
#   2 つを同時に走らせるときは、プロジェクトのフォルダごとコピーして、コピー側だけ CC_LUMPED = true にする (下の「同時に走らせる」)。
#
#   ケースの条件 (ic_cases.jl の C2 と同じ。短い計算 (:short、τ = 90 s) の C2 は res/transient_flat_ammonia_ic_C2/):
#       蒸気管の蒸気セル 11 (蒸発器側から)、グルーブの初期液 0.05、界面蒸発のゲートの過熱度 ΔT_onset = 0.55 K (実測の 1-b v2)。
#       ほかは transient_flat.jl の現在の値 (h_ev 35000、h_con 16500、足 G_foot_cc 0.128、孔内膨張 pore_expansion など)。
#   計算時間: :long2h は 4 時間程度 (C2lump は層が少ないので少し速い見込み。未測定)。
#   注意: transient_flat.jl は SWEEP に dT_onset、n_vap_vl_init、f_l_grv_init を受け付ける版 (11/3 版) であること。
#         そうでない版だと、SWEEP の値が黙って無視され、別のケースを走らせてしまう。下で、ファイルの中身を調べて止める。
#
#   同時に走らせる (2 つの端末):
#     1. プロジェクトのフォルダ (src/、scripts/、data/ を含む) を、隣に lhp_transient_model_lump としてコピーする。
#     2. コピー側の src/types.jl の `const CC_LUMPED = false` を true に書き換える。
#     3. 元のフォルダで julia ic_long2h.jl (C2)、コピー側の scripts/ で julia ic_long2h.jl (C2lump) を、別の端末で実行する。
#     4. 終わったら、コピー側の res/transient_flat_ammonia_2h_ic_C2lump/ を、元の res/ に移す。
# =============================================================================

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel

# ---------------------------- 設定 ----------------------------
ic2_t_end = 6910.0          # 計算の終わり [s] (:long2h の既定。加熱開始 10 s → τ = 6900 s)
ic2_SKIP_EXISTING = true    # 同じケースの summary.csv が保存済みなら飛ばす (途中で止まったときの再開用)。やり直すなら false か、フォルダを移動する
# ----------------------------------------------------------------

ic2_lump = LHPTransientModel.CC_LUMPED
ic2_name = ic2_lump ? "C2lump" : "C2"
println("CC_LUMPED = ", ic2_lump, ic2_lump ? "  (CC は液 1 層: 集中定数。ほぼ気液平衡)" : "  (CC は液 16 層: 成層・非平衡)", "  →  ケース ", ic2_name, " (:long2h)")

# transient_flat.jl が必要な SWEEP のキーを受け付けるかを確かめる (受け付けない版だと値が黙って無視される)
ic2_src = read(joinpath(@__DIR__, "transient_flat.jl"), String)
for key in ("dT_onset", "n_vap_vl_init", "f_l_grv_init")
    occursin(":" * key, ic2_src) ||
        error("transient_flat.jl が SWEEP の ", key, " を受け付けません (11/3 版の transient_flat.jl に、ic_cases.jl の対応と同じ修正が必要です)。このまま走らせると、C2 ではないケースになります")
end

ic2_tag = "_ic_" * ic2_name     # 結果フォルダ: res/transient_flat_ammonia_2h<tag>/ (stage = :long2h のとき。:short の C2 (res/transient_flat_ammonia_ic_C2/) とは別)
ic2_res = joinpath(@__DIR__, "..", "res")
ic2_existing = isdir(ic2_res) ? [d for d in readdir(ic2_res) if startswith(d, "transient_flat_ammonia_2h") && endswith(d, ic2_tag) &&
                                 isfile(joinpath(ic2_res, d, "summary.csv"))] : String[]
if ic2_SKIP_EXISTING && !isempty(ic2_existing)
    println("\n(飛ばす: 保存済み) ", join(ic2_existing, ", "), "\n  やり直すときは ic2_SKIP_EXISTING = false にするか、そのフォルダを別の場所へ移す。結果のフォルダは上書きしない")
else
    ic2_pars = (stage = :long2h, t_end = ic2_t_end, dT_onset = 0.55, f_l_grv_init = 0.05, n_vap_vl_init = 11)
    println("\n", "="^70, "\n  ケース ", ic2_name, "  ", ic2_pars, "\n", "="^70)
    global SWEEP = merge(ic2_pars, (tag = ic2_tag,))
    t0 = time()
    include(joinpath(@__DIR__, "transient_flat.jl"))
    println("\n", ic2_name, " 完了: ", sol.retcode, " (最終時刻 ", round(sol.t[end], digits = 2), " s)、計算時間 ", round((time() - t0) / 60, digits = 1), " 分")
end