using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# test_hev_long.jl … h_ev = 2e5 W/m²K、f_h_shah = 0.4、起動時の過熱度 (HB − EC1) = 0.43922 K で :long (816 s) を 1 本計算する
#   使い方: julia test_hev_long.jl      (scripts/ に置く。transient_flat.jl は今回更新したものに差し替える。引数なし)
#   しくみ: param_sweep.jl と同じ。Main.SWEEP を定義して transient_flat.jl を include する。stage は transient_flat.jl の :long のまま (816 s)。
#   起動時の過熱度の与え方 (解釈。違っていたら教えてください):
#     「起動時の過熱度」= 界面蒸発が始まる (ゲート n が開く) 過熱度 ΔT_onset。モデルでは ΔT_onset = 2σT/(ρ_v h_fg r_c_wick) なので、
#     指定した ΔT_onset から r_c_wick を逆算して渡す (22.6 °C の飽和アンモニアの物性。この式で r_c_wick = 1.85 µm → 0.79 K、検算済み)。
#     ただしモデルのゲートは、受熱面 P (蒸発器ケース) の温度 − 飽和温度で判定する。HB − EC1 は HB から蒸発器への接触の温度差を含むので、
#     この ΔT_onset を HB − EC1 の値として扱うことになる (実験の起動時の HB − EC1 を ΔT_onset に置く解釈)。
#   結果: res/transient_flat_ammonia<tag>/ (summary.csv, profiles.csv, ...)。ログの末尾の実験との比較表を見る。
# =============================================================================

# ---------------------------- 設定 ----------------------------
th_h_ev     = 2.0e5        # [W/m²K] (A_if = 4.28e-4 m² 基準)
th_f_h_shah = 0.4          # 凝縮器の Shah の式にかける係数 (廃止していた f_h_shah を 0.4 に戻す)
th_dT_onset = 0.43922      # 起動時の過熱度 (HB − EC1) [K] = ゲートの ΔT_onset
th_tag      = "_hev2e5_fsh0.4_dT0.439"    # 結果フォルダ名の接尾辞
# ΔT_onset から r_c_wick への換算 (22.6 °C の飽和アンモニア)
th_σ    = 0.021037         # 表面張力 [N/m]
th_ρ_v  = 7.2546           # 飽和蒸気の密度 [kg/m³]
th_h_fg = 1.175738e6       # 蒸発潜熱 [J/kg]
th_T    = 273.15 + 22.6    # [K]
# ----------------------------------------------------------------

th_r_c = 2 * th_σ * th_T / (th_ρ_v * th_h_fg * th_dT_onset)
println("ΔT_onset = ", th_dT_onset, " K → r_c_wick = ", round(th_r_c * 1e6, digits = 3), " µm (検算: 1.85 µm → ", round(2 * th_σ * th_T / (th_ρ_v * th_h_fg * 1.85e-6), digits = 3), " K)")
println("h_ev = ", th_h_ev, ",  f_h_shah = ", th_f_h_shah)

global SWEEP = (t_end = 816.0, tag = th_tag, h_ev = th_h_ev, f_h_shah = th_f_h_shah, r_c_wick = th_r_c)
include(joinpath(@__DIR__, "transient_flat.jl"))