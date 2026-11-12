using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# 平板型蒸発器 (金属3Dプリンター製、secondary wick 付き) の過渡計算。
# transient.jl (円筒型) と同じ流れで、ジオメトリだけが違う。
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel
using Serialization

# --- 結果の保存先 ---
# 名前だけ指定すると ../res/<名前>/ に保存される。絶対パスを指定すると、そこに保存される。
# 結果の保存先は、下の stage の設定の後で決める (stage と Q_cc_fluid に応じた名前。複数の計算を同時に走らせても上書きしない)

# 凝縮器プレートの格子 (既定 30 × 6 = 10 mm 角)。grid_sweep.jl が SWEEP に plate_nx, plate_ny を入れたときだけ置き換える (数値解像度の確認用。コールドプレートも同じ格子)
plate_nx_run = (isdefined(Main, :SWEEP) && hasproperty(Main.SWEEP, :plate_nx)) ? Int(Main.SWEEP.plate_nx) : 30
plate_ny_run = (isdefined(Main, :SWEEP) && hasproperty(Main.SWEEP, :plate_ny)) ? Int(Main.SWEEP.plate_ny) : 6
# 凝縮器の管の分割数 (既定 50 = 11.2 mm/セル)。cond_N_sweep.jl が SWEEP に cond_N を入れたときだけ置き換える (数値解像度の確認用)。チラー側 (ColdPlateGeometry の N) は変えない
cond_N_run = (isdefined(Main, :SWEEP) && hasproperty(Main.SWEEP, :cond_N)) ? Int(Main.SWEEP.cond_N) : 50

# --- 幾何 (実機) ---
# primary wick φ=0.30、secondary wick φ_s=0.70 (どちらも液で満たされている)。
geo = FlatLHPGeometry(
    FlatWickGeometry(N=3, φ=0.30, ρ_wick=7929.0, Cp_wick=470.0),         # primary wick (高さ方向 N 層)。固体は SUS304: ρ, Cp は physics.jl の ρ_sus304, c_sus304 の T_init (22.6 °C) での値
    SecondaryWickGeometry(φ_s=0.70, ρ_swick=7929.0, Cp_swick=470.0),      # secondary wick (1ノード、常に液で満たされる)
    FlatGrvGeometry(N=1),                                                 # グルーブ (平板)
    TLGeometry(N=20, L=1.17, D_i=0.00176, D_o=0.00318),                   # 蒸気管
    CondGeometry(N=cond_N_run, L=0.56, D_i=0.00176, D_o=0.00318),                 # 凝縮器
    TLGeometry(N=20, L=1.15, D_i=0.00099, D_o=0.00159),                   # 液管
    FlatCCGeometry(N=1),                                                  # CC (半円筒 D_i=27 mm, L=25 mm)
    # 凝縮器プレート 300 × 60 mm、厚さ 5 mm (A5052。コールドプレートも同寸で全面接触)。管は短辺から入り、折り返し 1 回 (管の間隔 30 mm、曲げ半径 15 mm) で戻る。格子は 10 mm 角
    CondPlateGeometry(L=0.300, W=0.060, t=0.005, nx=plate_nx_run, ny=plate_ny_run, pitch=0.030),
    # コールドプレート (A5052、同じ平面寸法・格子、厚さ 10 mm)。チラー配管 (1/8 in = 外径 3.175 mm、エチレングリコール) は LHP の管と反対側の短辺から入り、凝縮器と同じ経路 (間隔 30 mm、長さ 0.56 m) で折り返す (対向流)
    ColdPlateGeometry(t=0.010, D_o=3.175e-3, L_tube=0.56, N=50),
)

tables = init_fluid_lookups("ammonia", T_min=230.0, T_max=360.0, P_max=4.5e6, Ne=1600, Δδ=0.5)   # 表の格子を細かく (Ne 400→1600, Δδ 1.0→0.5): 飽和蒸気付近の二相の T の誤差 +0.5 K → +0.13 K。初回はキャッシュ構築に時間がかかる
LHPTransientModel.P_SMOOTH_EPS[] = 0.01   # P(ρ,e) の飽和液線の折れをなだらかにする幅 [kg/m³]。ε=0.05/0.02/0.01 は :short 完走で結果ほぼ同じ (センサ 0.03 K 以内)、0.005 は失敗。飽和線上の圧力のずれは最大 約 4.6e5·ε [Pa] (ε=0.01 で 4.6 kPa ≈ 0.15 K)

# --- 実験のセンサ温度 [°C] ---
# (1) t = 0 の値: 計算の初期条件に使う。
sensors_t0 = (HB = 22.748695, EC1 = 22.73939, EC2 = 22.717699, VL1 = 22.777954, VL2 = 22.498681, VL3 = 21.907645,
              CON1 = 20.33362, CON2 = 20.333598, CON3 = 20.372776, CON4 = 20.374002, CON5 = 20.583229,
              LL1 = 22.163272, LL2 = 22.884449, LL3 = 22.882322, RES1 = 22.82656, RES2 = 22.739005,
              cond = 20.35)   # cond: CON1〜4 の実測 (20.33〜20.37) に合わせた実効のシンク温度 (測定値 20.483 は CON より高い)
# (2) t = 290.76〜302.76 s (加熱前の7点) の平均: 加熱前の計算結果 (t = 303 s) との比較に使う。
sensors_ref = (HB = 22.6508, EC1 = 22.6472, EC2 = 22.5981, VL1 = 22.6013, VL2 = 22.4894, VL3 = 21.7551,
               CON1 = 20.3722, CON2 = 20.3682, CON3 = 20.3917, CON4 = 20.4093, CON5 = 20.6119,
               LL1 = 21.5653, LL2 = 22.3120, LL3 = 22.6705, RES1 = 22.6417, RES2 = 22.5700)

# 液管のセンサ位置 (凝縮器側の根元から、全長 1.15 m に対する比): LL1 = 凝縮器から約 10 cm、LL2 = CC から約 50 cm (推定)、LL3 = CC から約 4 cm (実測)
s_LL1, s_LL2, s_LL3 = 0.10 / 1.15, (1.15 - 0.50) / 1.15, (1.15 - 0.04) / 1.15

# --- 計算の段階 ---
#   :short … 加熱前 (290〜302 s) の実測の平均 sensors_ref から、10 s で 5 W を印加して 100 s まで (起動の速さ・HB/EC の昇温を実験と比べる試験)。
#            (t = 0 の値は、RES が加熱前より約 0.2 K、LL が 0.3〜0.6 K 高く、実験の加熱前の状態とは違う)
#   :pre   … t = 0 の実測から 303 s まで、加熱前の計算。終わると状態を保存する。
#   :heat  … 保存した t = 303 s の状態から、5 W を印加して 320 s まで。
# τ_nuc の感度試験 (tau_nuc_sweep.jl) が Main.SWEEP を定義したときだけ、短い計算 (:short、t_end = SWEEP.t_end) にして τ_nuc を SWEEP.τ_nuc に置き換える。単独で実行するときは何も変わらない
sweep = isdefined(Main, :SWEEP) ? Main.SWEEP : nothing
t_short_end = sweep === nothing ? 100.0 : sweep.t_end
# 感度試験で上書きできる値 (SWEEP に同名の項目があるときだけ置き換える): τ_nuc, r_c_wick, h_con, h_ev, f_h_shah, h_hb_ec, h_cp, T_amb, cc_adiabatic, f_ec_cc, stage (幾何は plate_nx, plate_ny。上の plate_nx_run 参照)。tag があれば結果フォルダ名の接尾辞にする
sweep_get(name::Symbol, default) = (sweep !== nothing && hasproperty(sweep, name)) ? getproperty(sweep, name) : default
stage_default = :long   # 【単独実行の段階】 :short / :long (816 s) / :long2h (6910 s)。ここを書き換える。感度試験・フィットなどが SWEEP に stage を入れたときはそちらが優先
stage = sweep_get(:stage, stage_default)   # :short / :long (816 s) / :long2h (6910 s)。2 つ同時に走らせるときは、このファイルをコピーして stage だけ変える (結果は stage ごとの別フォルダに保存される)
Q_cc_fluid_run = 0.0   # 【診断用】CC の流体への一定の熱 [W]。0 = なし (CC ケースを 2 節点にして、蒸発器ケースからの伝導で説明できるかを見る)
G_br_run = 0.01   # 熱橋のコンダクタンス [W/K]。承認済み (0.01)。実験の LL1→LL3 の昇温 (+1.33 K @90 s) に対してモデル +1.20 K。0.05 では +1.97 K で過大。別の値を試すときは結果のフォルダ名に _Gbr<値> が付く
gate_case = :P   # :P = 受熱面 P で判定 + r_c_wick 1.5 µm (v5 と同じ。v6a)、:wick = ウィック最外セルで判定 + r_c_wick 5.1 µm (v5b と同じ。v6b)。2 つ同時に走らせるときは、このファイルをコピーして gate_case だけ変える
run_tag = sweep !== nothing ? (hasproperty(sweep, :tag) ? string(sweep.tag) : "_tau" * string(sweep.τ_nuc)) : (gate_case === :wick ? "_v9b" : "_v9a") * "_hcon_geo2_foot_pore" * (LHPTransientModel.CC_LUMPED ? "_lump" : "")   # 結果フォルダ名の接尾辞 (前の計算を上書きしないため)。v9a/b = v8e/f + secondary wick を 2 節点 (A: 平面側 + 両端面、B: 半円環 = 蒸気が凝縮する壁。swick の状態が 2 要素になるので :pre の再実行が必要)。v8e/f = v8c/d + 壁面凝縮の温度差の比 g_T (cc_wall_dT)。v8c/d = v8a/b + 壁面凝縮 (cc_wall_cond = true)。v8a/b = v7 の CC を 気相 1 セル + 液 16 層 (ccl。界面で凝縮 ΔP_ne>0、液管は最下層に流入、壁面凝縮なし) に置き換えた版。v7 = v6 + 熱の流れを P → S → CC ケース (側面 S が断熱材と CC ケースに熱を出す) + 壁面凝縮 (cc_wall_cond。w_v と f_v_ref_ne の代わりに、液面と secondary wick 壁面の熱浸透率 × 面積の比。壁面に入る熱は swick へ)。v6 =v5 / v5b + C_hb を 45 → 47.2 J/K (A5052 ブロック − 取付穴 + ヒーター + ねじ)。v5b =v5 + 開始の判定の温度を primary wick の最外セルにし (onset_ref = :wick)、r_c_wick を実測の 5.1 µm に戻した版。v5 =v4b + 蒸発器ケースの 2 節点 (P 6.28・S 3.48 J/K、G_ps 0.93 W/K、C_ec_extra 廃止、EC のセンサは S、ゲートは P) + r_c_wick 1.5 µm。v4b =v4a + ヒーター立ち上がりの区切り点を 0.1 s 刻みに (types.jl)。v4a = v3 から f_v_ref_ne を 0.2 に戻した版 (v3 = r_c_wick 5.1 µm、f_l_grv_init = 0、h_con 9000 (f_con_area 1.0)、h_chiller = h_con) + CC ケース 2 節点 (cc_2node)、f_ec_cc = 1.0、C_ec_extra = 18、Q_cc_fluid = 0
result_name = (stage === :long2h ? "transient_flat_ammonia_2h" : "transient_flat_ammonia") * (G_br_run == 0.01 ? "" : "_Gbr" * string(G_br_run)) * run_tag
result_dir = isabspath(result_name) ? result_name : joinpath(@__DIR__, "..", "res", result_name)
#   :long2h … :long と同じ初期状態から 6910 s まで (加熱開始 +6900 s。実験は exp_5W_2h.csv = 2 時間の実測 (加熱開始 +6896 s まで)。計算に 4 時間ほどかかる)。
#   :long  … :short と同じ初期状態 (加熱前 290〜302 s の実測の平均) から、10 s で 5 W を印加して 816 s まで (12 分の実行。実験は exp_5W_12min.csv と比べる)。
#            :pre → :heat は、:pre の 303 s の間に液管末端の二相セル (f_v_ll_init) が消えて IDA が止まるため、使わない。
init_sensors = stage in (:short, :long, :long2h) ? merge(sensors_ref, (cond = sensors_t0.cond,)) : sensors_t0

# 初期状態: グルーブは蒸気 (x = 1)、蒸気管は蒸発器側が蒸気・凝縮器側が液 (前線は壁温から決まる: 壁が飽和温度より 0.15 K 以上低いセルが液 → セル 12〜20)、凝縮器・液管は液。
grv_init_state = sweep_get(:grv_init, :vapor)   # :liquid で満水 【ic_cases.jl で上書き可】
# 封入量 7.6 g。θ_tol: 蒸気管の壁が飽和温度より何 K 以上低いと液とみなすか (前線の位置を決める)。
# 実験の加熱前は VL1 ≈ T_sat (Θ ≈ 0)、VL2 は −0.12 K、VL3 は −0.85 K。VL3 は蒸気に触れていれば成り立たない温度 (液)。
# VL1・VL2 は蒸気でも矛盾しない (+22 s の P 上昇で飽和温度が +0.36 K 上がり、VL2 が +25 s に 22.48 → 22.96 に跳ねた)。θ_tol = 0.15 で前線は VL2 と VL3 の間 (セル 11/12)。
# 以前は θ_tol = 0.0 (蒸気管はすべて液)。蒸気管を蒸気にすると CC に割り当てる質量が増える (前線 1 セルで 0.085 g、グルーブ 0.2 g): 出力の「封入量の収支」を確認。
# n_vap_vl_init: 蒸気管の蒸気セル数 (蒸発器側から。VL1 はセル 3〜4、VL2 はセル 10、VL3 はセル 17 の位置)。nothing なら θ_tol で判定 (θ_tol=0.15 で 11 セル)。
n_vap_vl_init = sweep_get(:n_vap_vl_init, nothing)    # 【ic_cases.jl で上書き可】 実験データ (VL2 が +20〜25 s に飽和温度に追従して跳ねる) から、蒸気の範囲は θ_tol=0.15 の 11 セルのまま
# θ_sat_ll_init: 液管の壁温が T_sat − θ_sat_ll_init [K] 以上のセル (CC に近い側) を、飽和に近い液 (過冷却度 0.02 K) にする。0 なら使わない。
# 実験の LL3 は加熱前に 22.67 ℃ (RES 22.59〜22.64、T_sat 22.61) で、液なら過冷却のはずが、飽和以上を示すため。センサ誤差 0.1 K として 0.1。
θ_sat_ll_init = 0.1
# f_v_ll_init: 液管の CC 側の末端 n_ll_2ph_init セル (1 = LL3 を含む最終セル) を、蒸気の体積分率 f_v の飽和二相にする (0 なら使わない)。
# 【仮】0.10。LL3 が圧力の飽和温度に沿うのに必要な量は 8% (液の熱容量のみ) 〜 25% (壁を含む)。1 セル = 5.75 cm × 内径 0.99 mm = 44 mm³ なので、0.10 で蒸気スラグ約 6 mm。
f_v_ll_init = 0.25
# x_vap_vl_init: 蒸気管の初期蒸気セルの乾き度 (1 = 表の単相蒸気で飽和 +0.7 K。1 未満なら飽和の二相。まず 0.9995 で初期ログの T を確認する)
x_vap_vl_init = 0.9995
n_ll_2ph_init = 1
f_l_grv_init = sweep_get(:f_l_grv_init, 0.0)         # 【ic_cases.jl で上書き可】 グルーブの初期液体積分率 (0 = 蒸気のみ)。f_full_grv (0.08) より十分小さくし、加熱前の液の流出 (早い押し出し) を避ける。0.3 では t=1 s に P_cc が +2 kPa 跳ね、0.08 では s_full≈0.47 で +10 s に s=0.08 m
init = measured_initial(geo, tables, init_sensors; charge = 7.6e-3, θ_tol = 0.15, grv_init = grv_init_state,
                        n_vap_vl = n_vap_vl_init, f_l_grv = f_l_grv_init, θ_sat_ll = θ_sat_ll_init, f_v_ll = f_v_ll_init, n_ll_2ph = n_ll_2ph_init)

# 界面蒸発のゲートの過熱度 ΔT_onset [K] を直接与える (ic_cases.jl が SWEEP に dT_onset を入れたときだけ)。
# 本体のゲートは ΔT_onset = 2σT/(ρ_v h_fg r_c_wick) (上の核沸騰の診断と同じ式) なので、初期の飽和温度 T_init の物性で r_c_wick に換算して渡す。
# 値は実測 (1-b v2: 降下直前の受熱面 P の過熱度 0.53〜0.58 K) から。与えないときは従来どおり r_c_wick (1.85 µm、フィッティング)。
dT_onset_run = sweep_get(:dT_onset, NaN)
r_c_wick_run = if isnan(dT_onset_run)
    sweep_get(:r_c_wick, gate_case === :wick ? 5.1e-6 : 1.85e-6)
else
    let T0 = init.T_init, sp0 = LHPTransientModel.sat_props(init.T_init, tables)
        2 * sp0.σ * T0 / (sp0.ρ_v * sp0.h_fg * dT_onset_run)
    end
end
println("初期の気液分布: 蒸気管の蒸気セル = ", n_vap_vl_init === nothing ? "θ_tol で判定" : n_vap_vl_init, ",  グルーブの液体積分率 f_l_grv_init = ", f_l_grv_init,
        ",  ゲートの過熱度 ΔT_onset = ", isnan(dT_onset_run) ? "r_c_wick から" : string(dT_onset_run) * " K (実測から)", " → r_c_wick = ", round(r_c_wick_run * 1e6, digits = 3), " µm")

# 凝縮器: LHP の管 (SUS304) → 凝縮器プレート (A5052) → コールドプレート (チラー、20.6 °C 一定)。
# 実測の CON は境界条件に使わない (管壁・流体の初期値と、比較にだけ使う)。プレートはコールドプレートの温度から始める。
T_cold_plate = 20.4 + 273.15   # CON2〜4 の実測 (20.4 ℃) に合わせる (20.6 ℃ だと +0.2 K 高かった)
state_file = joinpath(result_dir, "state_t303.jls")

# 加熱の開始時刻と立ち上がり: :short は 10 s、:pre/:heat は実験の時刻 (303.8 s) に合わせる。
t_heat = stage in (:short, :long, :long2h) ? 10.0 : 303.8

# ヒーターから HB のセンサまでの遅れ [s]: 実験の HB は電力が入ってから約 4 s 動かない (電力 304.8 s → HB は 308.8 s から上昇)。
# 加熱開始 t_heat (= 実験の 303.8 s に対応) に対し、モデルの電力は heater_delay だけ遅れて入れる。比較表の時間軸は t_heat のまま。
heater_delay = 3.5

condition = LHPTransientCondition(;
    fluid = "ammonia",
    T_amb  = sweep_get(:T_amb, 21.5 + 273.15),        # 周囲温度 [K]。tamb_sweep.jl が 22.3 °C (室温の実測) と比べる
    Q_apply_init = 0.0,
    Q_apply = 5.0,
    t_step = t_heat + heater_delay, τ_step = 1.5,
    ΔT_startup = 0.0,
    # ヒーターブロック HB (A5052 53×33×10 mm ≈ 42 J/K) と蒸発器ケースの接触熱伝達。
    # 実験の起動後の HB−EVA1 = 0.85 K (5 W) から、G = h A ≈ 5〜6 W/K が必要。接触面積は蒸発器ケースの受熱面 (32×32 mm)。
    C_hb = 47.2,      # HB の熱容量 [J/K] (承認済み)。物理的な見積り: A5052 ブロック 53×33×10 mm = 42.2 J/K − ヒーターの取付穴 (Φ5×40 mm) のアルミ 1.9 + カートリッジヒーター (Φ5×40 mm、SUS シース+MgO) 2.5 + 固定用 M4 キャップスクリュー 4 本 4.4 (ご提示) ≈ 47.2 J/K。以前の 45 は加熱前半の当てはめ値だった。+15〜30 s の収支からの逆算は約 52 J/K で、差 (約 5 J/K) は未説明
    C_ec_extra = 0.0,  # (ec_2node = true のときは使わない) 蒸発器ケースの追加熱容量 [J/K]。1 節点のときの調整値は 18 J/K だった
    # 蒸発器ケースの 2 節点: P = 受熱面 (T_wall.evap。HB から入熱、ウィック・グルーブ・蒸気管・CC ケースへ出る。核沸騰のゲートの T_ec)、
    # S = 側面 (EC1・EC2 のセンサ。P から熱を受け、CC ケースと断熱材 (側面が接する) へ熱を出す: P → S → CC ケース)。
    # 熱容量は実測の見積り (ご提示): P 8.374 J/K (底板厚 2.0 mm、設計値 11/2。旧 6.28 は 1.5 mm 前提)、S 3.48 J/K (側壁 1.5 mm のまま。G_ec_ps も側壁の厚さだけに依存し底板厚には依存しないので 0.93 のまま)。
    # 固定用の M4 キャップスクリュー 4 本 (4.44 J/K) は先端が HB に留まっているので、HB 側 (C_hb = 47.2 J/K に含まれる) とみなして EC には足さない。
    # G_ec_ps: P–S の熱コンダクタンス = k A / L = 0.93 W/K (確定)。SUS304 k ≈ 15 W/mK、側壁の厚さ 1.5 mm、
    #   外周 32×2 + 33×2 = 130 mm → 肉厚中心の周長 2(32−1.5) + 2(33−1.5) = 124 mm、A = 124 × 1.5 = 186 mm²、
    #   受熱面からセンサまでの高さ L = 3 mm。外周 130 mm のまま計算すると 0.98 W/K。
    ec_2node = true, C_ec_P = 8.374, C_ec_S = 3.48, G_ec_ps = 0.93,
    # 以下の旧 CC 用の 4 つ (τ_ne_cc, K_ne_cc, f_v_min_ne, f_v_ref_ne) と cc_layer* / cc_meniscus / ne_liq_decouple / δ_layer_cc は、新しい CC (気相 1 + 液 16 層) では使われない (円筒モデル用に残してある)
    τ_ne_cc = 0.3,    # [旧 CC。未使用] CC の有限速度の凝縮の時定数 [s] (0 で無効)。大きいほど、液の流入による圧縮で CC の圧力が大きく上がる (目標: +42〜46 s で 8〜9 kPa)
    K_ne_cc = 1.2e6,  # 蒸気の体積弾性率 ≈ γP [Pa]
    f_v_min_ne = 0.05,
    f_v_ref_ne = 0.2,
    f_h_shah = sweep_get(:f_h_shah, 1.0),   # 【廃止: 1.0 = Shah の式そのまま (補正なし)。以前は 0.43 (Tsat−CON1 から) だったが、h_con と同じ CON1 を合わせるための 2 つ目の調整値だったので 1.0 に固定し、CON1 は h_con だけで合わせる】 元の説明: 蒸気管 (VL)・凝縮器・液管 (LL) の凝縮熱伝達 (Shah ≈ 2.4e4 W/m²K) にかける係数 (全体で統一)。実験の Tsat−CON1 (1 W あたり) から 0.25×1.72≈0.43。ε_phase_h を小さくした (静止した液体に凝縮の h が混ざらない) ので、0.43 でも収束するか確認する。LL にも同じ係数を掛ける (LL は従来 1.0)
    ε_phase_h = 0.002,   # 熱伝達の単相⇔二相のブレンド幅 (疑似乾き度。0.05 → 0.002 = 飽和から約 0.5 K)。0.05 だと過冷却度 2 K の静止した液体でも二相の重みが約 40% で、Shah の Nusselt 項 (≈4e4 W/m²K) が混ざり、f_h_shah>0.25 で t≈23 s に IDA が収束しなかった。摩擦のブレンド幅は変えない。おかしくなったら 0.05 に戻す (棄却)
    h_xps_mode = :P,   # 切り分け試験: :T にすると h のブレンドの疑似乾き度が従来 (T 基準) に戻る
    C_cc_extra = 1.9, # CC ケースの追加熱容量 [J/K]。幾何の V_w から 4.54 J/K、実測の見積り 6.4 J/K との差
    Q_cc_fluid = Q_cc_fluid_run,   # 【診断用・物理モデルではない】CC の流体へ一定の熱 [W] (t=0 から)。0.045 W = 2.4 K × G_sub (0.02 W/K)。遅い上昇 (800 s で +0.65 K、2 h で +2.4 K、τ≈2000 s) を再現するかを見る
    f_ec_cc = sweep_get(:f_ec_cc, 1.0),    # 蒸発器ケース→CC ケースの熱伝導にかける係数。1.0 = 幾何どおり (k_SUS304 と図面の断面・経路。FlatCCGeometry の W_o, t_shell, t_end)。以前の 0.1 は、CC ケースが 1 節点で、実機の熱の伝わる遅れ (拡散時間 約 1 分) を持たず、+20 s まで RES が動かない実験に合わせた値だった。2 節点にして遅れを持たせた
    pore_expansion = sweep_get(:pore_expansion, true),   # 孔内の液 (両 wick) の熱膨張を入れる (密度 = ρ_l,sat(T)。押し出された分は CC の液へ戻る)。false で従来
    G_foot_cc = sweep_get(:G_foot_cc, 0.128),   # 蒸発器ケースの足 (HB にボルト締め) を通る HB → CC ケースの熱コンダクタンス [W/K]。HB から CC ケースまでの全体の見積り (実測に基づく。ユーザー提示)。0 で従来
    cc_2node = true,  # CC ケースを 2 節点 (A: 蒸発器側 = C_cc_extra、B: 頂部側 = T_wall.cc、RES はこちら) にする
    h_hb_ec = sweep_get(:h_hb_ec, 5400.0), A_hb_ec = 0.032 * 0.032,   # G = 5.5 W/K
    # 蒸気管・液管と周囲の熱交換: :bare (裸管の自然対流) | :insulated (断熱材つき) | :adiabatic (断熱)
    amb_mode = :insulated, ins_k = 0.036, ins_D_in = 0.010, ins_D_out = 0.030,   # エアロフレックス (内径 10、外径 30 mm)
    # 熱橋: 平板型は CC・蒸発器・HB が一体で、液管の CC 側の末端が熱源に近い。熱源 (蒸発器ケース) から液管の CC 側の壁へ、
    # 末端から指数的に減る配分 (減衰長 br_decay セル) で n_br セルに入れる。
    # LL3 は CC から 4 cm (末端のセル)。液は CC へ向かうので、LL3 より上流のセルに入った熱が LL3 の液を温める (配分の約 9 割)。
    # 目安: LL3 を 22.5 °C にするのに、液の熱容量流量 0.018 W/K × 約 1.2 K ≈ 0.022 W、蒸発器ケース−液の温度差は約 1.7 K
    #  → 上流の G ≈ 0.013 W/K、全体では G_br ≈ 0.015 W/K (R ≈ 70 K/W)。
    # 実験の LL2 (CC から 50 cm) は +46 s に 21.3 °C (−1.0 K) まで冷え、LL3 (4 cm) は 22.5 °C。冷たい液が 46 cm で 1.2 K 温められる。
    # 起動中の液の流量は定常の約 3 倍 (0.02 mL/s) で、必要な熱は約 0.07 W → G_br を増やし、上流の約 30 cm (br_decay 5 セル) に配分する。
    G_br = G_br_run, n_br = 10, br_decay = 5.0, br_src = :ec,
    # 周囲との遅い熱経路 (residuals_housing.jl)。実験は +300 s 以降もループ全体が 0.9 mK/s で上がり続ける。CC は還流の過冷度で 0.02 W/K しか熱を吸収できず (50 K/W)、
    # 0.013 W の入熱の増加で T_cc が 0.67 K 動く。HB・蒸発器・CC の上側を包むエアロフレックス (厚さ 20 mm、時定数 数百 s) と、HB の下の PEEK ブロック (60 mm、拡散時間 2.6e4 s) が候補。
    #   【試行】 hi_on = peek_on = false で従来と同じ (熱を出さない)。面積 hi_A_*、ρc_p hi_rcp は目安 (要確認)。
    hi_on = true, cc_adiabatic = sweep_get(:cc_adiabatic, false), hi_t = 0.020, hi_rcp = 1.1e5,   # 面積 hi_A_hb = 2×32×10 mm、hi_A_ev = 2×32×5.5 mm、hi_A_cc = π×32/2×25 mm (側面だけ)
    peek_on = true, peek_L = 0.060, peek_k = 0.25,
    n_ll_ci = 2,   # 液管の CC 側の末端 2 セル (約 11.5 cm。LL3 のセンサは CC から 7 cm) は CC 側の断熱材と一続きの断熱材の中: 周囲温度を CC 側の断熱材の節点 (hi.T[3]) の温度にする (0 で従来)。要 hi_on = true
    # 蒸発の熱伝達係数 h_ev [W/m²K] = 50000 (据え置き。正確なフィッティングは今後)。界面積は A_ec_wick = フィン 18 本×0.75 mm×27 mm = 3.645e-4 m² に変えたので G = h_ev A = 18.2 W/K (旧: A = 4.28e-4 m² で 21.4 W/K。旧幾何の結果と比べるときは G で揃える: h_ev(新) = h_ev(旧) × 1.174)。以下は旧幾何 (A_if = 4.28e-4 m²) での決め方 (フィッティングではなく定常の熱収支):
    #   h_ev = Q / (A_if (T_EC − T_RES)),  Q = 5 W (印加熱)、T_EC = EC1・EC2 の平均、T_RES = RES1・RES2 の平均 (定常、τ≈6900 s で ΔT = 0.24 K)。
    #   → G = 20.7 W/K、h_ev = 48500 W/m²K を丸めて 50000。12 分 (τ≈800 s、ΔT = 0.19 K) では約 60000、Q = Q_evap (4.4 W) では約 43000〜53000。不確かさ ±20 % 程度 (HB・EC への影響は ±0.06 K)。
    #   同じ方法で h_hb_ec = 5400 も決まっている (Q = 5 W、A = 0.032²、T_HB − T_EC ≈ 0.92 K)。T_HB − T_RES で決めると HB→EC の抵抗が二重に入るので、h_ev には T_EC を使う。
    #   以前の注記: h_ev = 20000 (8.6 W/K) は、旧モデルで後半の EC−RES が 0.66 K (実験 約 0.15〜0.24 K) になり、熱橋の正のフィードバックで飽和温度の高さにも効いた。起動 (+20〜50 s) への影響は小さい。
    # h_ev = sweep_get(:h_ev, 50000.0), 
    h_ev = sweep_get(:h_ev, 35000.0),   # 【ic_hev_sweep.jl で上書き可。既定は従来どおり 35000】
    # 【試験用】h_ev の時間変化: h_ev_end を与えると h_ev (早い時刻) → h_ev_end (遅い時刻) へ τ = h_ev_tau_mid [s] を中心に滑らかに移る。既定 (NaN) は定数。
    h_ev_end   = sweep_get(:h_ev_end, NaN),
    h_ev_t_mid = t_heat + sweep_get(:h_ev_tau_mid, 50.0),
    h_ev_width = sweep_get(:h_ev_width, 15.0),
    h_con = sweep_get(:h_con, 16500.0),    # 凝縮器管 (SUS304) と凝縮器プレート (A5052) の間の熱伝達係数 [W/m²K]。基準面積は管の外周面積 × f_con_area。熱コンダクタンスは h_con × f_con_area = 12400 W/m²K (外周面積基準) で決まる。【計算値 (フィッティングではない)】SUS304 管 (k = 14.9 W/mK、D_o 3.18 mm、D_i 1.76 mm) の肉厚方向の厚肉円筒の伝導 (全周接触なら 15800 W/m²K) を、外周の下側約 3/4 だけが接触する断面の 2 次元定常熱伝導 (接触弧は等温・露出部は断熱、内面の h_i = 5000〜25000 W/m²K で ±2 %) で換算した 12100〜12600 W/m²K の代表値 12400 を、f_con_area = 0.75 で割ったもの (要望書_凝縮器の管壁コンダクタンスと凝縮熱伝達、2026-11-02)。含まれないもの: 拡管部の接触抵抗、プレート側の局所的な広がり抵抗 (どちらも実際を下げる方向)。以前の 20000 は平板近似 k/t (21000) に近く、約 26 % 大きかった
    f_con_area = 0.75,   # 管壁−プレートの接触面積 / 管の外周面積。管はプレートの上面側に埋め込まれ、上部はプレートから露出しており、接触は外周の下側約 3/4 (ユーザーの断面図からの判断値。正確な測定は困難)。CON1 の熱電対は管の上部の露出部に付いている (ユーザー確認、11/2)
    # --- グルーブの核沸騰 (グルーブが濡れている間は界面蒸発を止め、Chen の h で沸騰) ---
    grv_boiling = true,
    f_full_grv = 0.08,             # [0.3 → 0.08] 前の計算 (0.3): 起動後のグルーブの液割合が f_l = 0.26 (+25 s) → 0.217 (+50 s) で、乾かない (x = 0.043)。残る液は f_full − 約 0.04 (ε_full_grv = 0.03 の分) で、その後の減りは +25〜50 s で 0.043。+50 s に x = 1 (f_l → 0) にするには、+25 s の f_l ≈ 0.045 → f_full ≈ 0.085。成層出口で、液が出口のセルを f_l ≥ この値 満たすと液も出る (0.9 だと液体積分率 0.7 のまま溝に残った)。下げると液が抜けて溝が乾く
    grv_stratified = true,         # 成層した出口 (液は溝に残り、蒸気だけが出る)
    f_l_grv_init = f_l_grv_init,
    f_v_ll_init = f_v_ll_init, n_ll_2ph = n_ll_2ph_init,
    x_vap_vl_init = x_vap_vl_init,
    grv_init = grv_init_state,     # グルーブの初期状態 (:vapor = 蒸気 x = 1、:liquid = 満水)
    cc_layer = false,              # 【仮】一次遅れ版 (C_layer_cc, τ_layer_cc)。cc_layer_cond = true のときは使われない
    C_layer_cc = 2.0,  τ_layer_cc = 10.0,
    cc_layer_cond = true,          # 界面の液側の 1 次元熱伝導層 (面積 A_if は液深さから。フィット係数なし)
    f_ex_layer_cc = 0.035,         # 蒸気がなくなるときの ΔP_ne の消滅の尺度 (g_ex = f_v²/(f_v²+f_ex²))。実験の ΔP_ne の減衰から
    cc_meniscus = true,            # ΔP_ne の消滅の尺度を、蒸気キャップの高さ s と毛管長さ l_c から決める (g = s²/(s²+l_c²))。true のとき f_ex_layer_cc は使われない。物性値のみでフィット係数なし
    f_men_layer_cc = 1.0,          # g = 1/2 になる s / l_c (係数 1 = s が毛管長さに等しいとき)。f_ex_layer_cc = 0.035 は s ≈ 1.3 mm ≈ 0.7 l_c に相当
    # --- CC の液の混合 (residuals_flatcc.jl の cc_mixing)。:none = 従来。:jet / :plume / :both。SWEEP で上書き可 (scripts/cc_mix_run.jl) ---
    cc_mix = sweep_get(:cc_mix, :none),
    cc_mix_θ = sweep_get(:cc_mix_θ, 1.0),                 # プルームの出口温度の係数 (1 = 持ち上がる高さの上限。感度パラメータ)
    cc_mix_ε_T = sweep_get(:cc_mix_ε_T, 1.0e-3),          # 層の選択の滑らかさの温度幅 [K] (数値の幅。0.5〜2 mK で結果が変わらないことを確認する)
    cc_mix_D_port = sweep_get(:cc_mix_D_port, 1.76e-3),   # 【実測】流入部の管の内径 [m]
    cc_mix_z_port = sweep_get(:cc_mix_z_port, 4.5e-3),    # 【実測】ポートの高さ (secondary wick の内底面から) [m]
    cc_mix_x_jet = sweep_get(:cc_mix_x_jet, 0.025),       # 【仮定】噴流が壁に当たるまでの距離 [m]
    cc_mix_ent = sweep_get(:cc_mix_ent, 1.0),             # 巻き込みの係数 (1 = Schlichting の式そのまま)
    cc_wall_dT = true,             # 壁面凝縮に、駆動する温度差の比 g_T = (T_v − T_swick)/(T_v − T̄_液) をかける (T̄_液 = 液 16 層の体積平均)。false で v8c。壁が液より熱いほど壁面凝縮が減り、T_swick ≥ T_v では 0
    wall_radial_R = false,         # false: 管壁の肉厚方向の熱抵抗は h_con (≈1.65e4 W/m²K、接触 3/4 のとき。SUS304 の肉厚方向コンダクタンス + 周方向の熱の回り込みを含む) に含める。true にすると壁抵抗を別に直列に入れるので、h_con は接触だけの値にする (二重に数えない)
    cc_wall_cond = true,           # 壁面凝縮: 液面の凝縮に比例して、蒸気が接する壁 (A_vap) にも凝縮する。比 = 熱浸透率 √(kρc) × 面積 の比 (液 : secondary wick)。壁面の潜熱は secondary wick へ、凝縮水は層 1 へ。フィット係数なし。false で液面のみ (v8a)
    ne_liq_decouple = false,       # 【数値実験】液側が見る CC の圧力を P_cc − ΔP_ne にする (物理的根拠なし。ΔP_ne の影響の切り分け用)
    δ_layer_cc = 9.0e-3,           # 熱伝導層を解く深さ [m] (= 液深さ h_liq ≈ 9.2 mm。下端は一次ウィック側の壁)
    grv_evap_gate = :nucleation,   # 核生成にヒステリシス (開始 ΔT_onset、維持 f_hold·ΔT_onset)。x ラッチはフラッシュ蒸気で誤作動したので置換
    f_hold_onset = 0.0,            # ラッチ後の維持の閾値 = 0 (EC が飽和温度を下回ると沸騰をやめる)。以前の −0.2 (EC が Tsat より低くても続ける) は根拠が弱いので 0 にして影響を見る
    n_lat_nuc = 0.15,              # n がこれを超えたらラッチ (n → 1)。実験の EC−VL1 は +20 s に 0.34 K (≈ ΔT_onset) → +25 s に 0 へ崩れる (核生成の暴走)。以前の線形のヒステリシスでは部分開度 (n ≈ 0.3) で過熱度が 0.3 K に張り付いた
    ε_n_nuc = 0.03,
    x_hold_nuc = 0.003,            # グルーブの x がこれ以上の間は沸騰を続け、下回ると止める (成層出口の x は 0.005 程度まで。flash は ≈ 2e-3)
    ε_x_nuc = 5.0e-4,
    τ_nuc = sweep_get(:τ_nuc, 1.5),                   # 【数値的な滑らかさのための値。フィッティングパラメータではない】感度は tau_nuc_sweep.jl (0.5 / 1.5 / 3.0 s): 影響は τ = 21〜28 s に限られ、τ = 30 s 以降は HB・EC・RES・ΔP_grv とも差なし (τ = 24 s で EC 0.05 K 以内、ΔP_grv 約 5 kPa 以内)。核生成状態 nuc.n の時定数 [s]。0.5 → 1.5: 実験の EC−VL1 は閾値 (0.339 K) を 0.05 K 超えて (上昇率 0.036 K/s) +17〜19 s にピークをつけ、6〜7 s かけて崩れる (遅れ ≈ 1.4 s)。0.5 s では +16〜17 s の 1 s で ΔP が 6 kPa、LL 末端の速度が 41 mm/s に跳ねた。元の説明:
    cc_boiling = true,             # CC の核沸騰 (プール沸騰。Chen の h をそのまま使う)
    r_c_cc = 35.0e-6,              # CC の活性キャビティ開口半径 [m] (secondary wick の細孔半径 35 µm)
    ε_onset = 0.015,               # 核沸騰開始・界面蒸発の開始の滑らかさ [K]。0.04 では 0.3 K の過熱 (ΔT_onset=0.347 K の手前) で n が 0.2〜0.3 まで開き、+10〜14 s にグルーブの沸騰が始まって LL の流れが実験 (τ≈23 s) より約 10 s 早く始まり、EC−Tsat が 0.3 K に張り付いた。0.1 だと +5〜10 s で蒸発が始まった。硬くなりすぎるなら 0.025
    r_c_wick = r_c_wick_run,   # 【ic_cases.jl: dT_onset を与えたときは換算値】 ウィック界面蒸発のゲートの半径 [m]。:wick = 実測の細孔半径 5.1 µm (ΔT ≈ 0.29 K)。:P = 1.5 µm (ΔT ≈ 0.97 K、調整値。承認済み): 受熱面 P で判定すると、5.1 µm では約 10 s に開いてしまう (実験は約 20 s)
    onset_ref = gate_case === :wick ? :wick : :ec,       # 開始の判定の温度: :wick = primary wick の最外セル、:ec = 受熱面 P。維持は T_ec のまま
    ε_onset_wick = 0.015,          # 界面蒸発のゲートの幅 [K]。0.1 では +0〜10 s に微小な蒸発と押し出し (s = 0.06 m @+10 s) が出て、LL1・LL2・VL1 が動かない実験 (+21 s まで) と合わなかった
    r_c_grv = sweep_get(:r_c_grv, 4.3e-6),   # 【ic_cases.jl で上書き可。C2evap は 1e-8 m (ΔT ≈ 150 K) でグルーブの核沸騰を起こさない】 グルーブの気泡核生成 (蒸発器表面のキャビティ) の開口半径 [m]。実験の EC−VL1 が +20 s に 0.34 K → +25 s に 0 へ崩れる (ΔT_onset ≈ 0.34 K)
    t_init = stage === :heat ? 303.0 : 0.0,
    t_end  = stage === :short ? t_short_end : (stage === :long ? 816.0 : (stage === :long2h ? 6910.0 : (stage === :pre ? 303.0 : 1110.0))),   # :heat は 1110 s まで (加熱開始 +806 s。exp_5W_12min.csv は 1108.8 s まで)
    saveat = stage === :pre ? sort(unique(vcat(collect(0.0:5.0:303.0), collect(20.0:0.5:40.0)))) : (stage in (:short, :long, :long2h) ? sort(unique(round.(vcat((stage === :long2h ? vcat(collect(0.0:1.0:200.0), collect(210.0:10.0:6910.0)) : collect(0.0:1.0:(stage === :long ? 816.0 : t_short_end))), collect(9.0:0.1:15.0), (sweep === nothing ? Float64[] : collect(13.0:0.25:min(t_short_end, 45.0)))), digits = 6))) : 1.0),   # :short は t=9〜15 s を 0.1 s 刻みで保存 (加熱前後の圧力の跳ねの診断)
    init...,                       # T_init, T_cond, T_wall_init_*, init_phase_*, f_l_cc_init
    T_cond = T_cold_plate,         # コールドプレートの温度 [K] (init の T_cond を上書き)
    h_cp = sweep_get(:h_cp, 3000.0),                 # 凝縮器プレートとコールドプレートの接触面の熱伝達係数 [W/m²K] (グリス+ボルト締結)。コールドプレートを2次元の層にしたので、10 mm のアルミの内部の伝導と、チラー配管への広がりは、別に解く (以前の 1000〜3000 は、これらも含む実効値だった)。承認済みの 3000 のまま
    # h_chiller は指定しない (= h_con と同じ。5 W でのチラー配管側の温度降下は 0.1 K 以下で、独立のフィッティング値は不要)
)

println("凝縮器の管: N = ", geo.cond.N, " (1 セル ", round(geo.cond.L / geo.cond.N * 1e3, digits = 2), " mm)。CON1 は入口から 10 mm の位置を補間して読む。チラー側の管 N = ", geo.cplate.N, ",  プレート格子 ", plate_nx_run, "×", plate_ny_run)
println("CC ケース: 壁厚 t_shell = ", geo.cc.t_shell * 1e3, " mm, t_end = ", geo.cc.t_end * 1e3, " mm (外形 W_o = ", geo.cc.W_o * 1e3, " mm),  足 HB→CC ケース G_foot_cc = ", sweep_get(:G_foot_cc, 0.128), " W/K")
println("蒸発器の幾何: A_ec_wick (ケース-ウィック界面、h_ev の面積) = ", round(geo.grv.A_ec_wick * 1e6, digits = 2), " mm²,  A_ec_grv (ケース-溝蒸気の過熱面) = ", round(geo.grv.A_ec_grv * 1e6, digits = 2),
        " mm²,  溝長 L_grv = ", geo.grv.L_grv * 1e3, " mm (ウィック長 L = ", geo.grv.L * 1e3, " mm),  溝の総体積 = ", round(geo.grv.n_grv * geo.grv.A_grv_cross * geo.grv.L_grv * 1e6, digits = 3), " cm³,  G = h_ev A_ec_wick = ", round(condition.h_ev * geo.grv.A_ec_wick, digits = 2), " W/K,  C_ec_P = ", condition.C_ec_P, " J/K")
println("h のブレンドの疑似乾き度: h_xps_mode = ", condition.h_xps_mode, ",  ε_phase_h = ", condition.ε_phase_h, ",  f_h_shah = ", condition.f_h_shah)
u_init = nothing
if stage === :heat
    u_init = Serialization.deserialize(state_file)
    println("加熱段階: ", state_file, " の状態から開始")
end

# IDA の相対許容誤差。エネルギー残差が積分誤差かの確認用 (1e-3 → 1e-4 で残差が減れば積分誤差)。通常は 1e-3
reltol_run = 1e-3
# 1 ステップあたりの Newton 反復の上限 (IDA の既定 4)。蒸気管の蒸気と液の境界 (vl セル 12) で、飽和線をまたぐ前後に、4 回では収束せず、刻みが 1e-6 s まで縮んで止まる。
# 【試行】 10。以前の試験 (solverfuncs.jl) では 8, 15 で効果がなかった (別の状況)。効かなければ 4 に戻す。
newton_max = 10
t_wall = @elapsed (sol, p, u0) = solve_transient_lhp!(geo, condition, tables; u_init = u_init, reltol = reltol_run, max_order = 1, max_nonlinear_iters = newton_max)
println("計算時間: ", round(t_wall, digits = 1), " s (シミュレーション ", condition.t_end - condition.t_init, " s)")
if string(sol.retcode) ∉ ("Success", "Default")
    # IDA が途中で止まった: 最後の保存点の熱伝達 h の内訳を表示 (診断のみ)
    try
        h_failure_diag(sol, p, tables)
        failure_state_diag(sol, p, tables)
        failure_residual_diag(sol, p)
        failure_autopsy(sol, p, geo, condition, u0; tables = tables)   # 診断のみ (数秒)。前線まわりの残差の滑らかさと δ を調べる
    catch err
        @warn "h の診断に失敗" err
    end
end

energy_balance(sol, p; history = true)
try
    save_results(sol, p, tables, result_dir)
catch err
    @warn "結果の保存に失敗 (CSV を開いたままでは?)。表の出力は続けます" err
end

# ループの熱収支と凝縮器の二相長を budget.csv に保存 (診断のみ。plot_heat_budget.jl が読む)
include(joinpath(@__DIR__, "heat_budget_calc.jl"))
try
    write_heat_budget(sol, p, tables, result_dir; t_heat = t_heat)
catch err
    @warn "熱収支の保存に失敗" err
end

K = 273.15
cellpos(N) = [(i - 0.5) / N for i in 1:N]
at(Tw, pos) = pchip_interp(cellpos(length(Tw)), Tw, pos) - K
# CON1 の位置 (管の入口から 10 mm、ご提示 11/2) を、管の長さに対する比で表す。以前は at(.., 0.0) = 入口側の最初のセルの中心 (N=50 で 5.6 mm)。N を変えても同じ位置を読む
con1_pos = 0.010 / geo.cond.L
# 管壁の熱抵抗を入れた外面温度 [K] (出力用の後処理。モデルの計算は変えない。LHPTransientModel.cond_outer_wall_T)。プレートが無ければ壁の節点の温度
cond_out(ck) = something(LHPTransientModel.cond_outer_wall_T(ck, p), ck.T_wall.cond)

if stage === :pre
    string(sol.retcode) in ("Success", "Default") || error("IDA が t = 303 s まで到達しませんでした (最終時刻 " * string(sol.t[end]) * " s)。state_t303.jls は保存しません")
    mkpath(result_dir)
    Serialization.serialize(state_file, collect(sol.u[end]))
    println("t = 303 s の状態を保存しました: ", state_file, "  (stage = :heat で続きを実行)")
    c = components_at(sol, p, length(sol.t))

    # --- 加熱前の計算結果 (t = 303 s) と、実験 (t = 290〜302 s の平均) の比較 ---
    println("\n加熱前 (t = 303 s) の壁温度 [°C]:   センサ | 実験 | 計算 | 差")
    rows = [
        ("HB",   sensors_ref.HB,  c.T_wall.hb[1] - K),
        ("EC",   (sensors_ref.EC1 + sensors_ref.EC2) / 2, (condition.ec_2node ? c.eca.T[1] : c.T_wall.evap[1]) - K),
        ("RES",  (sensors_ref.RES1 + sensors_ref.RES2) / 2, c.T_wall.cc[1] - K),
        ("VL1",  sensors_ref.VL1, at(c.T_wall.vl, 1/6)), ("VL2", sensors_ref.VL2, at(c.T_wall.vl, 1/2)),
        ("VL3",  sensors_ref.VL3, at(c.T_wall.vl, 5/6)),
        ("CON1", sensors_ref.CON1, at(c.T_wall.cond, con1_pos)), ("CON2", sensors_ref.CON2, at(c.T_wall.cond, 0.25)),
        ("CON3", sensors_ref.CON3, at(c.T_wall.cond, 0.5)), ("CON4", sensors_ref.CON4, at(c.T_wall.cond, 0.75)),
        ("CON5", sensors_ref.CON5, at(c.T_wall.cond, 1.0)),
        ("CON1外", sensors_ref.CON1, at(cond_out(c), con1_pos)), ("CON2外", sensors_ref.CON2, at(cond_out(c), 0.25)),   # 外面 (管壁の抵抗を入れた推定)
        ("LL1",  sensors_ref.LL1, at(c.T_wall.ll, s_LL1)), ("LL2", sensors_ref.LL2, at(c.T_wall.ll, s_LL2)),
        ("LL3",  sensors_ref.LL3, at(c.T_wall.ll, s_LL3)),
    ]
    for (name, Te, Tm) in rows
        println("  ", rpad(name, 5), " ", round(Te, digits = 3), "  ", round(Tm, digits = 3), "  ", round(Tm - Te, digits = 3))
    end
    println("CC の圧力 = ", round(c.cc.P[1] / 1e3, digits = 2), " kPa,  T_cc = ", round(c.cc.T[1] - K, digits = 3), " °C")
else
    # --- 加熱開始 (t_heat) からの経過時間ごとの壁温度 [°C] ---
    # 実験 (グラフからの読み取り、5 W): 起動は加熱開始から約 48 s。起動後 (約 +50 s) は
    #   HB ≈ 24.4、EVA1 ≈ 23.55 (起動時に 23.65 まで上がって下がる)、RES1 ≈ 23.4、VL1 ≈ 23.4 (一時 23.65)、
    #   CON1 ≈ 22.4 (起動時に 21.8 へ跳ねる)、LL3 ≈ 22.5 (加熱開始 +20 s に 22.9 まで上がる)、CON2, CON3 ≈ 20.4
    println("\n加熱開始からの時間 [s], 壁温度 [°C]:  HB | EC | RES | VL1 | CON1 | CON2 | LL3 |  P_cc [kPa] | Q_apply, Q_hb_ec [W]")
    for τ in (0, 10, 20, 30, 40, 48, 60, 75, 90)
        t = t_heat + τ
        t > condition.t_end && continue
        k = argmin(abs.(sol.t .- t))
        ck = components_at(sol, p, k)
        Qhb = condition.h_hb_ec * condition.A_hb_ec * (ck.T_wall.hb[1] - ck.T_wall.evap[1])
        println("  +", lpad(τ, 3), "  ", round(ck.T_wall.hb[1] - K, digits = 2), "  ", round(ck.T_wall.evap[1] - K, digits = 2),
                "  ", round(ck.T_wall.cc[1] - K, digits = 2), "  ", round(at(ck.T_wall.vl, 1/6), digits = 2),
                "  ", round(at(ck.T_wall.cond, con1_pos), digits = 2), "  ", round(at(ck.T_wall.cond, 0.25), digits = 2),
                "  ", round(at(ck.T_wall.ll, s_LL3), digits = 2), "  ", round(ck.cc.P[1] / 1e3, digits = 1),
                "  ", round(applied_heat(sol.t[k], condition), digits = 2), "  ", round(Qhb, digits = 2))
    end

    # --- 核沸騰の判定の診断: 蒸発の開始が過熱度 ΔT_onset で起きているか、x の維持が働いているか ---
    #   EC−Tsat は EC の壁温度とグルーブ圧力の飽和温度の差 (これが ΔT_onset ≈ 0.37 K を超えたら開始)。
    if condition.grv_evap_gate === :nucleation
        println("\n核沸騰の診断 (加熱開始からの時間 τ):  P−Tsat(P_grv) [K] (ゲートに使う受熱面の節点) (S: 側面 = センサ側、W: primary wick 最外セル = 開始の判定) | n (セル平均, 最大) | x_grv (最大) | ΔT_onset (界面蒸発のゲート) [K] | f_l (平均, 出口セル) | h_ev A [W/K]")
        for τ in (0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 90)
            t = t_heat + τ
            t > condition.t_end && continue
            k = argmin(abs.(sol.t .- t))
            ck = components_at(sol, p, k)
            Ts = LHPTransientModel.Tsat(ck.grv.P[1], tables)
            spk = LHPTransientModel.sat_props(Ts, tables)
            dTon = 2 * spk.σ * Ts / (spk.ρ_v * spk.h_fg * (condition.r_c_wick > 0 ? condition.r_c_wick : condition.r_c_grv))
            println("  +", lpad(τ, 3), "   ", lpad(round(ck.T_wall.evap[1] - Ts, digits = 3), 7), " (S:", lpad(round((condition.ec_2node ? ck.eca.T[1] : ck.T_wall.evap[1]) - Ts, digits = 3), 7), " W:", lpad(round(ck.wick.T[end] - Ts, digits = 3), 7), ")",
                    "   ", lpad(round(sum(ck.nuc.n) / length(ck.nuc.n), digits = 3), 6), " ", lpad(round(maximum(ck.nuc.n), digits = 3), 6), " nc=", lpad(round(maximum(ck.nuc.nc), digits = 3), 5),
                    "   ", lpad(round(maximum(ck.grv.x), digits = 5), 8), "   ", round(dTon, digits = 3),
                    "   ", lpad(round(sum((1 - ck.grv.x[i]) * ck.grv.ρm[i] / tables.ρ_l_sat(ck.grv.T[i]) for i in eachindex(ck.grv.x)) / length(ck.grv.x), digits = 3), 6),
                    " ", lpad(round((1 - ck.grv.x[end]) * ck.grv.ρm[end] / tables.ρ_l_sat(ck.grv.T[end]), digits = 3), 6),
                    "   ", round(LHPTransientModel.h_ev_at(condition, t) * geo.grv.A_ec_wick, digits = 1))
        end
    end

    # --- 熱の行き先: 加熱開始から各時刻までに、印加熱がどの部位の蓄積エネルギーになったか ---
    # 実験の EC は +20 s で +0.37 K (モデル +0.61 K)。差の約 2 J がどこへ入っているかを見る
    # (流体の欄は顕熱と相変化の潜熱の合計。ウィック・グルーブ・CC・蒸気管の液が大きいかどうかを比べる)。
    energy_partition(sol, p; t0 = t_heat, ts = [t_heat + τ for τ in (10, 15, 20, 25, 30, 35, 40, 50) if t_heat + τ <= condition.t_end], tables = tables)
    ll_end_budget(sol, p, tables; ts = [t_heat + τ for τ in (-9, -5, 0, 5, 10, 14, 16, 18, 20, 25) if condition.t_init <= t_heat + τ <= condition.t_end])

    # --- 実験 (5 W、加熱開始 303.8 s) との時系列の比較: 全センサ ---
    # data/ の実験の csv (exp_5W_startup.csv, exp_5W_12min.csv, exp_5W_2h.csv) の値を、加熱開始からの時間 τ で線形補間して、モデルと比べる。
    exp_file = joinpath(@__DIR__, "..", "data", stage === :long2h ? "exp_5W_2h.csv" : (stage in (:heat, :long) ? "exp_5W_12min.csv" : "exp_5W_startup.csv"))
    t_exp_onset = 303.8
    if isfile(exp_file)
        lines = filter(!isempty, readlines(exp_file))
        hdr = [strip(replace(String(x), "\ufeff" => "", "\"" => "")) for x in split(lines[1], ",")]   # BOM・引用符・空白を除く
        vals = [(y = tryparse(Float64, strip(x)); y === nothing ? NaN : y) for l in lines[2:end] for x in split(l, ",")]   # 空欄 (測定なし) は NaN
        Dm = permutedims(reshape(vals, length(hdr), :))
        function ecol(nm)
            # 列名は "HB" でも、記録装置の生の列名 ("101_HB"、"116_RES_2" など) でもよい: 先頭の番号と "_" を除いて比べる
            normh(x) = replace(replace(x, r"^\d+_" => ""), "_" => "")
            j = findfirst(==(nm), hdr)
            j === nothing && (j = findfirst(x -> normh(x) == normh(nm), hdr))
            j === nothing && (println("警告: 実験 csv に列 \"", nm, "\" がありません。列名: ", hdr); return fill(NaN, size(Dm, 1)))
            return Dm[:, j]
        end
        function lin1(xs, ys, x)
            x <= xs[1] && return ys[1]
            x >= xs[end] && return ys[end]
            kk = searchsortedlast(xs, x)
            return ys[kk] + (ys[kk+1] - ys[kk]) * (x - xs[kk]) / (xs[kk+1] - xs[kk])
        end
        te = ecol("t")
        sens = [
            ("HB",   ["HB"],           ck -> ck.T_wall.hb[1] - K),
            ("EC",   ["EC1", "EC2"],   ck -> (condition.ec_2node ? ck.eca.T[1] : ck.T_wall.evap[1]) - K),
            ("RES",  ["RES1", "RES2"], ck -> ck.T_wall.cc[1] - K),
            ("VL1",  ["VL1"],  ck -> at(ck.T_wall.vl, 1/6)), ("VL2", ["VL2"], ck -> at(ck.T_wall.vl, 1/2)),
            ("VL3",  ["VL3"],  ck -> at(ck.T_wall.vl, 5/6)),
            ("CON1", ["CON1"], ck -> at(ck.T_wall.cond, con1_pos)),  ("CON2", ["CON2"], ck -> at(ck.T_wall.cond, 0.25)),
            ("CON3", ["CON3"], ck -> at(ck.T_wall.cond, 0.5)),  ("CON4", ["CON4"], ck -> at(ck.T_wall.cond, 0.75)),
            ("CON5", ["CON5"], ck -> at(ck.T_wall.cond, 1.0)),
            ("CON1外", ["CON1"], ck -> at(cond_out(ck), con1_pos)), ("CON2外", ["CON2"], ck -> at(cond_out(ck), 0.25)),
            ("CON3外", ["CON3"], ck -> at(cond_out(ck), 0.5)),      ("CON4外", ["CON4"], ck -> at(cond_out(ck), 0.75)),
            ("CON5外", ["CON5"], ck -> at(cond_out(ck), 1.0)),
            ("LL1",  ["LL1"],  ck -> at(ck.T_wall.ll, s_LL1)),  ("LL2", ["LL2"], ck -> at(ck.T_wall.ll, s_LL2)),
            # LL3 は比較から外す (熱電対の値が流体の温度を表していないと判断。LL2 で比較する)
        ]
        τs = [τ for τ in (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 75, 90, 120, 180, 240, 300, 400, 500, 600, 700, 800, 1200, 1800, 2400, 3600, 4800, 6000, 6800) if t_heat + τ <= condition.t_end]
        cks = [components_at(sol, p, argmin(abs.(sol.t .- (t_heat + τ)))) for τ in τs]
        println("\n実験との比較 (加熱開始からの時間 τ [s]):  各センサで モデル / 実験 / 差 (モデル − 実験) [°C or K]")
        println(rpad("", 9), join([lpad("+" * string(τ), 7) for τ in τs]))
        for (nm, cols, f) in sens
            mo = [f(ck) for ck in cks]
            ex = [sum(lin1(te, ecol(cn), t_exp_onset + τ) for cn in cols) / length(cols) for τ in τs]
            println(rpad(nm * " 計算", 9), join([lpad(string(round(v, digits = 2)), 7) for v in mo]))
            println(rpad(nm * " 実験", 9), join([lpad(string(round(v, digits = 2)), 7) for v in ex]))
            println(rpad(nm * " 差",   9), join([lpad(string(round(v, digits = 2)), 7) for v in (mo .- ex)]))
        end

        # --- CC 気液界面の診断 (ΔP_ne の減衰の原因を調べる。診断のみ) ---
        #   VL1−RES (加熱前からの変化) は圧力を使わない ΔP_ne の目安 (実験は +26〜42 s に 0.4 K で平坦、+58 s に 0)
        dT_model(ck) = at(ck.T_wall.vl, 1/6) - (ck.T_wall.cc[1] - K)
        dT_exp(τ) = (lin1(te, ecol("VL1"), t_exp_onset + τ) - lin1(te, ecol("VL1"), t_exp_onset)) -
                    ((lin1(te, ecol("RES1"), t_exp_onset + τ) + lin1(te, ecol("RES2"), t_exp_onset + τ)) / 2 -
                     (lin1(te, ecol("RES1"), t_exp_onset) + lin1(te, ecol("RES2"), t_exp_onset)) / 2)
        cc_interface_diag(sol, p, tables; t0 = t_heat, ts = collect(10:2:70), dT_vl_res = dT_model, dT_vl_res_exp = dT_exp)

        # --- 追加の診断 (診断のみ): 起動の立ち上がり、液管の温度の前線、前線通過後のループの熱収支 ---
        startup_onset_diag(sol, p, tables; t0 = t_heat, ts = collect(10:1:28))
        ll_thermal_diag(sol, p, tables; t0 = t_heat, ts = collect(16:2:60), s_LL1 = s_LL1, s_LL2 = s_LL2)
        loop_heat_diag(sol, p, tables; t0 = t_heat, ts = [10, 20, 30, 40, 44, 48, 50, 52, 56, 60, 66, 75, 90, 120, 180, 240, 300, 400, 500, 600, 700, 800, 1200, 1800, 3600, 6000] |> v -> filter(<=(condition.t_end - t_heat), v))

        # --- 液管 (LL) の押し出し量と前線の到達時刻 ---
        #   s(τ) = ∫u dt: 液管の中央のセルの速度を積算した、液が進んだ距離 [m] (液管の全長 1.15 m = 約 0.9 mL)。
        #   前線 = 凝縮器出口の液が、センサ位置に着く時刻 (s がセンサの位置 [凝縮器から] に達する τ)。
        #   壁温 0.1 K 低下 = 加熱前 (τ = 0) から 0.1 K 下がる時刻 (計算・実験)。実験の LL3 が下がらないなら NaN。
        ks   = [k for k in eachindex(sol.t) if sol.t[k] >= t_heat]
        tk   = sol.t[ks]
        cks2 = [components_at(sol, p, k) for k in ks]
        uk   = [ck.ll.u[10] for ck in cks2]
        sk   = vcat(0.0, cumsum((uk[2:end] .+ uk[1:end-1]) ./ 2 .* diff(tk)))
        L_ll = 1.15
        println("\n液管の押し出し量 s(τ) [m] (全長 ", L_ll, " m):  ",
                join(["τ=" * string(τ) * ": " * string(round(lin1(tk, sk, t_heat + τ), digits = 2)) for τ in (10, 20, 30, 40, 50, 60, 90) if t_heat + τ <= condition.t_end], "  "))
        println("液管のセンサ位置への前線の到達 / 壁温 0.1 K 低下の時刻 τ [s]:  前線(計算) | 壁温(計算) | 壁温(実験)")
        for (nm, xs, sf) in (("LL1", 0.10, s_LL1), ("LL2", L_ll - 0.50, s_LL2))
            i_f = findfirst(>=(xs), sk)
            τ_front = i_f === nothing ? NaN : tk[i_f] - t_heat
            Twk = [at(ck.T_wall.ll, sf) for ck in cks2]
            i_m = findfirst(<(Twk[1] - 0.1), Twk)
            τ_mod = i_m === nothing ? NaN : tk[i_m] - t_heat
            y = ecol(nm); y0 = lin1(te, y, t_exp_onset)
            i_e = findfirst(kk -> te[kk] >= t_exp_onset && y[kk] < y0 - 0.1, eachindex(te))
            τ_exp = i_e === nothing ? NaN : te[i_e] - t_exp_onset
            println("  ", nm, "   ", lpad(round(τ_front, digits = 1), 6), "  ", lpad(round(τ_mod, digits = 1), 6), "  ", lpad(round(τ_exp, digits = 1), 6))
        end
    else
        println("\n(実験との比較: ", exp_file, " がありません。実験の csv を data/ に置いてください)")
    end
end