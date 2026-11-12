# =============================================================================
# residuals_grv.jl (グルーブ)
#
# グルーブは、ウィックの蒸発界面から蒸気が湧き出し、蒸発器ケースから
# 過熱され、蒸気管へ流れ出す、軸方向の流路。n_grv 本の溝を、全溝の断面積の
# 合計を持つ1本の管として扱い、輸送管の核 (residual_transport!) を、
# 質量・エネルギーのソースと、蒸発器ケースからの熱を与えて呼ぶ。
#
#   質量:       ∂ρm/∂t + ∂(ρm u)/∂z = S_m                       (S_m: ウィックの蒸発)
#   エネルギー: ρ ∂Em/∂t + ρ u ∂φ/∂z − (P/ρ)∂ρm/∂t = q + S_E − φ S_m
#                q   : 蒸発器ケースからグルーブ蒸気への熱 (式45)
#                S_E : 蒸発界面を通過するエネルギー (wick_groove_source)
#
# 【幾何】 A = n_grv A_grv_cross (全溝の合計)、水力直径 D_h = 2wh/(w+h)。
# 溝の断面は矩形 (幅 w、深さ h)。矩形ダクトの層流摩擦係数は、円管の 64/Re
# ではなく、アスペクト比に依存する (正方形で約 57/Re)。ここでは円管の式を、
# 水力直径で使っている (誤差は約12%)。グルーブの摩擦圧損は小さい
# (蒸気 ṁ=1.7e-5 kg/s で、溝1本の全長 35 mm で約 2 Pa) ので、影響は軽微。
#
# 【境界】 入口 (溝の端) は閉端 (u=0)。ゴーストは、セル1と同じ値にして、
# 入口側の勾配をゼロにする。出口は蒸気管の最初のセル
# (transport_ghost_down(u_vl, p_vl, p_grv) で作る)。
#
# 【蒸発器ケースからグルーブ蒸気への熱 Q_ec-grv (式45)】
#   Q_ec-grv,i = h_ec-grv,i (A_ec-grv / N) (T_ec − T_i)
# 過熱面は蒸発器ケースとグルーブ蒸気が接する上面だけ (A_ec_grv = width L n_grv)。
#
#   h_ec-grv は、ドキュメントの式72 (水平平板の強制対流) に従う:
#     Re_L = G L / μ_v,   G = (ṁ_evap/n_grv)/A_grv_cross
#     Nu_L = 0.664 Re_L^(1/2) Pr^(1/3)    (Re_L < 5e5)
#          = 0.0296 Re_L^(4/5) Pr^(1/3)   (それ以外)
#     h = Nu_L k_v / L                    (L: 溝の長さ。Re と h の代表長さを L に統一)
#   ṁ_evap はウィックの蒸発量 (全体)。溝1本あたりに直すため n_grv で割る。
#   物性 (μ_v, k_v, Pr) はグルーブのセルの状態から取る (ρ は Re の中で相殺)。
#
# 【旧コードとの違い】
#   - 旧 calc_h_ec_gr は、Re を溝の深さ (height) 基準、h を Nu k/L (溝の長さ L で
#     割る) としており、代表長さが混在していた (h が約1/6 になる)。
#     本実装は、ドキュメントに従い、Re も h も溝の長さ L に統一する。
#   - 旧は Re≤2300 で層流/乱流を切り替え、乱流側は 0.0371 Pr Re^0.8 だった。
#     本実装は式72 (遷移 Re=5e5、乱流側 0.0296 Re^0.8 Pr^(1/3))。
#   - 旧 steadyfuncs.jl は溝1本あたりの流量、過渡側は 12本まとめて計算していた。
#     本実装は、全体の ṁ_evap を n_grv で割って溝1本あたりにする。
#
# 【判断事項】
#   1. ウィックの蒸発界面の圧力 P_grv は、グルーブの圧力の平均
#      (groove_pressure)。ウィックは半径方向の1次元で、軸方向の位置を
#      区別しないため。
# =============================================================================

# -----------------------------------------------------------------------------
# 1. 幾何
# -----------------------------------------------------------------------------

"""
    grv_params(geo::Union{GrvGeometry,FlatGrvGeometry}) -> NamedTuple

グルーブの幾何量。pipe_params と同じフィールド (N, dl, D, A, V, θ, elev0) を
持つので、residual_transport! や transport_ghost_up/down にそのまま渡せる。
A, V は全溝の合計。D は矩形溝の水力直径。
"""
grv_params(geo::Union{GrvGeometry,FlatGrvGeometry}) = (
    N = geo.N, dl = geo.dl, L = groove_length(geo),
    D = 2 * geo.width * geo.height / (geo.width + geo.height),
    A = geo.n_grv * geo.A_grv_cross,
    V = geo.n_grv * geo.A_grv_cross * geo.dl,
    θ = geo.θ, elev0 = geo.elev0,
    n_grv = geo.n_grv, A_one = geo.A_grv_cross, A_ec_grv = geo.A_ec_grv,
    width = geo.width, height = geo.height,
)

"""
    groove_pressure(u_grv) -> P_grv

ウィックの蒸発界面温度 T_sat(P_grv) に使う、グルーブの代表圧力
(セルの圧力の平均)。N=1 なら、そのセルの圧力。
"""
groove_pressure(u_grv) = sum(u_grv.P) / length(u_grv.P)

# -----------------------------------------------------------------------------
# 2. 蒸発器ケースからグルーブ蒸気への熱
# -----------------------------------------------------------------------------

"""
    groove_h_ec(u, p, ṁ_evap, g, tables; δ_m=1e-9, ε_Re=2e4) -> h [W/m²K] (セルごと)

蒸発器ケース→グルーブ蒸気の熱伝達係数 (式72、水平平板の強制対流)。
  Re_L = G L/μ_v,  G = (|ṁ_evap|/n_grv)/A_grv_cross
  h = Nu_L k_v / L   (Nu_L は Nu_flatplate_smooth: 遷移 Re=5e5)
μ_v, k_v, Pr はグルーブのセルの状態 (ρm, Em, u) から transport_props で取る。
ṁ_evap の符号は soft_abs で処理する (δ_m は ṁ の典型値 1e-5 kg/s より十分小さく)。
入口は閉端 (u=0)。
"""
function groove_h_ec(u, p, ṁ_evap, g, tables::FluidTables; δ_m=1.0e-9, ε_Re=2.0e4)
    h = Vector{Float64}(undef, p.N)
    G = soft_abs(ṁ_evap, δ_m) / (p.n_grv * p.A_one)
    for i in 1:p.N
        u_c = 0.5 * ((i == 1 ? 0.0 : u.u[i-1]) + u.u[i])
        em = internal_energy(u.Em[i], u_c, cell_elevation(p, i), g)
        tp = transport_props(u.ρm[i], em, tables)
        Re = G * p.L / tp.μ
        h[i] = h_ec_grv(Re, tp.Pr, p.L, tp.k; ε_Re = ε_Re)
    end
    return h
end

"""
    groove_evap_heat(u, T_ec, h, p) -> Q [W] (セルごと)

蒸発器ケースからグルーブの各セルの蒸気へ渡す熱 (式45)。蒸発器ケースが
熱ければ正。蒸発器ケースの残差 (residual_evapcase) では sum(Q) を失う熱と
して、グルーブの残差 (residual_grv!) ではこの値を熱源として使う
(同じ値を共有するので、エネルギーが厳密に一致する)。
"""
groove_evap_heat(u, T_ec, h, p) = [h[i] * (p.A_ec_grv / p.N) * (T_ec - u.T[i]) for i in 1:p.N]

# -----------------------------------------------------------------------------
# 2b. グルーブの濡れ (液の有無) と核沸騰 (ドキュメント 2.8.1)
# -----------------------------------------------------------------------------
#
# 質量の扱い: グルーブは均質二相 (ρm, Em) の流路なので、グルーブ内の液が沸騰して蒸気に
# なる変化は、壁からの熱 Q_ec-grv が Em を増やし、x が上がるだけで表現される。
# 質量ソースは足さない (ドキュメント式86の Q_evap は熱の内訳で、質量源ではない)。
# グルーブの質量が変わるのは、ウィック界面の蒸発 ṁ (S_m) と、蒸気管への流出だけ。
# グルーブが濡れている間は、界面の ṁ を 0 にする (wick_interface の w_wet)。

"""
    groove_wet_weight(u, p, g, tables; f_th=1e-3, ε_wet=5e-4) -> Vector (セルごと, 0〜1)

グルーブのセルが濡れている (液がある) 重み。液の体積分率
    f_l = (1 − x) ρm / ρ_l(T)
を使う:  w_wet = smooth_step(f_l − f_th, ε_wet)。液が残っていれば 1、液がなくなれば 0。
(ドキュメントの「x_grv = 1 で起動判定」に対応する。)

x ではなく f_l を使う理由: 乾いたと判定されると、界面から蒸気が入る (ṁ>0)。蒸気が入ると
x が上がるので、x で判定すると、x↑ → 乾き → ṁ↑ → x↑ の正のフィードバックになり、
判定の境目で不安定になる。f_l (液の量) は、蒸気を足しても変わらないので、この
フィードバックがない。
"""
function groove_wet_weight(u, p, g, tables::FluidTables; f_th=1.0e-3, ε_wet=5.0e-4)
    w = Vector{Float64}(undef, p.N)
    for i in 1:p.N
        f_l = (1 - u.x[i]) * u.ρm[i] / tables.ρ_l_sat(u.T[i])
        w[i] = smooth_step(f_l - f_th, ε_wet)
    end
    return w
end

"""
    groove_subcool_weight(u, p, tables; ΔT_th=0.1, ε=0.1) -> Vector (セルごと, 0〜1)

グルーブのセルがサブクール (飽和温度に達していない) の重み:
    w_sub = smooth_step(T_sat(P) − T − ΔT_th, ε)
サブクール液のあいだは 1 (界面蒸発を止める)、飽和に達して沸騰が始まれば 0 (蒸発を許す)。
蒸気が入ると T→T_sat になるので、一度開いたゲートは閉じない (正のフィードバックにならない)。
"""
function groove_subcool_weight(u, p, tables::FluidTables; ΔT_th=0.1, ε=0.1)
    w = Vector{Float64}(undef, p.N)
    for i in 1:p.N
        w[i] = smooth_step(Tsat(u.P[i], tables) - u.T[i] - ΔT_th, ε)
    end
    return w
end

"""
    evap_stop_weight(u, p, c, tables, w_wet_c) -> (w_stop, w_sub_c)

界面蒸発を止める重み (wick_interface の w_wet に渡す)。
  :dry        … 液がある間 (w_wet の平均)
  :saturation … 「サブクールの液がある」間だけ = mean(w_sub × w_liq)。w_liq はほぼ液だけ (x が小さい)
                のセルで 1。蒸気が増えたセルは、サブクール判定に関係なく蒸発を許す。
"""
function evap_stop_weight(u, p, c, tables::FluidTables, w_wet_c; T_ec = nothing, n_nuc = nothing)
    if c.grv_boiling && c.grv_evap_gate === :nucleation
        # 沸騰開始のヒステリシスつきの判定 (状態変数 nuc.n。グルーブの沸騰の h と同じ n)。
        # n = 0 の間は界面蒸発を止め、n = 1 (沸騰中) で開く。x のラッチは使わない
        # (加熱直後の flash の蒸気 x ≈ 2e-3 でラッチが開いてしまうため)。
        n_nuc === nothing && error("evap_stop_weight: grv_evap_gate = :nucleation には n_nuc (状態 nuc.n) が必要")
        w_on = [min(max(n_nuc[i], 0.0), 1.0) for i in 1:p.N]
        return sum(1 .- w_on) / p.N, 1 .- w_on
    end
    if c.grv_boiling && c.grv_evap_gate === :onset
        # グルーブの核沸騰と同じ判定 (nucleation_onset_weight)。液だけのセル (x が小さい) では、
        # T_ec が T_sat + ΔT_onset を超えるまで界面蒸発を止める。蒸気が増えたセル (x > x_liq_gate) では
        # 開いたままにする (沸騰開始のヒステリシス: 一度沸騰が始まれば、過熱度が下がっても続く)。
        T_ec === nothing && error("evap_stop_weight: grv_evap_gate = :onset には T_ec が必要")
        w_on = [nucleation_onset_weight(T_ec, Tsat(u.P[i], tables), sat_props(Tsat(u.P[i], tables), tables),
                                        c.r_c_grv, c.ε_onset) for i in 1:p.N]
        w_liq = [1 - smooth_step(u.x[i] - c.x_latch_onset, c.ε_x_latch) for i in 1:p.N]
        w = w_liq .* (1 .- w_on)
        return sum(w) / length(w), 1 .- w_on
    end
    if c.grv_boiling && c.grv_evap_gate === :saturation
        w_sub = groove_subcool_weight(u, p, tables; ΔT_th = c.ΔT_sat_th, ε = c.ε_sat)
        # サブクールの判定は、ほぼ液だけのセル (x が小さい) にだけ使う。蒸気が増えたセルでは、
        # 状態表の補間誤差などで T が Tsat(P) より 0.1〜0.2 K 低く出ることがあり、ゲートが
        # 半分閉じて ṁ が揺れた (200 s の計算の 196〜200 s)。x が大きくなる側へゲートは開くだけで、
        # x↑ → ゲート閉 の正のフィードバックにはならない。
        w_liq = [1 - smooth_step(u.x[i] - c.x_liq_gate, c.ε_x_gate) for i in 1:p.N]
        w = w_sub .* w_liq
        return sum(w) / length(w), w_sub
    end
    return sum(w_wet_c) / length(w_wet_c), zeros(p.N)
end

"""
    nucleation_onset_weight(T_ec, Ts, sp, r_c, ε_onset) -> w_on (0〜1)

核沸騰開始の重み: w_on = smooth_step(T_ec − T_sat − ΔT_onset, ε_onset)、
ΔT_onset = 2σ T_sat/(ρ_v λ r_c) (式82)。グルーブの沸騰の h (groove_h_ec_boiling) と、
界面蒸発のゲート (evap_stop_weight, grv_evap_gate = :onset) の両方がこの関数を使う
(同じ判定にするため)。
"""
function nucleation_onset_weight(T_ec, Ts, sp, r_c, ε_onset)
    ΔT_onset = 2 * sp.σ * Ts / (sp.ρ_v * sp.h_fg * r_c)
    return smooth_step(T_ec - Ts - ΔT_onset, ε_onset)
end

"""
    nucleation_target(T_ec, Ts, sp, r_c, ε_onset, f_hold, n) -> n_target (0〜1)

核沸騰の活性度 n の目標値 (ラッチつき):
    n_target = smooth_step(T_ec − T_sat − ΔT_eff, ε_onset)
    ΔT_eff = ΔT_hold + (ΔT_onset − ΔT_hold)(1 − s_n),   s_n = smooth_step(n − n_lat, ε_n),   ΔT_hold = f_hold ΔT_onset
n < n_lat (沸騰していない) では開始の閾値 ΔT_onset、n > n_lat (沸騰中) では維持の閾値 ΔT_hold。
n は [0, 1] に丸めて使う。
"""
function nucleation_target(T_ec, Ts, sp, r_c, ε_onset, f_hold, n; x = 0.0, x_hold = 2.0, ε_x = 1.0e-3,
    n_lat = 0.15, ε_n = 0.03, T_onset = nothing)
    ΔT_onset = 2 * sp.σ * Ts / (sp.ρ_v * sp.h_fg * r_c)
    n_c = min(max(n, 0.0), 1.0)
    # ラッチ (双安定): n が n_lat を超えたら、閾値を一気に維持側 ΔT_hold = f_hold ΔT_onset へ切り替える。
    # (線形のヒステリシスだと、n = Q_ev/(G ΔT) の部分開度で過熱度が ΔT_onset に張り付くので、ラッチにした。
    #  いったん n が n_lat を超えると n → 1 へ進み、過熱度が崩れる = 気泡の核生成の暴走に相当。)
    s_n = smooth_step(n_c - n_lat, ε_n)
    ΔT_eff = f_hold * ΔT_onset + (1 - f_hold) * ΔT_onset * (1 - s_n)
    # T_onset を与えたとき (onset_ref = :wick): 開始の判定 (n < n_lat) はその温度 (ウィック最外セル) で、
    # 維持の判定 (n > n_lat) は従来どおり T_ec で行う。蒸発が始まるとウィック最外セルは飽和温度以下に下がる
    # (冷たい液が上がってくる) ので、維持までウィックの温度で判定すると、開いた直後に閉じてしまう。
    T_ref = T_onset === nothing ? T_ec : s_n * T_ec + (1 - s_n) * T_onset
    w_T = smooth_step(T_ref - Ts - ΔT_eff, ε_onset)
    # 維持: グルーブに蒸気 (気泡) が x_hold 以上ある間は、過熱度に関係なく沸騰を続ける。
    # x が x_hold を下回ると維持が外れ、再び過熱度 (開始の閾値) だけで決まる。
    # x_hold = 2 (既定) なら w_x = 0 で、この項は効かない。
    # 維持の項に n_c を掛ける: すでに沸騰している (n ≈ 1) セルだけを維持し、開始はしない。
    # (掛けないと、加熱直後の flash の蒸気 x ≈ 2〜3e-3 だけで n が立ち上がり、過熱度 ΔT_onset に
    #  届く前に蒸発が始まってしまった。)
    w_x = smooth_step(x - x_hold, ε_x) * n_c
    return w_T + w_x - w_T * w_x          # なめらかな OR
end

"""
    groove_h_ec_boiling(u, p, ṁ_evap, T_ec, w_wet, g, tables; r_c, ε_onset=0.1) -> h [W/m²K] (セルごと)

濡れ (w_wet) を考えた、蒸発器ケース→グルーブの熱伝達係数:

    h = w_wet h_wet + (1 − w_wet) h_dry
    h_dry  = groove_h_ec (蒸気側の強制対流、式81。従来と同じ)
    h_wet  = (1 − w_on) h_nc + w_on h_chen
    w_h    = w_wet · smooth_step(f_l − f_th_h, ε_h)   (液膜が薄くなったら h_dry へ。h = w_h h_wet + (1−w_h) h_dry)
    h_nc   : 液への自然対流 (式83)。Nu_cc_smooth(Ra) k_l/H (H: 溝の深さ。Ra が小さいので伝導に近い)
    h_chen : Chen の式 (h_chen_smooth)。サブクール沸騰・飽和沸騰とも同じ h を使う (式85, 96)
    w_on   = smooth_step(T_ec − T_sat(P) − ΔT_onset, ε_onset)
    ΔT_onset = 2σ T_sat/(ρ_v λ r_c)    (式82、r_c: 活性キャビティ開口半径)

どれも、状態 (ρm, Em, P, T, x, u) と T_ec から決まる (if なし)。
"""
function groove_h_ec_boiling(u, p, ṁ_evap, T_ec, w_wet, g, tables::FluidTables;
    r_c, ε_onset=0.1, f_th_h=2.0e-2, ε_h=1.0e-2, δ_m=1.0e-9, ε_Re=2.0e4, n_nuc=nothing)

    # n_nuc (状態 nuc.n、セルごと) を与えると、沸騰開始の重み w_on にそれを使う (ヒステリシスつき)
    h_dry = groove_h_ec(u, p, ṁ_evap, g, tables; δ_m = δ_m, ε_Re = ε_Re)
    h = Vector{Float64}(undef, p.N)
    for i in 1:p.N
        P  = u.P[i]
        Ts = Tsat(P, tables)
        sp = sat_props(Ts, tables)

        # 液への自然対流 (溝の深さ H で Ra)
        ΔT  = soft_abs(T_ec - u.T[i], 1.0e-3)
        ν_l = sp.μ_l / sp.ρ_l
        α_l = sp.k_l / (sp.ρ_l * sp.cp_l)
        Ra  = g * sp.β_l * ΔT * p.height^3 / (ν_l * α_l)
        h_nc = Nu_cc_smooth(Ra) * sp.k_l / p.height

        # Chen (質量流束 G は小さな下限つき。h_mic が主)
        u_c = 0.5 * ((i == 1 ? 0.0 : u.u[i-1]) + u.u[i])
        G = soft_floor(u.ρm[i] * soft_abs(u_c, 1.0e-4), 1.0e-3)
        h_chen = h_chen_smooth(T_ec, Ts, Psat(T_ec, tables), P, tables.Pcrit, u.x[i], G, p.D, sp, sp.h_fg)

        w_on = n_nuc === nothing ? nucleation_onset_weight(T_ec, Ts, sp, r_c, ε_onset) :
                                   min(max(n_nuc[i], 0.0), 1.0)
        h_wet = (1 - w_on) * h_nc + w_on * h_chen

        # 液膜が薄いと沸騰の h は効かない: f_l が f_th_h 付近で、蒸気側 h_dry へ移す
        f_l = (1 - u.x[i]) * u.ρm[i] / tables.ρ_l_sat(u.T[i])
        w_h = w_wet[i] * smooth_step(f_l - f_th_h, ε_h)
        h[i] = w_h * h_wet + (1 - w_h) * h_dry[i]
    end
    return h
end

# -----------------------------------------------------------------------------
# 3. グルーブの残差
# -----------------------------------------------------------------------------

"""
    residual_grv!(resid, du, u, down, Q_ec_grv, S_m_total, S_E_total, p, g, tables; kwargs...)

グルーブ (N セル) の流体残差を resid に書き込む。resid, du, u は
グルーブの ComponentVector (ρm,u,Em,P,T,x)。

  down       : 下流 (蒸気管の最初のセル) のゴースト値。
               transport_ghost_down(u_vl, p_vl, p) で作る
  Q_ec_grv   : groove_evap_heat の戻り値 [W] (セルごとのベクトル)
  S_m_total  : ウィックから蒸発して入る質量流量 [kg/s] (wick_groove_source の第1戻り値)
  S_E_total  : ウィックから入るエネルギー [W] (wick_groove_source の第2戻り値)
  p          : grv_params(geo.grv)

S_m, S_E は、全セルに一様に分配する (蒸発は溝の全長で一様に起きるとする)。
kwargs は residual_transport! へ渡される (ε_u, δ_u, ε_Re, ε_phase)。
"""
function residual_grv!(resid, du, u, down, Q_ec_grv, S_m_total, S_E_total, p, g,
    tables::FluidTables; strat=nothing, kwargs...)

    N = p.N
    # 入口は閉端: u=0、ゴーストはセル1と同じ値 (入口側の勾配をゼロにする)
    up = (ρ = u.ρm[1], Em = u.Em[1], P = u.P[1], u = 0.0, dl = p.dl)

    V_total = N * p.V
    S_m = fill(S_m_total / V_total, N)
    S_E = fill(S_E_total / V_total, N)
    q = Q_ec_grv ./ p.V

    # 成層した出口: 出口面の流束から液の成分を引いた分を、最後のセルへ戻す (groove_stratified_face)
    if strat !== nothing
        S_m[N] += strat.S_m_corr
        S_E[N] += strat.S_E_corr
    end

    return residual_transport!(resid, du, u, up, down, q, p, g, tables;
                               S_m = S_m, S_E = S_E, kwargs...)
end

"""
    groove_stratified_face(u, p, down, g, tables; ε_u, f_full=0.9, ε_full=0.03)
        -> (ρ_eff, Em_eff, S_m_corr, S_E_corr)

グルーブの成層した出口 (重力で、液は底=蒸発器ケース側に溜まり、蒸気だけが出る)。
均質モデルの出口は、セルの平均密度で流出するので、蒸気を流すと、残った液も一緒に流して
しまう (流出が速いほど乾きが速くなる、正のフィードバック)。出口面の流束を、液の成分
(1−x)ρm を除いたものに置き換える:

    ρ_eff = ρm − (1 − s_full) (1−x) ρm,   φ_eff = (ρm φ − (1 − s_full) (1−x) ρm φ_l)/ρ_eff
    φ_l = h_l,sat(T) + Em − em    (液の成分の全エンタルピー)
    s_full = smooth_step(f_l − f_full, ε_full)    (f_l = (1−x)ρm/ρ_l: 液で満たされると 1 → 液も出る)

出口面の流束は、グルーブと蒸気管の最初のセルで共有するので、両方に同じ流束を使う:
  ・蒸気管へ渡すゴースト (上流) は (ρ_eff, Em_eff) に置き換える (flatmain で up_vl を差し替える)
  ・グルーブには、均質の流束 (Φ_hom) と成層の流束 (Φ_str) の差を、最後のセルの質量・
    エネルギーのソースとして戻す (S_m_corr, S_E_corr [kg/m³/s, W/m³])
面の流束は residual_transport! と同じ (風上ブレンド w = smooth_step(u, ε_u)、
F = ρ_f u、Φ = F φ_f)。質量・エネルギーは、グルーブと蒸気管の間で厳密に保存する。
逆流 (u<0) のときは、w≈0 で、グルーブ側の ρ_eff は流束に効かない。
"""
function groove_stratified_face(u, p, down, g, tables::FluidTables;
    ε_u=1.0e-4, f_full=0.9, ε_full=0.03)

    i = p.N
    u_N = u.u[i]
    u_c = 0.5 * ((i == 1 ? 0.0 : u.u[i-1]) + u_N)
    em = internal_energy(u.Em[i], u_c, cell_elevation(p, i), g)
    T_i = u.T[i]
    ρ_i = u.ρm[i]
    φ_i = u.Em[i] + u.P[i] / ρ_i
    φ_l = tables.h_l_sat(T_i) + u.Em[i] - em

    f_l = (1 - u.x[i]) * ρ_i / tables.ρ_l_sat(T_i)
    s_full = smooth_step(f_l - f_full, ε_full)
    m_rm = (1 - s_full) * (1 - u.x[i]) * ρ_i              # 出口から除く液の成分 [kg/m³]

    ρ_eff = soft_floor(ρ_i - m_rm, 1.0e-3)
    φ_eff = (ρ_i * φ_i - m_rm * φ_l) / ρ_eff
    Em_eff = φ_eff - u.P[i] / ρ_eff

    w = smooth_step(u_N, ε_u)
    φ_R = down.Em + down.P / down.ρ
    F_hom = (w * ρ_i + (1 - w) * down.ρ) * u_N
    F_str = (w * ρ_eff + (1 - w) * down.ρ) * u_N
    Φ_hom = F_hom * (w * φ_i + (1 - w) * φ_R)
    Φ_str = F_str * (w * φ_eff + (1 - w) * φ_R)

    return (ρ_eff = ρ_eff, Em_eff = Em_eff,
            S_m_corr = (F_hom - F_str) / p.dl, S_E_corr = (Φ_hom - Φ_str) / p.dl)
end