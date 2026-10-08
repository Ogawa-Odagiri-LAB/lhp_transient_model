# =============================================================================
# types.jl
#
# このファイルは3種類の「型」を定義する。
#
#   1. ジオメトリ構造体      … 既存設計をほぼ踏襲
#   2. 運転条件構造体        … 新しい変数セット (ρm,u,Em,P,T,x,T_wall) 用に
#                              スケーリング定数を更新
#   3. 状態ベクトルの組み立て … ComponentArray として状態を表現するための
#                              コンストラクタ関数、および微分変数/代数変数
#                              マスクのビルダー
#
# 【設計方針】
# 状態 (u, du) はカスタム struct ではなく ComponentArray として表現する。
# 理由:
#   - IDA (Sundials.jl) に渡す際、getdata(u)/getaxes(u) で平坦なVector⇔
#     階層構造を相互変換できる必要があるため (DAEProblemはVector{Float64}
#     しか受け付けない)。
#   - 既存コードもすでに ComponentArrays を利用しており (var.wick, var.grv.T
#     のようなドット記法)、この慣習をそのまま引き継げる。
#   - 「TransportState」のようなカスタム struct を別途定義すると、状態の
#     実体 (ComponentArray) と型定義が二重管理になり、同期ミスの元になる。
#
# 【状態ベクトルの階層構造 (initialize_state が組み立てるもの)】
#
#   u.vl.ρm, u.vl.u, u.vl.Em, u.vl.P, u.vl.T, u.vl.x   … 蒸気管 (N要素)
#   u.cond.*  (同上)                                     … 凝縮器
#   u.ll.*    (同上)                                     … 液管
#   u.grv.*   (同上)                                     … グルーブ
#   u.cc.ρm, u.cc.u, u.cc.Em, u.cc.P, u.cc.T, u.cc.x     … CC
#       (輸送管・グルーブと全く同じ6変数構成。気液分離幾何の角度θ_liqは
#        DAEの未知数に含めず、液体積分率 f_l から θ_liq を求める1次元
#        テーブル (以前議論したケプラー方程式型の関係) をオフラインで
#        構築し、residuals_cc.jl / heattransfer.jl 側で都度補間して
#        使う。これにより、以前懸念した「端点(空・満液)での特異性が
#        Newton反復に直接影響する」リスクを、DAE本体から切り離せる)
#   u.wick.E, u.wick.T                                   … ウィック (N要素)
#       (E: 式46の内部エネルギー, T: 温度。E=微分変数, T=代数変数)
#   u.T_wall.evap                                        … 蒸発器ケース (長さ1の配列)
#       (式41の T_ec そのもの。T_ec = u.T_wall.evap[1]。他の壁と同じ配列に
#        しているのは、ComponentArray のスカラー成分への書き込みや、
#        ./= による一括更新が、配列成分と挙動が違うため)
#   u.T_wall.vl, u.T_wall.cond, u.T_wall.ll              … 各輸送管の壁 (N要素)
#       (流体側と同じ分割・同じ添字に対応する設計)
#   u.T_wall.cc                                          … CCケースの壁 (N要素、通常N=1)
#
#   ※ グルーブ・ウィックには独立した壁構造体を持たせない
#     (グルーブは蒸発器ケース内部にあり、ウィックは式46で固体マトリクスと
#      液を一体化したエネルギーとして扱っているため)。
#
# =============================================================================

using Parameters: @with_kw
using ComponentArrays: ComponentArray

# -----------------------------------------------------------------------------
# 1. ジオメトリ構造体
# -----------------------------------------------------------------------------

# ウィック
@with_kw struct WickGeometry
    N::Int = 5                     # 半径方向の分割数

    n_grv::Int = 12                # グルーブ本数
    L::Float64 = 0.035             # 軸方向長さ [m]

    D_i::Float64 = 0.003           # 内径 [m] (CC側)
    D_o::Float64 = 0.006           # 外径 [m] (グルーブ部を除く)
    dD::Float64 = (D_o - D_i) / N  # 半径方向の分割幅 [m]

    width::Float64 = 0.001         # グルーブ溝幅 [m]
    height::Float64 = 0.001        # グルーブ溝深さ [m]

    # --- 式46 (ウィックエネルギー) で必要な材料物性 ---
    # デフォルト値は与えていない: 誤った物理定数が気づかれないまま
    # 使われることを避けるため、呼び出し側で必ず指定する。
    φ::Float64                     # 空隙率 (0〜1)
    ρ_wick::Float64                # 固体マトリクス密度 [kg/m3]
    Cp_wick::Float64               # 固体マトリクス比熱 [J/(kg・K)]
    # 固相の実効熱伝導率 k_eff は、入力にしない。ウィックのステンレス部の熱伝導率
    # (physics.jl の k_sus304) と、液の熱伝導率 (飽和テーブル) から、空隙率 φ を使って、
    # マクスウェルの複合則 (式53) で計算する (residuals_wick.jl の k_eff_maxwell)。
end

# グルーブ (蒸発器ケース内部、輸送管と同じ6変数構成で扱う)
@with_kw struct GrvGeometry
    N::Int = 1
    n_grv::Int = 12
    L::Float64 = 0.035
    dl::Float64 = L / N
    D_i::Float64 = 0.008           # ケース内径 [m]
    D_o::Float64 = 0.01            # ケース外径 [m]
    width::Float64 = 0.001
    height::Float64 = 0.001
    A_ec_wick::Float64 = width * L * n_grv   # 蒸発器ケース-ウィック界面積 [m2]
    A_ec_grv::Float64 = A_ec_wick            # 蒸発器ケース-グルーブ界面積 [m2]
    A_grv_cross::Float64 = width * height    # グルーブ1本あたりの断面積 [m2]
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
end

# 輸送管 (蒸気管・液管共用)
@with_kw struct TLGeometry
    N::Int
    L::Float64
    dl::Float64 = L / N
    D_i::Float64
    D_o::Float64
    A_cross::Float64 = (π / 4) * D_i^2       # 流路断面積 [m2]
    S_inner::Float64 = π * D_i * dl          # 1セルあたりの内壁表面積 [m2]
    V::Float64 = A_cross * dl                # 1セルあたりの流体体積 [m3]
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
end

# 凝縮器
@with_kw struct CondGeometry
    N::Int
    L::Float64
    dl::Float64 = L / N
    D_i::Float64
    D_o::Float64
    A_cross::Float64 = (π / 4) * D_i^2
    S_inner::Float64 = π * D_i * dl
    V::Float64 = A_cross * dl
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
end

# CC (気液分離幾何を扱うため、半径 R も明示的に保持しておく)
@with_kw struct CCGeometry
    N::Int
    L::Float64                                # 軸方向長さ [m] (geometry_cc.jl の Lc)
    dl::Float64 = L / N
    D_i::Float64                              # 内径 [m]
    D_o::Float64                              # 外径 [m]
    R::Float64 = D_i / 2                      # 内半径 [m] (水平円筒の気液分離計算用)
    A_cross::Float64 = (π / 4) * D_i^2
    S_inner::Float64 = π * D_i * dl
    V::Float64 = A_cross * dl
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
end

@with_kw struct LHPGeometry
    wick::WickGeometry
    grv::GrvGeometry
    vl::TLGeometry
    cond::CondGeometry
    ll::TLGeometry
    cc::CCGeometry
end

# -----------------------------------------------------------------------------
# 2. 運転条件構造体
# -----------------------------------------------------------------------------

@with_kw struct LHPTransientCondition
    # --- 計算条件 ---
    Q_apply::Float64                   # ステップ後の印加熱量 [W] (最終値)
    t_init::Float64                    # 計算開始時刻 [s]
    t_end::Float64                     # 計算終了時刻 [s]
    saveat::Float64 = 1.0              # 結果を保存する時間間隔 [s]
                                        # (IDAの内部刻み幅とは独立。出力間引き用)

    # --- 印加熱量のステップアップ (滑らかなステップ) ---
    # Q(t) = Q_apply_init + (Q_apply − Q_apply_init) s((t − t_step)/τ_step)
    #   s(x) = 0 (x≤0),  3x²−2x³ (0<x<1),  1 (x≥1)         (3次の smoothstep)
    # t ≤ t_step では厳密に Q_apply_init、t ≥ t_step+τ_step では厳密に Q_apply で、
    # その間を C¹ (1階微分が連続) につなぐ。不連続なステップは IDA の積分を難しく
    # するので、有限の立ち上がり時間 τ_step を持たせる。tanh 形と違い、t=t_init で
    # 厳密に Q_apply_init から始まり、区切り点 (t_step, t_step+τ_step) も決まるので、
    # IDA の tstops に渡せる (applied_heat_breakpoints)。
    # Q_apply_init == Q_apply なら、ステップなしの一定値になる。
    Q_apply_init::Float64 = Q_apply    # ステップ前の印加熱量 [W] (既定: ステップなし)
    t_step::Float64 = t_init           # ステップ (立ち上がり) の開始時刻 [s]
    τ_step::Float64 = 1.0              # 立ち上がりにかかる時間 [s]

    # --- 作動流体・境界条件 ---
    fluid::String = "ammonia"
    T_cond::Float64 = 15.0 + 273.15    # 凝縮器境界温度 [K]
    T_amb::Float64 = 25.0 + 273.15     # 周囲環境温度 [K]
    T_init::Float64 = 25.0 + 273.15    # 初期の基準温度 [K]。CC・ウィック・蒸発器ケースの温度、
                                        # ウィックの液密度の基準温度にも使う
    ΔT_startup::Float64 = 0.0          # 起動判定: |T_wick,out側 − T_CC| がこの値を超えたら蒸発を許す [K]
                                        # (0 なら判定なし。旧コードの fluid_wick の起動判定に相当)
    ε_startup::Float64 = 0.1           # 起動判定のゲート (smooth_step) の幅 [K]
    g::Float64 = 9.8                   # 重力加速度 [m/s2]

    # --- 整合初期条件 (initialcondition.jl) ---
    # 定常状態は求めず、静止 (u=0) の平衡に近い状態から始める。液・蒸気は、表の補間
    # 誤差の影響を受けやすい飽和線の上を避けるため、基準温度 T_init から少しずらす
    # (液は過冷却、蒸気は過熱)。壁は、接する流体と同じ温度から始める。
    ΔT_sub_init::Float64 = 1.0         # 液 (凝縮器・液管) の過冷却度 [K]。正の値
    ΔT_sup_init::Float64 = 1.0         # 蒸気 (グルーブ・蒸気管) の過熱度 [K]。正の値
    f_l_cc_init::Float64 = 0.5         # CCの初期の液体積分率 (0〜1)。CCは T_init の気液二相
    cond_init_phase::Symbol = :liquid  # 凝縮器の初期状態 (:liquid または :vapor)

    # --- 熱伝達係数 (定数) ---
    h_ev::Float64  = 10000.0           # 蒸発熱伝達係数 (Kaya[14]を参考) [W/m2K]
    h_con::Float64 = 500.0             # 凝縮器プレートとの接触熱伝達係数 (式51) [W/m2K]

    # --- 風上ブレンドの幅 ---
    ε_u::Float64 = 1.0e-4              # 輸送管の速度の判定幅 [m/s] (流速の典型値より十分小さく)

    # --- スケーリング定数 ---
    # IDAのNewton反復・誤差評価は無次元の残差ノルムに基づくため、
    # 変数間の桁違い (P~1e6, x~O(1), T~1e2 など) を解消しておくと
    # 収束性が大きく改善する。state_scaleのみ用い、scale_var/restore_var
    # (以前のTrustRegion実装) と同じ考え方を新変数セットに合わせて再定義する。
    ρm_scale::Float64 = 50.0           # [kg/m3]  (気液平均密度の代表値)
    u_scale::Float64  = 10.0           # [m/s]
    Em_scale::Float64 = 4.0e5          # [J/kg]   (全エネルギーの代表値)
    P_scale::Float64  = 1.0e6          # [Pa]
    T_scale::Float64  = 300.0          # [K]
    x_scale::Float64  = 1.0            # [-]      (乾き度はもともとO(1))

    # --- 反復対象コンポーネントの列挙 (ループ処理の共通化に利用) ---
    fluid_comp::NTuple{5,Symbol} = (:grv, :vl, :cond, :ll, :cc)  # TransportState系
    wall_comp::NTuple{5,Symbol}  = (:evap, :vl, :cond, :ll, :cc) # T_wall系
end

# -----------------------------------------------------------------------------
# 印加熱量 Q_apply(t)
# -----------------------------------------------------------------------------

"""
    smoothstep3(x)

3次の smoothstep: x≤0 で 0、x≥1 で 1、その間は 3x²−2x³ (C¹ 連続)。
"""
smoothstep3(x) = x <= 0 ? 0.0 : (x >= 1 ? 1.0 : x^2 * (3 - 2x))

"""
    applied_heat(t, condition) -> Q_apply(t) [W]

時刻 t の印加熱量。Q_apply_init から Q_apply へ、t_step から τ_step かけて
滑らかにステップアップする (LHPTransientCondition のコメント参照)。
引数 t は時刻 (状態変数ではないので、ここでの区間分けの分岐は、残差の
状態依存の不連続を生まない)。
"""
applied_heat(t, c::LHPTransientCondition) =
    c.Q_apply_init + (c.Q_apply - c.Q_apply_init) * smoothstep3((t - c.t_step) / c.τ_step)

"""
    applied_heat_breakpoints(condition) -> Vector{Float64}

Q_apply(t) の区切り点 (立ち上がりの開始・終了時刻)。IDA に tstops として渡すと、
ここで積分が確実に刻まれ、C¹ の接続点を刻み幅が飛び越えない。
Q_apply_init == Q_apply (ステップなし) のときは空。計算区間外の点は除く。
"""
function applied_heat_breakpoints(c::LHPTransientCondition)
    c.Q_apply_init == c.Q_apply && return Float64[]
    pts = [c.t_step, c.t_step + c.τ_step]
    return filter(t -> c.t_init < t <= c.t_end, pts)
end

# -----------------------------------------------------------------------------
# 3. 状態ベクトルの組み立て
# -----------------------------------------------------------------------------

"""
    _zeros_transport(N)

輸送管・グルーブ用の6変数 (ρm,u,Em,P,T,x) をすべてゼロ初期化した
NamedTuple を返す。ComponentArray の入れ子として使う。

要素型は Float64 に固定している。疎行列パターン (jac_prototype) は
自動検出ツールに頼らず手動で構築する方針としたため、Dual数対応は
不要と判断した。
"""
_zeros_transport(N) = (ρm=zeros(N), u=zeros(N), Em=zeros(N), P=zeros(N), T=zeros(N), x=zeros(N))

"""
    initialize_state(geo::LHPGeometry) -> ComponentArray

与えられたジオメトリに対応する形状の状態ベクトルを、すべてゼロで
初期化して返す。実際の値は initialcondition.jl 側の
build_consistent_ic で、環境温度・封入量などから物理的に設定する。
"""
function initialize_state(geo::LHPGeometry)
    return ComponentArray(
        vl   = _zeros_transport(geo.vl.N),
        cond = _zeros_transport(geo.cond.N),
        ll   = _zeros_transport(geo.ll.N),
        grv  = _zeros_transport(geo.grv.N),
        cc   = _zeros_transport(geo.cc.N),   # θ_liqは状態に持たず、residuals_cc.jl内で
                                              # f_l → θ_liq の1次元テーブル補間として
                                              # その場で計算する(DAEの未知数には含めない)
        wick = (E = zeros(geo.wick.N), T = zeros(geo.wick.N)),
        T_wall = (
            evap = zeros(1),
            vl   = zeros(geo.vl.N),
            cond = zeros(geo.cond.N),
            ll   = zeros(geo.ll.N),
            cc   = zeros(geo.cc.N),
        ),
    )
end

"""
    initialize_differential_vars(geo::LHPGeometry) -> ComponentArray{Bool}

initialize_state と同じ形状で、各スロットが微分変数 (true) か
代数変数 (false) かを示す Bool の ComponentArray を返す。
DAEProblem の differential_vars にそのまま渡す。

対応関係 (これまでの議論で確定した分類):
    ρm, u, Em        … 微分変数 (質量・運動量・エネルギー保存則)
                       (ただしCCのuは代数変数。CCは貯水池として扱い、運動量式を使わない)
    P, T, x          … 代数変数 (CoolPropのflash拘束)
    wick.E           … 微分変数 (式47のエネルギー保存則)
    wick.T           … 代数変数 (式46とCoolPropの整合性拘束)
    T_wall.*         … すべて微分変数 (式20/41、代数変数を伴わない純粋なODE)
"""
function initialize_differential_vars(geo::LHPGeometry)
    transport_mask(N) = (ρm=trues(N), u=trues(N), Em=trues(N),
                          P=falses(N), T=falses(N), x=falses(N))

    return ComponentArray(
        vl   = transport_mask(geo.vl.N),
        cond = transport_mask(geo.cond.N),
        ll   = transport_mask(geo.ll.N),
        grv  = transport_mask(geo.grv.N),
        # CCは貯水池 (0次元) として扱う (residuals_cc.jl)。運動量式は使わず、
        # 出口面の速度 u は代数拘束 u = ṁ_out/(ρ_l A) で決まるため、
        # ρm, Em だけが微分変数で、u は P, T, x と同じ代数変数。
        cc   = merge(transport_mask(geo.cc.N), (u = falses(geo.cc.N),)),
        wick = (E = trues(geo.wick.N), T = falses(geo.wick.N)),
        T_wall = (
            evap = trues(1),
            vl   = trues(geo.vl.N),
            cond = trues(geo.cond.N),
            ll   = trues(geo.ll.N),
            cc   = trues(geo.cc.N),
        ),
    )
end

"""
    scale_state(u, condition::LHPTransientCondition) -> ComponentArray

状態ベクトルを無次元化する。IDAに渡す前に適用する。
対応する復元関数 restore_state とペアで使うこと。
"""
function scale_state(u, condition::LHPTransientCondition)
    scaled = copy(u)
    for comp in condition.fluid_comp
        c = getproperty(scaled, comp)
        c.ρm ./= condition.ρm_scale
        c.u  ./= condition.u_scale
        c.Em ./= condition.Em_scale
        c.P  ./= condition.P_scale
        c.T  ./= condition.T_scale
        c.x  ./= condition.x_scale
    end
    scaled.wick.E ./= condition.Em_scale   # ウィックのEも同じエネルギースケールを流用
    scaled.wick.T ./= condition.T_scale
    for comp in condition.wall_comp
        c = getproperty(scaled.T_wall, comp)
        c ./= condition.T_scale
    end
    return scaled
end

"""
    restore_state(u_scaled, condition::LHPTransientCondition) -> ComponentArray

scale_state の逆変換。
"""
function restore_state(u_scaled, condition::LHPTransientCondition)
    restored = copy(u_scaled)
    for comp in condition.fluid_comp
        c = getproperty(restored, comp)
        c.ρm .*= condition.ρm_scale
        c.u  .*= condition.u_scale
        c.Em .*= condition.Em_scale
        c.P  .*= condition.P_scale
        c.T  .*= condition.T_scale
        c.x  .*= condition.x_scale
    end
    restored.wick.E .*= condition.Em_scale
    restored.wick.T .*= condition.T_scale
    for comp in condition.wall_comp
        c = getproperty(restored.T_wall, comp)
        c .*= condition.T_scale
    end
    return restored
end