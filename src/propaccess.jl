# =============================================================================
# propaccess.jl (Level 2: 物性テーブルへの統一アクセス関数)
#
# proptables.jl で構築した FluidTables を、残差関数 (residuals_*.jl) から
# 直接呼びやすい形でラップする。残差関数は、このファイルの関数だけを
# 呼べばよく、FluidTables の内部構造 (補間オブジェクトの種類など) を
# 一切気にしなくてよいようにする。
#
# 【性能上の注意】
# FluidTables の各補間オブジェクトは、構造体の定義を簡潔に保つため
# 型を Any で保持している (proptables.jl 参照)。そのため
# `tables.P_of_ρe(ρ,e)` のような呼び出しは、Juliaの型推論が効かず、
# 動的ディスパッチ (実行時に型を調べてから呼ぶ) になる。
# CoolPropの直接呼び出しに比べれば十分高速なはずだが、もしプロファイル
# の結果この部分がボトルネックになるようであれば、FluidTablesを
# 補間オブジェクトの型でパラメトリックにする (例: FluidTables{F1,F2,...})
# ことで動的ディスパッチを解消できる。現時点では可読性を優先している。
# =============================================================================

# -----------------------------------------------------------------------------
# 1. テーブル取得
# -----------------------------------------------------------------------------

"""
    get_tables(fluid_name::String) -> FluidTables

FLUID_TABLES から該当流体のテーブルを取り出す。init_fluid_lookups で
事前に構築・登録されていない場合は、対処法を示したエラーを出す。
"""
function get_tables(fluid_name::String)
    haskey(FLUID_TABLES, fluid_name) ||
        error("流体 \"$fluid_name\" のテーブルが未構築です。先に " *
              "init_fluid_lookups(\"$fluid_name\"; T_min=…, T_max=…, P_max=…) " *
              "を呼んでテーブルを構築・登録してください。")
    return FLUID_TABLES[fluid_name]
end

# -----------------------------------------------------------------------------
# 2. (ρ,e) flash への窓口
# -----------------------------------------------------------------------------

"""
    flash_ρe(ρ, e, tables::FluidTables) -> (P, T, x)

(ρ,e) から (P,T,x) を2次元補間テーブルを通じて求める。
残差関数が式4-6 (P=P(ρm,em), T=T(ρm,em), x=x(ρm,em)) の代数拘束を
評価する際、このテーブル検索を直接呼ぶのではなく、必ずこの関数を
経由すること (テーブルの内部構造が変わっても残差関数側を変更せずに
済むようにするため)。
"""
# P(ρ,e) は、飽和液線 (δ = ρ − ρ_sl(e) = 0) で傾きが約 1.4 万倍に跳ねる折れを持つ
# (二相側 dP/dρ ≈ 129、液側 ≈ 1.84e6 Pa/(kg/m³))。前線セルがこの折れをまたぐと、Newton が収束しない。
# |δ| < ε の範囲だけ、両側の傾きを保ったまま 2 次曲線でつなぎ、傾きを連続にする。
# 飽和線から ε 以上離れた領域は変わらない。飽和線上の P は最大で (s_液 − s_二相) ε/4 だけ高くなる
# (ε = 0.002 kg/m³ で約 0.9 kPa)。ε ≤ 0 で無効 (従来の折れ)。
const P_SMOOTH_EPS = Ref(0.002)

function P_smoothed(ρ, e, tables::FluidTables)
    f = tables.P_of_ρe
    ε = P_SMOOTH_EPS[]
    δ = ρ - f.ρ_ref(e)
    (ε > 0 && abs(δ) < ε) || return f(ρ, e)
    P_a = f.itp(e, -ε)
    P_0 = f.itp(e, zero(ε))
    P_b = f.itp(e, ε)
    s_a = (P_0 - P_a) / ε
    s_b = (P_b - P_0) / ε
    return P_a + s_a * (δ + ε) + (s_b - s_a) * (δ + ε)^2 / (4ε)
end

flash_ρe(ρ, e, tables::FluidTables) = (
    P = P_smoothed(ρ, e, tables),
    T = tables.T_of_ρe(ρ, e),
    x = tables.x_of_ρe(ρ, e),
)

"""
    transport_props(ρ, e, tables::FluidTables) -> (μ, k, cp, Pr)

(ρ,e) における単相域の輸送物性 (粘度・熱伝導率・比熱) を、
flash_ρe と同じ(ρ,e)2次元テーブルから求める。Prandtl数も
その場で組み立てて一緒に返す (Pr = μ*cp/k)。

sat_props (飽和線上のみ) では単相域の値が得られないため、単相の
強制対流 (式50) や摩擦係数の計算で局所粘度などが必要な場合は、
必ずこちらを使うこと。
"""
function transport_props(ρ, e, tables::FluidTables)
    # 表の補間・外挿で、物性が負になっても、Pr^n (小数べき) が定義域エラーに
    # ならないよう、下限をつける。下限は典型値より桁違いに小さく (μ: 蒸気1e-5,
    # 液1e-4 Pa s に対し 1e-8、k: 0.02〜0.5 W/mK に対し 1e-4、cp: 1e3〜1e4 J/kg/K
    # に対し 10) 通常の状態の値は変わらない (相対変化 1e-5 以下)。Newton 反復の途中で、状態が
    # 表の範囲外 (線形外挿) に出たときの安全装置。
    μ  = soft_floor(tables.μ_of_ρe(ρ, e), 1.0e-8)
    k  = soft_floor(tables.k_of_ρe(ρ, e), 1.0e-4)
    cp = soft_floor(tables.cp_of_ρe(ρ, e), 10.0)
    return (μ = μ, k = k, cp = cp, Pr = μ * cp / k)
end

"""
    internal_energy(Em, u, elev, g=9.8) -> em

全内部エネルギー Em (式9: Em = em + u^2/2 + g z sinθ) から、
比内部エネルギー em を逆算する。輸送管・グルーブ・CCで共通して
使う変換であり、P を一切使わずに計算できる点が重要
(以前議論した「循環参照の回避」。h ではなく e を使う設計のため、
 この変換にPは登場しない)。

位置エネルギー項は、式9の z sinθ を「そのセルの標高 elev [m]」として
受け取る。Em は流れによって運ばれる量なので、この項はコンポーネント
をまたいで連続な標高でなければならない (各管が別々の z の原点を持つと、
コンポーネント境界で Em が不連続になり、見かけのエネルギー生成が出る)。
直管なら elev = elev0 + (セル中心までの距離) sinθ (residuals_transport.jl の
cell_elevation 参照)。傾斜角 θ は運動量式の重力項 g sinθ にだけ現れる。
"""
internal_energy(Em, u, elev, g=9.8) = Em - u^2 / 2 - g * elev

"""
    flash_transport_state(Em, u, ρ, elev, tables::FluidTables; g=9.8) -> (P, T, x)

輸送管・グルーブ・CCの状態 (ρ,u,Em) から、直接 (P,T,x) を求める。
internal_energy による Em→em の変換と flash_ρe の呼び出しをまとめた、
残差関数側で最も頻繁に使うことになる便利関数。
"""
function flash_transport_state(Em, u, ρ, elev, tables::FluidTables; g=9.8)
    em = internal_energy(Em, u, elev, g)
    return flash_ρe(ρ, em, tables)
end

# -----------------------------------------------------------------------------
# 3. 1次元飽和物性への窓口
# -----------------------------------------------------------------------------

"""
    sat_props(T, tables::FluidTables)

Lockhart-Martinelli型の相関式、Chenの式、CC自然対流のRa計算などで
必要となる、飽和液・飽和蒸気の物性をまとめて返す。
"""
function sat_props(T, tables::FluidTables)
    h_l = tables.h_l_sat(T)
    h_v = tables.h_v_sat(T)
    μ_l = tables.μ_l_sat(T)
    μ_v = tables.μ_v_sat(T)
    k_l = tables.k_l_sat(T)
    k_v = tables.k_v_sat(T)
    cp_l = tables.cp_l_sat(T)
    cp_v = tables.cp_v_sat(T)
    return (
        μ_l = μ_l, μ_v = μ_v,
        ρ_l = tables.ρ_l_sat(T), ρ_v = tables.ρ_v_sat(T),
        k_l = k_l, k_v = k_v,
        cp_l = cp_l, cp_v = cp_v,
        Pr_l = μ_l * cp_l / k_l,   # Shah/Chenの式で必要な飽和液のPrandtl数
        Pr_v = μ_v * cp_v / k_v,   # Dougall-Rohsenowの式で必要な飽和蒸気のPrandtl数
        h_l = h_l, h_v = h_v,
        h_fg = h_v - h_l,
        σ = tables.σ_sat(T),
        β_l = tables.β_l_sat(T),
        β_v = tables.β_v_sat(T),
    )
end

"""
    Psat(T, tables::FluidTables)

飽和温度Tにおける飽和圧力。Chenの式(60)の P_sat(T_wall) や、
凝縮・沸騰の判定 (T_wall と T_sat(P) の比較) で使う。
"""
Psat(T, tables::FluidTables) = tables.Psat_T(T)

"""
    Tsat(P, tables::FluidTables)

飽和圧力Pにおける飽和温度。ウィック最外層 (式48周辺、
groove圧力から飽和温度を求める場面) で使う。
"""
Tsat(P, tables::FluidTables) = tables.Tsat_P(P)

# -----------------------------------------------------------------------------
# 3.5 圧縮液 (ウィック内の液など) の物性
# -----------------------------------------------------------------------------

"""
    liquid_props(T, P, tables::FluidTables) -> (ρ, e, h)

温度T・圧力Pの液の密度・比内部エネルギー・比エンタルピー。
(ρ,e)の2次元テーブルは、(T,P)から引く向きには使えないため、
飽和液の1次元テーブルから近似する。

    ρ = ρ_l,sat(T)                      (液の圧縮率は小さく、圧力依存は無視)
    e = h_l,sat(T) - P_sat(T)/ρ_l,sat(T)  (液の内部エネルギーは圧力にほぼ依存しない)
    h = e + P/ρ

圧力Pの影響は h に P/ρ としてだけ入る。ウィック内のように圧力が一様な領域では、
隣り合うセルの h の差にはほぼ影響しない。Pは飽和圧力より高い (圧縮液) か、
近い値を想定している。
"""
function liquid_props(T, P, tables::FluidTables)
    ρ = tables.ρ_l_sat(T)
    e = tables.h_l_sat(T) - tables.Psat_T(T) / ρ
    return (ρ = ρ, e = e, h = e + P / ρ)
end

# -----------------------------------------------------------------------------
# 4. 疑似乾き度 (単相⇔二相の滑らかな接続用)
# -----------------------------------------------------------------------------

"""
    pseudo_quality(em, ρ, P, T, tables::FluidTables) -> x_pseudo

熱伝達係数のディスパッチ (Level4, heattransfer.jl) で、単相⇔二相の
相関式をtanhブレンドするために使う「疑似乾き度」。

flash_ρe が返す x (CoolPropの二相flashそのものの乾き度) は、単相域
では0または1にクランプされた値になる (proptables.jl の
flash_ρe_once 参照) ため、この値をそのままブレンド変数として使うと
単相域の内部でブレンド重みが変化せず、不自然になる。

pseudo_quality はその場の T における飽和液・飽和蒸気エンタルピーを
基準に、単相域でも連続的に0未満・1超へ滑らかに伸びる値を与える
(以前議論した「過熱蒸気・過冷液でも滑らかに繋がる疑似乾き度」)。

flash_ρe で P,T が既に得られていることを前提とするため、
h_m = em + P/ρ の計算にPを使っても循環参照にはならない。
"""
function pseudo_quality(em, ρ, P, T, tables::FluidTables)
    h_m = em + P / ρ
    h_l = tables.h_l_sat(T)
    h_v = tables.h_v_sat(T)
    return (h_m - h_l) / (h_v - h_l)
end

"""
    pseudo_quality_P(em, ρ, P, tables) -> x_ps

熱伝達の単相⇔二相のブレンド用の疑似乾き度。`pseudo_quality` はそのセル自身の温度 T での
飽和エンタルピーを基準にするため、単相液では過冷却度に関係なく ≈0、単相蒸気では ≈1 になる
(過冷却液にも二相の h が半分掛かる)。こちらは圧力の飽和温度 Tsat(P) を基準にするので、
過冷却液は x_ps < 0 (過冷却度 1 K で約 −0.004)、過熱蒸気は x_ps > 1 に伸びる。
"""
function pseudo_quality_P(em, ρ, P, tables::FluidTables)
    h_m = em + P / ρ
    Ts = Tsat(P, tables)
    h_l = tables.h_l_sat(Ts)
    h_v = tables.h_v_sat(Ts)
    return (h_m - h_l) / (h_v - h_l)
end