# =============================================================================
# proptables.jl
#
# 2次元 (ρ,e) flashテーブル (単相・二相を同じ入力ペアで扱う) と、
# 1次元飽和物性テーブル (T依存) を構築する。
#
# 【入力ペアについて】
# 以前の設計検討の過程で、エネルギーの状態変数を全エンタルピー
# (E=h_m+u^2/2+gzsinθ) として扱う案と、全内部エネルギー
# (E_m=e_m+u^2/2+gzsinθ) として扱う案を比較した。前者は貯蔵項に
# 本来不要な∂P/∂t項が必要になる問題があったため、最終的に後者
# (全内部エネルギーEmを状態変数とし、e_m=Em-u^2/2-gzsinθ を
# flashの入力とする) で確定している (設計ドキュメント式4-6,8-9)。
# したがって、このテーブルの入力ペアは (ρ,h) ではなく (ρ,e) である。
#
# 【重要な訂正】
# これまでの会話で提示してきたJuliaコード例の多くは
#     state = CoolProp.AbstractState("HEOS", fluid)
#     CoolProp.update(state, CoolProp.DmassUmass_INPUTS, ρ, e)
#     P = CoolProp.PropsSI(state, "P")
# という書き方をしていましたが、これは Python版CoolPropのAPIであり、
# Julia版 (CoolProp.jl) の実際の低レベルインターフェースとは異なります。
# 正しいJulia側のAPIは次の通りです (ハンドルはInt、関数名も異なる):
#
#     handle = CoolProp.AbstractState_factory(backend, fluid)   # ハンドル取得
#     CoolProp.AbstractState_update(handle, input_pair_index, v1, v2)
#     val = CoolProp.AbstractState_keyed_output(handle, param_index)
#     CoolProp.AbstractState_free(handle)                       # 必ず解放する
#
# input_pair_index / param_index は、それぞれ
#     CoolProp.get_input_pair_index("DmassUmass_INPUTS")
#     CoolProp.get_param_index("P")
# のように文字列から取得する。本ファイルはこの正しいAPIに統一している。
#
# 【2次元テーブルの座標: (ρ,e) ではなく (e, δ=ρ−ρ_sl(e)) 】
# 当初は (ρ,e) の矩形格子でテーブル化したが、圧縮液で精度が全く足りなかった。
# 液体は圧縮性が極端に小さく、圧力が 1→4 MPa と変わっても密度は 0.5% しか
# 変わらない。蒸気から液までを1つの矩形で覆う格子は、液体の圧力変化の全幅が
# 格子1マスに収まってしまい、圧縮液の P の誤差が 25〜180%、∂P/∂ρ が厳密値の
# 0.4〜0.9倍になった (test_table_accuracy.jl で確認)。
# そこで、飽和液線からの密度のずれ
#     δ = ρ − ρ_sl(e)
# (ρ_sl(e): 飽和液線上で内部エネルギーが e となる飽和液の密度) を座標にとる。
#   - 圧縮液は 0 ≤ δ ≲ 数 kg/m³ の薄い帯で、P はほぼ δ に線形
#     (P − P_sat ≈ δ/(ρ κ_T))。粗い格子でも線形補間が成立する。
#   - 気液の境界 (飽和液線) がちょうど δ=0 の格子線になるので、境界を
#     またぐセルがない。境界での圧力の傾きの急変 (液体の硬さと二相の
#     柔らかさの差) を、補間が正確に表せる。
# 圧縮液の ∂P/∂ρ の比 (テーブル/厳密) は 0.4〜0.9 から 0.993〜1.006 に、
# T の誤差は 0.35 K から 0.000 K になる。二相・蒸気の精度は保たれる。
# 呼び出し側の API (tables.P_of_ρe(ρ,e) など) は変わらない
# (RhoEInterp が (ρ,e)→(e,δ) の変換を内部で行う)。
#
# 【設計方針の再確認】
#   - (ρ,e) 入力は単相・二相を区別せず一意に状態を決定できる
#     (以前確認した HSU_D_flash / HSU_D_flash_twophase に対応。
#      内部エネルギーUも対応する入力として扱われる)。
#   - CoolProp組み込みのTTSE/Bicubicテーブルは (D,U) 入力に非対応のため、
#     自前でJulia側にテーブルを構築する。
#   - 2次元テーブルは「オフラインで1回だけ構築し、以降は補間だけで済ます」
#     という方針。構築コストが高いため、Serialization でキャッシュする。
#   - 1次元飽和テーブルは T のみに依存する軽量な計算なので、毎回
#     高レベルAPI (PropsSI) で素直に計算してよい。
# =============================================================================

using Interpolations: interpolate, scale, extrapolate, Line, BSpline, Linear, Gridded
using Serialization: serialize, deserialize

# -----------------------------------------------------------------------------
# 0. テーブル本体の型
# -----------------------------------------------------------------------------

"""
    FluidTables

作動流体の物性テーブル。補間オブジェクトのフィールドは、型パラメータで具体的な型を
持たせている (RE, RD, RT: 範囲、I2: (e,δ) の2次元補間 (RhoEInterp)、I1: 1次元の補間、
IG: 非等間隔 (Gridded) の1次元補間)。

【型パラメータにする理由】 以前は、補間オブジェクトのフィールドを `::Any` にしていた
(型が長大になるため)。しかし、`tables.μ_l_sat(T)` のように、`Any` 型のフィールドを関数として
呼ぶと、Julia は、コンパイル時に型が分からず、呼び出しのたびに動的ディスパッチと、結果の
ヒープ確保 (箱詰め) を行う。残差評価の 1 回で、補間を数千回呼ぶので、プロファイルでは、
NamedTuple の生成だけで 12% の自己時間、メモリ確保が 590 KiB/回 (650 変数) あった。
型パラメータにすると、型が確定し、動的ディスパッチと箱詰めがなくなる。

型パラメータの対応 (全て、build 時に同じ型になる):
  RE, RD, RT … e_range, δ_range, T_sat_range (StepRangeLen)
  I2 … P/T/x/μ/k/cp_of_ρe (RhoEInterp{補間, ρ_ref}) の 6 個
  I1 … Psat_T と、*_sat の 13 個 (make_interp1)
  IG … ρ_ref と Tsat_P (make_interp1_nonuniform)
"""
struct FluidTables{RE,RD,RT,I2,I1,IG}
    fluid_name::String
    Pcrit::Float64
    Tcrit::Float64

    # --- 2次元 flashテーブル (単相・二相共通)。座標は (e, δ=ρ−ρ_sl(e))。
    #     呼び出しは (ρ, e) のまま (RhoEInterp が座標変換する) ---
    e_range::RE       # e: 比内部エネルギー [J/kg]
    δ_range::RD       # δ: 飽和液線からの密度のずれ [kg/m³] (0 が格子点)
    ρ_ref::IG         # ρ_sl(e): 飽和液線上の密度 (e の1次元補間)
    P_of_ρe::I2
    T_of_ρe::I2
    x_of_ρe::I2
    # --- 単相域の輸送物性 (同じ(ρ,e)格子上。sat_propsは飽和線上のみ
    #     なので、単相の強制対流・摩擦係数の計算にはこちらを使う) ---
    μ_of_ρe::I2
    k_of_ρe::I2
    cp_of_ρe::I2

    # --- 1次元飽和テーブル (T依存) ---
    T_sat_range::RT
    Psat_T::I1
    Tsat_P::IG
    μ_l_sat::I1; μ_v_sat::I1
    ρ_l_sat::I1; ρ_v_sat::I1
    k_l_sat::I1; k_v_sat::I1
    cp_l_sat::I1; cp_v_sat::I1
    h_l_sat::I1; h_v_sat::I1
    σ_sat::I1
    β_l_sat::I1   # 定圧膨張係数 (CC自然対流のRa計算用、液側)
    β_v_sat::I1   # 定圧膨張係数 (CC自然対流のRa計算用、蒸気側)
end

const FLUID_TABLES = Dict{String,FluidTables}()

"""
    RhoEInterp(itp, ρ_ref)

(e, δ) 座標の2次元補間 itp を、(ρ, e) で呼べるようにする小さな構造体。
呼び出し f(ρ, e) は itp(e, ρ − ρ_ref(e)) を返す。FluidTables の
P_of_ρe, T_of_ρe, x_of_ρe, μ_of_ρe, k_of_ρe, cp_of_ρe がこの型で、
propaccess.jl や各テストは座標の違いを意識せずに (ρ, e) で呼べる。
型パラメータで補間オブジェクトの型を保持するので、呼び出しは型安定。
"""
struct RhoEInterp{I,R}
    itp::I
    ρ_ref::R
end
function (f::RhoEInterp)(ρ, e)
    # NaN/Inf を Interpolations に渡すと、添字計算が不正になり @inbounds の範囲外読み出し
    # (EXCEPTION_ACCESS_VIOLATION = Julia ごと落ちる) になる。必ずここで止める。
    (isfinite(ρ) && isfinite(e)) ||
        throw(DomainError((ρ, e), "RhoEInterp: 非有限の (ρ, e) が入力された"))
    return f.itp(e, ρ - f.ρ_ref(e))
end

# キャッシュ (.jls) の中身や、FluidTables のフィールド構成を変えたら、この値を上げる。
# 古いキャッシュは、ファイル名に含めて自動的に無効化する (構築し直される)。
#   v5 まで: FluidTables (型 + 補間オブジェクト) を、そのまま serialize していた。
#            Serialization は、型が定義されたモジュールの名前 (LHPTransientModel.FluidTables
#            など) を保存するので、モジュールの中 (transient.jl) で作ったキャッシュを、
#            ファイルを直接 include するテストが読めず (KeyError: key LHPTransientModel
#            not found)、逆も同様だった。
#   v6 以降: 素のデータ (配列・範囲・数値だけ) を保存し、読み込み後に補間オブジェクトと
#            FluidTables を組み立てる。モジュールに依存する型が入らないので、
#            どこで作ったキャッシュも、どこからでも読める。
const TABLE_VERSION = "v6"

# -----------------------------------------------------------------------------
# 1. 低レベルAPIの薄いラッパー
# -----------------------------------------------------------------------------

"""
    with_abstract_state(f, backend, fluid)

AbstractStateのハンドルを生成し、関数 f(handle) を実行した後、
例外が起きても確実にハンドルを解放する (try/finally)。
2次元テーブル構築のように同じ流体に対して大量の呼び出しを行う場合、
毎回ハンドルを生成・解放するコスト (高レベルPropsSIの文字列パース
オーバーヘッドの主因) を避けるため、1つのハンドルを使い回す。
"""
function with_abstract_state(f, backend::String, fluid::String)
    handle = CoolProp.AbstractState_factory(backend, fluid)
    try
        return f(handle)
    finally
        CoolProp.AbstractState_free(handle)
    end
end

"""
    flash_ρe_once(handle, ρ, e, fluid_name, idx_DU, idx_P, idx_T, idx_Q, idx_V, idx_L, idx_C)

(ρ,e) から (P,T,x,μ,k,cp) を1点だけflashする。e は比内部エネルギー。
μ,k,cp (粘度・熱伝導率・比熱) は、単相域の強制対流・摩擦係数の
計算に必要な輸送物性で、sat_props (飽和線上のみ) では代用できない
ため、flash本体と同じ(ρ,e)格子上でまとめて取得する。

単相域でのQの扱いはCoolPropのバックエンド・バージョンによって
異なりうる (例外を投げる場合、意味のない値を返す場合など) ため、
Qの取得に頼らず「Qが[0,1]の範囲に収まらなければ単相とみなし、
Tと飽和温度の大小関係から0/1を確定させる」という、バックエンドに
依存しない判定方法を採用している。

【二相域の μ,k,cp】 CoolProp が二相の状態に返す値は、熱力学的に意味が
定まらず、特に c_p は、二相では T が P だけで決まる (飽和) ため、定義
(∂h/∂T)_P が不定で、符号すら定まらない (アンモニア 298 K で、x=0.1 で
c_p = −669 J/kg/K、x=0.5 で 4953 など)。この値を表にそのまま入れると、
(a) Pr = μ c_p/k が負になり、Nu の Pr^n が DomainError になる、
(b) 飽和液線の近くで、単相の正しい値が、隣の二相の無意味な値と補間で
    混ざって狂う、という問題が起きた。
そこで、二相のときは、同じ温度の飽和液・飽和蒸気の値の、乾き度による
線形混合 (1−x)φ_l + xφ_v に置き換える。全域で正値で、単相の端
(x=0 で飽和液、x=1 で飽和蒸気) に滑らかにつながる。二相域では単相用の
相関式の重みはほぼ0なので、必要なのは滑らかさと正値性だけである。
"""
function flash_ρe_once(handle, ρ, e, fluid_name::String,
    idx_DU, idx_P, idx_T, idx_Q, idx_V, idx_L, idx_C)
    CoolProp.AbstractState_update(handle, idx_DU, ρ, e)
    P = CoolProp.AbstractState_keyed_output(handle, idx_P)
    T = CoolProp.AbstractState_keyed_output(handle, idx_T)
    μ = CoolProp.AbstractState_keyed_output(handle, idx_V)
    k = CoolProp.AbstractState_keyed_output(handle, idx_L)
    cp = CoolProp.AbstractState_keyed_output(handle, idx_C)

    x = try
        q = CoolProp.AbstractState_keyed_output(handle, idx_Q)
        (0.0 <= q <= 1.0) ? q : NaN
    catch
        NaN
    end

    if isnan(x)
        # 単相: 飽和温度との大小関係で0(過冷液)/1(過熱蒸気)を確定させる
        Tsat_here = CoolProp.PropsSI("T", "P", P, "Q", 0.0, fluid_name)
        x = (T < Tsat_here) ? 0.0 : 1.0
    else
        # 二相: μ,k,cp は意味が定まらないので、飽和液・飽和蒸気の値の線形混合にする
        μ  = (1 - x) * CoolProp.PropsSI("V", "T", T, "Q", 0.0, fluid_name) +
             x * CoolProp.PropsSI("V", "T", T, "Q", 1.0, fluid_name)
        k  = (1 - x) * CoolProp.PropsSI("L", "T", T, "Q", 0.0, fluid_name) +
             x * CoolProp.PropsSI("L", "T", T, "Q", 1.0, fluid_name)
        cp = (1 - x) * CoolProp.PropsSI("C", "T", T, "Q", 0.0, fluid_name) +
             x * CoolProp.PropsSI("C", "T", T, "Q", 1.0, fluid_name)
    end

    return P, T, x, μ, k, cp
end

# -----------------------------------------------------------------------------
# 2. 失敗セルの最近傍補完 (BFS)
# -----------------------------------------------------------------------------

"""
    fill_invalid_cells!(A::Matrix{Float64})

NaNセル (flash失敗箇所) を、最も近い有効セルの値でその場埋めする。
臨界点近傍などで一部の格子点のflashが失敗しても、シミュレーション
実行中にその失敗が伝播しないようにするための処置
(以前議論した「オフライン構築時にtry-catchで検出し、近傍補完する」
という方針の実装)。
"""
function fill_invalid_cells!(A::Matrix{Float64})
    nx, ny = size(A)
    valid = .!isnan.(A)
    all(valid) && return A
    any(valid) || error("すべてのグリッド点でflashが失敗しました。(ρ,h)の範囲設定を見直してください。")

    dist = fill(typemax(Int), nx, ny)
    src  = Array{Tuple{Int,Int}}(undef, nx, ny)
    queue = Vector{Tuple{Int,Int}}()
    sizehint!(queue, nx * ny)

    for i in 1:nx, j in 1:ny
        if valid[i, j]
            dist[i, j] = 0
            src[i, j] = (i, j)
            push!(queue, (i, j))
        end
    end

    head = 1
    dirs = ((1, 0), (-1, 0), (0, 1), (0, -1))
    while head <= length(queue)
        i, j = queue[head]
        head += 1
        for (di, dj) in dirs
            ni, nj = i + di, j + dj
            if 1 <= ni <= nx && 1 <= nj <= ny && dist[ni, nj] > dist[i, j] + 1
                dist[ni, nj] = dist[i, j] + 1
                src[ni, nj] = src[i, j]
                push!(queue, (ni, nj))
            end
        end
    end

    for i in 1:nx, j in 1:ny
        if !valid[i, j]
            si, sj = src[i, j]
            A[i, j] = A[si, sj]
        end
    end
    return A
end

# -----------------------------------------------------------------------------
# 3. (ρ,e) 格子範囲の見積もり
# -----------------------------------------------------------------------------

"""
    estimate_ρe_bounds(fluid_name, T_min, T_max, P_min, P_max; margin=0.15)

想定する運転範囲 (T_min〜T_max, P_min〜P_max) から、(ρ,e) 格子が
カバーすべき範囲を見積もる。(T,P)の矩形領域の四隅に加え、
二相ドーム上の飽和液・飽和蒸気の値も候補に含めることで、
運転中に二相領域へ入っても格子の範囲外に出ないようにしている。
marginは安全マージン (格子の端ぎりぎりで外挿に頼らないため)。
"""
function estimate_ρe_bounds(fluid_name::String, T_min, T_max, P_min, P_max; margin=0.15)
    ρs = Float64[]
    es = Float64[]

    for T in (T_min, T_max), P in (P_min, P_max)
        try
            push!(ρs, CoolProp.PropsSI("Dmass", "T", T, "P", P, fluid_name))
            push!(es, CoolProp.PropsSI("Umass", "T", T, "P", P, fluid_name))
        catch
            # 物理的に無効な(T,P)の組み合わせはスキップ
        end
    end

    Tcrit = CoolProp.PropsSI("Tcrit", "", 0, "", 0, fluid_name)
    Ttriple = CoolProp.PropsSI("Ttriple", "", 0, "", 0, fluid_name)
    T_dome_lo = max(T_min, Ttriple + 1.0)
    T_dome_hi = min(T_max, Tcrit - 1.0)

    if T_dome_hi > T_dome_lo
        for T in range(T_dome_lo, T_dome_hi, length=50), Q in (0.0, 1.0)
            try
                push!(ρs, CoolProp.PropsSI("Dmass", "T", T, "Q", Q, fluid_name))
                push!(es, CoolProp.PropsSI("Umass", "T", T, "Q", Q, fluid_name))
            catch
            end
        end
    end

    isempty(ρs) && error("(ρ,e)範囲の見積もりに失敗しました。T_min/T_max/P_min/P_maxを確認してください。")

    ρ_lo, ρ_hi = extrema(ρs)
    e_lo, e_hi = extrema(es)
    Δρ, Δe = ρ_hi - ρ_lo, e_hi - e_lo

    return (max(ρ_lo - margin * Δρ, 1.0e-3), ρ_hi + margin * Δρ,
            e_lo - margin * Δe, e_hi + margin * Δe)
end

# -----------------------------------------------------------------------------
# 4. 2次元 (ρ,e) テーブルの構築
# -----------------------------------------------------------------------------

"""
    build_ρe_table(fluid_name, ρ_range, e_range; verbose=true)

【診断用】(ρ,e) の矩形格子の全点についてflashを行い、P,T,x,μ,k,cpの
2次元配列を返す (失敗セルは最近傍補完済み)。

本体のテーブル (build_fluid_tables) は (e, δ=ρ−ρ_sl(e)) 座標の
build_eδ_table を使う。この矩形格子版は、圧縮液で精度が足りない
(ファイル冒頭参照) ため、本体では使わず、diagnose_grid.jl や
test_proptables.jl のような診断・動作確認用に残してある。
1つのAbstractStateハンドルを使い回すことで、ハンドル生成コストを
格子点数ぶん払わずに済む。
"""
function build_ρe_table(fluid_name::String, ρ_range, e_range; verbose::Bool=true)
    Nρ, Ne = length(ρ_range), length(e_range)
    P_grid  = fill(NaN, Nρ, Ne)
    T_grid  = fill(NaN, Nρ, Ne)
    x_grid  = fill(NaN, Nρ, Ne)
    μ_grid  = fill(NaN, Nρ, Ne)
    k_grid  = fill(NaN, Nρ, Ne)
    cp_grid = fill(NaN, Nρ, Ne)

    idx_DU = CoolProp.get_input_pair_index("DmassUmass_INPUTS")
    idx_P  = CoolProp.get_param_index("P")
    idx_T  = CoolProp.get_param_index("T")
    idx_Q  = CoolProp.get_param_index("Q")
    idx_V  = CoolProp.get_param_index("V")   # 粘度
    idx_L  = CoolProp.get_param_index("L")   # 熱伝導率 (Lambda)
    idx_C  = CoolProp.get_param_index("C")   # 定圧比熱

    n_fail = 0
    with_abstract_state("HEOS", fluid_name) do handle
        for (i, ρ) in enumerate(ρ_range), (j, e) in enumerate(e_range)
            try
                P, T, x, μ, k, cp = flash_ρe_once(handle, ρ, e, fluid_name,
                    idx_DU, idx_P, idx_T, idx_Q, idx_V, idx_L, idx_C)
                P_grid[i, j]  = P
                T_grid[i, j]  = T
                x_grid[i, j]  = x
                μ_grid[i, j]  = μ
                k_grid[i, j]  = k
                cp_grid[i, j] = cp
            catch
                n_fail += 1
            end
        end
    end

    if verbose
        println("  flash失敗: $(n_fail) / $(Nρ * Ne) グリッド点 (最近傍補完を実施)")
    end

    fill_invalid_cells!(P_grid)
    fill_invalid_cells!(T_grid)
    fill_invalid_cells!(x_grid)
    fill_invalid_cells!(μ_grid)
    fill_invalid_cells!(k_grid)
    fill_invalid_cells!(cp_grid)

    return P_grid, T_grid, x_grid, μ_grid, k_grid, cp_grid
end

"""
    build_eδ_table(fluid_name, e_range, δ_range, ρ_ref; ρ_floor=1e-3, verbose=true)

(e, δ) 格子の全点についてflashを行い、P,T,x,μ,k,cpの2次元配列
(サイズ: length(e_range) × length(δ_range)) を返す。
格子点 (e_i, δ_j) の密度は ρ = ρ_ref(e_i) + δ_j。

  ρ_floor 未満の密度になる点は、物理的に存在しない (飽和液密度から
          大きく下がった領域の端) ので、flashせずに飛ばす。
  flashに失敗した点は、最近傍補完する (fill_invalid_cells!)。
          これらの点は、実際の計算では問い合わせない領域である。

μ,k,cp (粘度・熱伝導率・比熱) は、単相域の強制対流・摩擦係数の計算に
必要な輸送物性で、同じ格子・同じハンドルの流れの中でまとめて取得する。
"""
function build_eδ_table(fluid_name::String, e_range, δ_range, ρ_ref;
    ρ_floor::Float64=1.0e-3, verbose::Bool=true)

    Ne, Nδ = length(e_range), length(δ_range)
    P_grid  = fill(NaN, Ne, Nδ)
    T_grid  = fill(NaN, Ne, Nδ)
    x_grid  = fill(NaN, Ne, Nδ)
    μ_grid  = fill(NaN, Ne, Nδ)
    k_grid  = fill(NaN, Ne, Nδ)
    cp_grid = fill(NaN, Ne, Nδ)

    idx_DU = CoolProp.get_input_pair_index("DmassUmass_INPUTS")
    idx_P  = CoolProp.get_param_index("P")
    idx_T  = CoolProp.get_param_index("T")
    idx_Q  = CoolProp.get_param_index("Q")
    idx_V  = CoolProp.get_param_index("V")   # 粘度
    idx_L  = CoolProp.get_param_index("L")   # 熱伝導率 (Lambda)
    idx_C  = CoolProp.get_param_index("C")   # 定圧比熱

    n_fail = 0
    n_skip = 0
    with_abstract_state("HEOS", fluid_name) do handle
        for (i, e) in enumerate(e_range)
            ρr = ρ_ref(e)
            for (j, δ) in enumerate(δ_range)
                ρ = ρr + δ
                if ρ < ρ_floor
                    n_skip += 1
                    continue
                end
                try
                    P, T, x, μ, k, cp = flash_ρe_once(handle, ρ, e, fluid_name,
                        idx_DU, idx_P, idx_T, idx_Q, idx_V, idx_L, idx_C)
                    P_grid[i, j]  = P
                    T_grid[i, j]  = T
                    x_grid[i, j]  = x
                    μ_grid[i, j]  = μ
                    k_grid[i, j]  = k
                    cp_grid[i, j] = cp
                catch
                    n_fail += 1
                end
            end
        end
    end

    if verbose
        println("  格子 $(Ne)×$(Nδ) = $(Ne * Nδ) 点: ρ<$(ρ_floor) で対象外 $(n_skip) 点, " *
                "flash失敗 $(n_fail) 点 (これらは最近傍補完)")
    end

    fill_invalid_cells!(P_grid)
    fill_invalid_cells!(T_grid)
    fill_invalid_cells!(x_grid)
    fill_invalid_cells!(μ_grid)
    fill_invalid_cells!(k_grid)
    fill_invalid_cells!(cp_grid)

    return P_grid, T_grid, x_grid, μ_grid, k_grid, cp_grid
end

# -----------------------------------------------------------------------------
# 5. 1次元飽和テーブルの構築
# -----------------------------------------------------------------------------

"""
    build_saturation_tables(fluid_name, T_sat_range)

Lockhart-Martinelli型の相関式やCC自然対流のRa計算で必要となる、
飽和液・飽和蒸気の物性をTのみの関数として計算する。
1次元・呼び出し回数も少ないため、高レベルAPI (PropsSI) をそのまま使う。
"""
function build_saturation_tables(fluid_name::String, T_sat_range)
    n = length(T_sat_range)
    μ_l = zeros(n); μ_v = zeros(n)
    ρ_l = zeros(n); ρ_v = zeros(n)
    k_l = zeros(n); k_v = zeros(n)
    cp_l = zeros(n); cp_v = zeros(n)
    h_l = zeros(n); h_v = zeros(n)
    e_l = zeros(n)   # 飽和液の比内部エネルギー (飽和液線の座標 ρ_sl(e) の構築用)
    σ = zeros(n)
    β_l = zeros(n); β_v = zeros(n)
    Psat = zeros(n)

    for (i, T) in enumerate(T_sat_range)
        μ_l[i]  = CoolProp.PropsSI("V", "T", T, "Q", 0.0, fluid_name)
        μ_v[i]  = CoolProp.PropsSI("V", "T", T, "Q", 1.0, fluid_name)
        ρ_l[i]  = CoolProp.PropsSI("D", "T", T, "Q", 0.0, fluid_name)
        ρ_v[i]  = CoolProp.PropsSI("D", "T", T, "Q", 1.0, fluid_name)
        k_l[i]  = CoolProp.PropsSI("L", "T", T, "Q", 0.0, fluid_name)
        k_v[i]  = CoolProp.PropsSI("L", "T", T, "Q", 1.0, fluid_name)
        cp_l[i] = CoolProp.PropsSI("C", "T", T, "Q", 0.0, fluid_name)
        cp_v[i] = CoolProp.PropsSI("C", "T", T, "Q", 1.0, fluid_name)
        h_l[i]  = CoolProp.PropsSI("H", "T", T, "Q", 0.0, fluid_name)
        e_l[i]  = CoolProp.PropsSI("U", "T", T, "Q", 0.0, fluid_name)
        h_v[i]  = CoolProp.PropsSI("H", "T", T, "Q", 1.0, fluid_name)
        σ[i]    = CoolProp.PropsSI("SURFACE_TENSION", "T", T, "Q", 0.0, fluid_name)
        β_l[i]  = CoolProp.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T, "Q", 0.0, fluid_name)
        β_v[i]  = CoolProp.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T, "Q", 1.0, fluid_name)
        Psat[i] = CoolProp.PropsSI("P", "T", T, "Q", 0.0, fluid_name)
    end

    return (; μ_l, μ_v, ρ_l, ρ_v, k_l, k_v, cp_l, cp_v, h_l, h_v, e_l, σ, β_l, β_v, Psat)
end

# -----------------------------------------------------------------------------
# 6. 全体の組み立て
# -----------------------------------------------------------------------------

"""
    build_raw_tables(fluid_name; T_min, T_max, P_min=1e5, P_max,
                     Ne=400, Δδ=1.0, δ_liq_max=20.0, N_T_sat=300,
                     margin=0.15, verbose=true) -> NamedTuple

テーブルの素のデータ (配列・範囲・数値) を計算する。時間がかかる (CoolProp を数十万回
呼ぶ) のは、この部分。戻り値は、Float64, Vector, Matrix, 範囲, String だけで、
モジュールに依存する型を含まない (キャッシュに、そのまま保存できる)。補間オブジェクトと
FluidTables への組み立ては、assemble_tables が行う。

  Ne        … e方向の格子数 (e_range は運転範囲から見積もる)
  Δδ        … δ=ρ−ρ_sl(e) 方向の格子間隔 [kg/m³]。δ=0 (飽和液線) は必ず
               格子点になる。圧縮液は P がほぼ δ に線形なので粗くてよいが、
               粗くするほど液の ∂P/∂ρ の誤差が増える (Δδ=1:±1%以内,
               2:±1.2%, 4:±2.5%。Pythonで確認)。
  δ_liq_max … 圧縮液の側 (δ>0) の上限 [kg/m³]。P_max までの圧縮液で
               δ は数 kg/m³ 程度。余裕を持たせて既定 20。
  N_T_sat   … 飽和テーブルの温度方向の点数 (飽和液線 ρ_sl(e) の精度にも効く)

構築の順序: 飽和テーブル → 飽和液線 ρ_sl(e) → 2次元テーブル。
飽和液線は、2次元テーブルの座標の基準 (δ=0 の線) になる。
"""
function build_raw_tables(fluid_name::String;
    T_min::Float64, T_max::Float64,
    P_min::Float64=1.0e5, P_max::Float64,
    Ne::Int=400, Δδ::Float64=1.0, δ_liq_max::Float64=20.0, N_T_sat::Int=300,
    margin::Float64=0.15, verbose::Bool=true)

    Tcrit   = CoolProp.PropsSI("Tcrit", "", 0, "", 0, fluid_name)
    Ttriple = CoolProp.PropsSI("Ttriple", "", 0, "", 0, fluid_name)
    Pcrit   = CoolProp.PropsSI("Pcrit", "", 0, "", 0, fluid_name)

    # --- 1) 飽和テーブル ---
    T_sat_lo = max(T_min, Ttriple + 1.0)
    T_sat_hi = min(T_max, Tcrit - 1.0)
    T_sat_lo >= T_sat_hi && error("T_min/T_maxと流体の三重点・臨界点の関係から、" *
                                   "飽和テーブル用の有効なT範囲が確保できません。")
    T_sat_range = range(T_sat_lo, T_sat_hi, length=N_T_sat)

    verbose && println("[$fluid_name] 1次元飽和テーブル構築中...")
    sat = build_saturation_tables(fluid_name, T_sat_range)

    # --- 2) 飽和液線 ρ_sl(e): 2次元テーブルの座標の基準 ---
    # e_l(T) は T とともに単調に増える。Gridded補間には狭義単調増加の格子が必要。
    issorted(sat.e_l; lt=<=) ||
        error("飽和液の内部エネルギー e_l(T) が単調増加ではありません。T範囲を確認してください。")
    ρ_ref = make_interp1_nonuniform(sat.e_l, sat.ρ_l)

    # --- 3) 2次元テーブルの格子 ---
    verbose && println("[$fluid_name] 格子の範囲を見積もり中...")
    ρ_lo, ρ_hi, e_lo, e_hi = estimate_ρe_bounds(fluid_name, T_min, T_max, P_min, P_max; margin=margin)
    e_range = range(e_lo, e_hi, length=Ne)

    # δ の格子: δ=0 を必ず格子点にする。負側は ρ≈ρ_lo (ほぼ真空) まで。
    ρ_ref_max = maximum(ρ_ref.(e_range))
    n_neg = ceil(Int, (ρ_ref_max - ρ_lo) / Δδ)
    n_pos = round(Int, δ_liq_max / Δδ)
    δ_range = range(-n_neg * Δδ, step=Δδ, length=n_neg + n_pos + 1)

    verbose && println("[$fluid_name] 2次元 (e,δ) テーブル構築中 " *
                       "(e: $(Ne) 点 × δ: $(length(δ_range)) 点)...")
    P_grid, T_grid, x_grid, μ_grid, k_grid, cp_grid =
        build_eδ_table(fluid_name, e_range, δ_range, ρ_ref; verbose=verbose)

    return (fluid_name = fluid_name, Pcrit = Pcrit, Tcrit = Tcrit,
            e_range = e_range, δ_range = δ_range, T_sat_range = T_sat_range,
            P_grid = P_grid, T_grid = T_grid, x_grid = x_grid,
            μ_grid = μ_grid, k_grid = k_grid, cp_grid = cp_grid,
            sat = sat)
end

# --- 補間オブジェクトの組み立てヘルパー ---
# 等間隔格子 (e_range, δ_range, T_sat_range はすべて range なので scale+BSpline が使える。
# Tsat_P と ρ_ref は格子が非等間隔 (Psat や e_sl の値が格子) なので Gridded を使う。
make_interp1(grid, vals) = extrapolate(scale(interpolate(vals, BSpline(Linear())), grid), Line())
make_interp2(g1, g2, vals) = extrapolate(scale(interpolate(vals, (BSpline(Linear()), BSpline(Linear()))), g1, g2), Line())
make_interp1_nonuniform(xs, vals) = extrapolate(interpolate((xs,), vals, Gridded(Linear())), Line())

"""
    assemble_tables(raw) -> FluidTables

素のデータ (build_raw_tables の戻り値、またはキャッシュから読んだもの) から、補間
オブジェクトと FluidTables を組み立てる。一瞬で終わる。戻り値は FLUID_TABLES にも登録する。
"""
function assemble_tables(raw)
    sat = raw.sat
    T_sat_range = raw.T_sat_range
    e_range, δ_range = raw.e_range, raw.δ_range
    ρ_ref = make_interp1_nonuniform(sat.e_l, sat.ρ_l)

    # FluidTables は、同じ役割のフィールドを、同じ型パラメータで束ねている (I2: 2次元の6個、
    # I1: 1次元の14個、IG: Gridded の2個)。型が揃っていないと、コンストラクタが、分かりにくい
    # MethodError になるので、先に調べて、分かりやすいエラーにする。
    I1s = unique(typeof.([make_interp1(T_sat_range, getfield(sat, k)) for k in
                          (:Psat, :μ_l, :μ_v, :ρ_l, :ρ_v, :k_l, :k_v, :cp_l, :cp_v, :h_l, :h_v, :σ, :β_l, :β_v)]))
    length(I1s) == 1 ||
        error("FluidTables: 1次元の補間の型が揃っていません (型パラメータ I1)。proptables.jl の FluidTables の定義と、" *
              "make_interp1 の戻り値の型を見直してください: $(I1s)")
    IGs = unique(typeof.([ρ_ref, make_interp1_nonuniform(sat.Psat, collect(T_sat_range))]))
    length(IGs) == 1 ||
        error("FluidTables: Gridded の補間の型が揃っていません (型パラメータ IG): $(IGs)")
    I2s = unique(typeof.([RhoEInterp(make_interp2(e_range, δ_range, g), ρ_ref) for g in
                          (raw.P_grid, raw.T_grid, raw.x_grid, raw.μ_grid, raw.k_grid, raw.cp_grid)]))
    length(I2s) == 1 ||
        error("FluidTables: 2次元の補間の型が揃っていません (型パラメータ I2): $(I2s)")

    tables = FluidTables(
        raw.fluid_name, raw.Pcrit, raw.Tcrit,
        e_range, δ_range, ρ_ref,
        RhoEInterp(make_interp2(e_range, δ_range, raw.P_grid), ρ_ref),
        RhoEInterp(make_interp2(e_range, δ_range, raw.T_grid), ρ_ref),
        RhoEInterp(make_interp2(e_range, δ_range, raw.x_grid), ρ_ref),
        RhoEInterp(make_interp2(e_range, δ_range, raw.μ_grid), ρ_ref),
        RhoEInterp(make_interp2(e_range, δ_range, raw.k_grid), ρ_ref),
        RhoEInterp(make_interp2(e_range, δ_range, raw.cp_grid), ρ_ref),
        T_sat_range,
        make_interp1(T_sat_range, sat.Psat),
        make_interp1_nonuniform(sat.Psat, collect(T_sat_range)),
        make_interp1(T_sat_range, sat.μ_l), make_interp1(T_sat_range, sat.μ_v),
        make_interp1(T_sat_range, sat.ρ_l), make_interp1(T_sat_range, sat.ρ_v),
        make_interp1(T_sat_range, sat.k_l), make_interp1(T_sat_range, sat.k_v),
        make_interp1(T_sat_range, sat.cp_l), make_interp1(T_sat_range, sat.cp_v),
        make_interp1(T_sat_range, sat.h_l), make_interp1(T_sat_range, sat.h_v),
        make_interp1(T_sat_range, sat.σ),
        make_interp1(T_sat_range, sat.β_l),
        make_interp1(T_sat_range, sat.β_v),
    )
    FLUID_TABLES[raw.fluid_name] = tables
    return tables
end

"""
    build_fluid_tables(fluid_name; kwargs...) -> FluidTables

build_raw_tables (素のデータの計算) と assemble_tables (組み立て) をまとめて行う。
引数は build_raw_tables と同じ。
"""
function build_fluid_tables(fluid_name::String; verbose::Bool=true, kwargs...)
    raw = build_raw_tables(fluid_name; verbose=verbose, kwargs...)
    tables = assemble_tables(raw)
    verbose && println("[$fluid_name] テーブル構築完了。")
    return tables
end

# -----------------------------------------------------------------------------
# 7. キャッシュ (Serialization) とエントリポイント
# -----------------------------------------------------------------------------

"""
    save_raw_tables(raw, path)

素のデータ (build_raw_tables の戻り値) をファイルに保存する。FluidTables は保存しない
(型がモジュールに依存するため。ファイル冒頭の TABLE_VERSION の説明参照)。
"""
function save_raw_tables(raw, path::String)
    mkpath(dirname(path))
    open(path, "w") do io
        serialize(io, raw)
    end
end

"""
    load_fluid_tables!(path) -> FluidTables

保存済みの素のデータを読み込み、FluidTables を組み立てて、FLUID_TABLES にも登録して返す。
"""
function load_fluid_tables!(path::String)
    raw = open(path, "r") do io
        deserialize(io)
    end
    return assemble_tables(raw)
end

"""
    init_fluid_lookups(fluid_name; T_min=230.0, T_max=360.0, P_min=1e5, P_max=4.5e6,
                        Ne=400, Δδ=1.0, δ_liq_max=20.0, N_T_sat=300,
                        cache_dir=..., force_rebuild=false, verbose=true)

既存コードの `init_fluid_lookups` と互換なエントリポイント。
キャッシュファイルが存在すればロードするだけで済ませ、なければ
構築してからキャッシュに保存する。
"""
function init_fluid_lookups(fluid_name::String;
    T_min::Float64=230.0, T_max::Float64=360.0,
    P_min::Float64=1.0e5, P_max::Float64=4.5e6,
    Ne::Int=400, Δδ::Float64=1.0, δ_liq_max::Float64=20.0, N_T_sat::Int=300,
    cache_dir::String=joinpath(@__DIR__, "..", "fluid_tables_cache"),
    force_rebuild::Bool=false, verbose::Bool=true)

    cache_path = joinpath(cache_dir,
        "$(TABLE_VERSION)_$(fluid_name)_T$(T_min)-$(T_max)_P$(P_min)-$(P_max)_Ne$(Ne)_dd$(Δδ)_dl$(δ_liq_max)_NT$(N_T_sat).jls")

    if !force_rebuild && isfile(cache_path)
        verbose && println("[$fluid_name] キャッシュからテーブルを読み込みます: $cache_path")
        return load_fluid_tables!(cache_path)
    end

    raw = build_raw_tables(fluid_name;
        T_min=T_min, T_max=T_max, P_min=P_min, P_max=P_max,
        Ne=Ne, Δδ=Δδ, δ_liq_max=δ_liq_max, N_T_sat=N_T_sat, verbose=verbose)

    save_raw_tables(raw, cache_path)
    verbose && println("[$fluid_name] キャッシュを保存しました: $cache_path")
    tables = assemble_tables(raw)
    verbose && println("[$fluid_name] テーブル構築完了。")
    return tables
end