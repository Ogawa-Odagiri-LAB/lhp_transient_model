module LHPTransientModel

# ループヒートパイプ (LHP) の過渡解析モデル。
# 完全陰的な DAE ソルバー IDA (Sundials.jl) で、残差 F(t, u, du) = 0 を時間積分する。
# 作動流体の物性は、CoolProp から事前にテーブル化して、実行時は補間のみ行う。

using CoolProp, ComponentArrays, Interpolations, LinearAlgebra, Parameters, Random,
      Serialization, SparseArrays, StaticArrays
using CSV, DataFrames
using Sundials
using DiffEqBase: DAEProblem, DAEFunction, solve

export LHPGeometry, WickGeometry, GrvGeometry, TLGeometry, CondGeometry, CCGeometry
export FlatLHPGeometry, FlatWickGeometry, SecondaryWickGeometry, FlatGrvGeometry, FlatCCGeometry
export LHPTransientCondition
export init_fluid_lookups, FluidTables
export solve_transient_lhp!, components_at, applied_heat
export energy_balance
export save_results, results_summary, results_profiles, results_wick, results_geometry

# プロット (plot.jl) は、計算とは独立の単体のファイルで、このモジュールには含めない。
# 保存した CSV を、plot_results.jl で、いつでも眺められる。

# ファイルを役割ごとに分割して読み込む (依存の順)
include(joinpath(@__DIR__, "smoothing.jl"))            # Level0: 平滑化ユーティリティ
include(joinpath(@__DIR__, "types.jl"))                # 構造体・状態ベクトル
include(joinpath(@__DIR__, "physics.jl"))              # SUS304 の物性
include(joinpath(@__DIR__, "proptables.jl"))           # 流体物性のテーブル化
include(joinpath(@__DIR__, "propaccess.jl"))           # Level2: 物性へのアクセス
include(joinpath(@__DIR__, "correlations.jl"))         # Level1: 相関式
include(joinpath(@__DIR__, "heattransfer.jl"))         # Level3/4: 熱伝達・摩擦
include(joinpath(@__DIR__, "geometry_cc.jl"))          # CC の気液分離幾何
include(joinpath(@__DIR__, "geometry_flatcc.jl"))      # 平板型CC (半円筒) の気液分離幾何
include(joinpath(@__DIR__, "residuals_transport.jl"))  # 輸送管の流体
include(joinpath(@__DIR__, "residuals_wall.jl"))       # 輸送管の壁
include(joinpath(@__DIR__, "residuals_wick.jl"))       # ウィック
include(joinpath(@__DIR__, "residuals_cc.jl"))         # CC
include(joinpath(@__DIR__, "residuals_evapcase.jl"))   # 蒸発器ケース
include(joinpath(@__DIR__, "residuals_grv.jl"))        # グルーブ
include(joinpath(@__DIR__, "residuals_main.jl"))       # 全体の残差 (円筒型)
include(joinpath(@__DIR__, "residuals_flatwick.jl"))   # 平板型: primary / secondary wick
include(joinpath(@__DIR__, "residuals_flatcc.jl"))     # 平板型: CC
include(joinpath(@__DIR__, "residuals_flatevapcase.jl")) # 平板型: 蒸発器ケース
include(joinpath(@__DIR__, "residuals_flatmain.jl"))   # 平板型: 全体の残差
include(joinpath(@__DIR__, "initialcondition.jl"))     # 整合初期条件
include(joinpath(@__DIR__, "sparsity.jl"))             # ヤコビアンの疎構造
include(joinpath(@__DIR__, "solverfuncs.jl"))          # IDA による時間積分
include(joinpath(@__DIR__, "postprocess.jl"))          # 後処理 (エネルギー収支)
include(joinpath(@__DIR__, "savefuncs.jl"))            # 結果の保存 (CSV)

end