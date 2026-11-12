# =============================================================================
# ic_Nsweep.jl … CC の液の層数 N_CC_LAYER = 1, 2, 4, 8, 16 を、C2 の条件の :short (τ = 90 s まで) で走らせる (2026-11-05。未実行)
#   目的: 起動直後の圧力の山・界面の層 1 の温度 (T_l1 − T_lN)・τ = 100 s の T_cc が、層数でどう変わるかを見る
#         (C2 と C2lump の差の経路 A: 凝縮の潜熱が層 1 に集中すること、の確認)。
#   層の厚さは A: 等比 r = CC_LAYER_RATIO = 1.23 を固定 (層数だけ変える。層 1 の厚さは N で変わる: N=16 で約 0.08 mm、N=8 で約 0.4 mm など)。
#   使い方: プロジェクトのルート (src/、scripts/、data/、Project.toml がある所) で
#       julia --project=. scripts/ic_Nsweep.jl
#   しくみ: N_CC_LAYER は src/types.jl の const (状態の長さを決める) なので、1 回の Julia の中では変えられない。
#           N ごとにプロジェクトを <ルート>_Nsweep/N<n>/ にコピーし、コピー側の src/types.jl だけ書き換えて、別の Julia のプロセスで走らせる。
#           元の src/ は変更しない。結果は、元の res/ に transient_flat_ammonia_ic_C2_N<n>/ として戻す (同名があれば飛ばす = 上書きしない)。
#   段階: (1) 各コピーを 1 つずつ事前コンパイル → (2) N_JOBS 本まで並列に計算 → (3) 結果を元の res/ に戻し、要約を表示。
#   注意: transient_flat.jl は SWEEP の dT_onset、n_vap_vl_init、f_l_grv_init を受け付ける版 (11/3 版) であること (下で調べる)。
# =============================================================================

using Pkg; Pkg.activate(".")
using Printf, Statistics
using CSV, DataFrames      # 要約の表示用 (プロジェクトの環境に入っている)

# ---------------------------- 設定 ----------------------------
NS        = [1, 2, 4, 8, 16]    # 層数
N_JOBS    = 3                   # 同時に走らせる数 (CPU とメモリに合わせて変える)
T_END     = 100.0               # 計算の終わり [s] (加熱開始は 10 s なので τ = 90 s)
TAG_BASE  = "_ic_C2_N"          # 結果フォルダ: res/transient_flat_ammonia<tag>/ (tag = TAG_BASE * N)
CASE      = (dT_onset = 0.55, f_l_grv_init = 0.05, n_vap_vl_init = 11)   # ic_cases.jl の C2 と同じ
KEEP_COPIES = true              # false なら、終わったあとコピー (<ルート>_Nsweep/) を消す
# ----------------------------------------------------------------

root = abspath(pwd())
isfile(joinpath(root, "src", "types.jl")) || error("ルート (src/types.jl がある所) で実行してください。いまの場所: ", root)
res  = joinpath(root, "res")
work = root * "_Nsweep"
tagof(n) = TAG_BASE * string(n)

# transient_flat.jl が SWEEP のキーを受け付けるかを確かめる (受け付けない版だと値が黙って無視される)
tf = read(joinpath(root, "scripts", "transient_flat.jl"), String)
for key in ("dT_onset", "n_vap_vl_init", "f_l_grv_init")
    occursin(":" * key, tf) || error("transient_flat.jl が SWEEP の ", key, " を受け付けません (11/3 版が必要です)")
end

result_dir(n) = joinpath(res, "transient_flat_ammonia" * tagof(n))
todo = [n for n in NS if !isfile(joinpath(result_dir(n), "summary.csv"))]
for n in NS
    n in todo || println("(飛ばす: 保存済み) ", result_dir(n))
end
isempty(todo) && (println("すべて保存済みです。やり直すときは、結果フォルダを別の場所へ移してください (上書きしません)"); exit(0))

# ---------- (1) コピーして types.jl を書き換え、事前コンパイル ----------
SKIP = Set(["res", ".git", basename(work), "Nsweep_logs"])
julia = Base.julia_cmd()
function make_copy(n)
    d = joinpath(work, "N" * string(n))
    isdir(d) && rm(d; recursive = true, force = true)
    mkpath(d)
    for e in readdir(root)
        e in SKIP && continue
        cp(joinpath(root, e), joinpath(d, e); force = true)
    end
    mkpath(joinpath(d, "res"))
    p = joinpath(d, "src", "types.jl")
    s = read(p, String)
    s2 = replace(s, r"^const N_CC_LAYER[^\r\n]*"m => "const N_CC_LAYER = $(n)   # ic_Nsweep.jl が書き換えた")
    s2 = replace(s2, r"^const CC_LUMPED[^\r\n]*"m => "const CC_LUMPED = false   # ic_Nsweep.jl が書き換えた (層数は N_CC_LAYER で決める)")
    s2 == s && error("types.jl の書き換えが効いていません (const N_CC_LAYER の行が見つかりません)")
    occursin("const N_CC_LAYER = $(n) ", s2) || error("N_CC_LAYER = ", n, " の書き換えに失敗しました")
    write(p, s2)
    runner = """
    using Pkg; Pkg.activate(".")
    push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
    using LHPTransientModel
    get(ENV, "NSWEEP_PRECOMPILE", "") == "1" && (println("事前コンパイル済み"); exit(0))
    LHPTransientModel.N_CC_LAYER == $(n) || error("N_CC_LAYER が $(n) になっていません: ", LHPTransientModel.N_CC_LAYER)
    println("N_CC_LAYER = ", LHPTransientModel.N_CC_LAYER, "  CC_LUMPED = ", LHPTransientModel.CC_LUMPED)
    global SWEEP = (stage = :short, t_end = $(T_END), dT_onset = $(CASE.dT_onset), f_l_grv_init = $(CASE.f_l_grv_init),
                    n_vap_vl_init = $(CASE.n_vap_vl_init), tag = "$(tagof(n))")
    include(joinpath(@__DIR__, "transient_flat.jl"))
    open(joinpath(@__DIR__, "..", "done.txt"), "w") do io
        println(io, sol.retcode, " ", sol.t[end])
    end
    """
    write(joinpath(d, "scripts", "_run_Nsweep.jl"), runner)
    return d
end

println("コピー先: ", work)
dirs = Dict{Int,String}()
for n in todo
    dirs[n] = make_copy(n)
    println("N = ", n, ": コピーして事前コンパイル中 ...")
    run(addenv(Cmd(`$julia --project=. scripts/_run_Nsweep.jl`; dir = dirs[n]), "NSWEEP_PRECOMPILE" => "1"))   # 本番と同じ読み込み経路で 1 つずつ事前コンパイル
end

# ---------- (2) 並列に計算 ----------
logdir = joinpath(res, "Nsweep_logs"); mkpath(logdir)
sem = Base.Semaphore(N_JOBS)
status = Dict{Int,Any}()
@sync for n in todo
    @async begin
        Base.acquire(sem)
        try
            logf = joinpath(logdir, "ic_C2_N$(n).log")
            println("N = ", n, ": 計算開始 (ログ: ", logf, ")")
            t0 = time()
            ok = try
                run(pipeline(Cmd(`$julia --project=. scripts/_run_Nsweep.jl`; dir = dirs[n]); stdout = logf, stderr = logf, append = false))
                true
            catch e
                println("N = ", n, ": 失敗: ", sprint(showerror, e)); false
            end
            status[n] = (ok = ok, minutes = (time() - t0) / 60)
            println("N = ", n, ": ", ok ? "完了" : "失敗", " (", round((time() - t0) / 60, digits = 1), " 分)")
        finally
            Base.release(sem)
        end
    end
end

# ---------- (3) 結果を元の res/ に戻す ----------
for n in todo
    haskey(status, n) && status[n].ok || continue
    src = joinpath(dirs[n], "res")
    found = [d for d in readdir(src) if startswith(d, "transient_flat_ammonia") && endswith(d, tagof(n))]
    if isempty(found)
        println("N = ", n, ": 結果フォルダが見つかりません (", src, ")"); continue
    end
    for d in found
        dst = joinpath(res, d)
        if ispath(dst)
            println("N = ", n, ": ", dst, " があるので戻しません (コピー側に残してあります: ", joinpath(src, d), ")")
        else
            cp(joinpath(src, d), dst)
            println("N = ", n, ": 戻しました → ", dst)
        end
    end
end

# ---------- 要約 ----------
try
    println("\n  N | 層1の厚さ[mm] | P_cc の山 [kPa] (τ[s]) | max(T_l1−T_lN) [K] (τ[s]) | ΔT_cc(τ=100 s) [K] | max(Q_nc_l) [mW]")
    r = 1.23; hliq = 9.0
    for n in NS
        f = joinpath(result_dir(n), "summary.csv")
        isfile(f) || continue
        df = CSV.read(f, DataFrame)
        t = df.t .- 10.0
        pre = t .< 0
        any(pre) || (pre = t .< 0.5)
        dP = (df.P_cc .- mean(df.P_cc[pre])) ./ 1e3
        m = (t .>= 15) .& (t .<= 100)
        ip = argmax(ifelse.(m, dP, -Inf))
        dTs = df.T_l1 .- df.T_lN
        it = argmax(ifelse.(m, dTs, -Inf))
        i100 = argmin(abs.(t .- 100))
        dTcc = df.T_cc[i100] - mean(df.T_cc[pre])
        g1 = n == 1 ? 1.0 : (r - 1) / (r^n - 1)
        @printf("%3d | %8.3f | %8.2f (%5.1f) | %8.3f (%5.1f) | %8.3f | %8.1f\n", n, g1 * hliq, dP[ip], t[ip], dTs[it], t[it], dTcc, 1e3 * maximum(df.Q_nc_l[m]))
    end
    println("(層 1 の厚さは h_liq ≈ 9 mm から見積もった値。ΔT_cc は加熱前の平均からの上昇)")
catch e
    println("要約の表示に失敗しました (結果の CSV は res/ にあります): ", sprint(showerror, e))
end

KEEP_COPIES || rm(work; recursive = true, force = true)
println("\n完了。結果: res/transient_flat_ammonia", TAG_BASE, "<N>/ (summary.csv、profiles.csv、budget.csv)。ログ: res/Nsweep_logs/")