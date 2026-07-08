# lhp_transient_model
LHPの過渡モデルです．流体は陰解法，管壁は陽解法で解いています．ソルバーはJuliaのNonlinearSolve.jlを用いていて，アルゴリズムは信頼区間法を採用しています．
初期条件として，定常解析の結果を用いています．
## 使い方
scriptsの中のdebug.jl，plots.jl，steady.jl，transient.jlのいずれかを実行します．
- debug.jl
debug.jlはその名の通り作った関数のデバッグを行うための場所で，いろいろ実験しているのでぐちゃぐちゃになっています．
- plots.jl
プロット関数がまとめられています．定常解析も過渡解析もここにある関数でプロットできます．
- steady.jl
定常解析を行います．
- transient.jl
過渡解析を行います．


