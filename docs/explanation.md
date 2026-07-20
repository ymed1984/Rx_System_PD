# OMA-to-BER実装と実電流スケールにおける最適しきい値

## 1. この資料の目的

この資料は、NRZ/OOK受信器の2つの電流レベルをガウス分布で近似したときの、
単一判定しきい値とBERの関係を説明するものです。特に、µAやnAの実電流値を
そのまま数値最適化した場合に生じる精度問題と、本コードベースで採用している
無次元化した解析的なしきい値計算を扱います。

対象実装は `src/oma_ber/ber.py` の `optimum_threshold()` と、
`src/oma_ber/time_domain/statistical.py` のGaussian-mixture最適化です。

前半は、単一の電流比較器によるNRZ/OOK判定を想定した簡略化されたスカラー・
ガウスモデルです。後半では決定論的ISIと解析的ガウス雑音をTIA電圧のGaussian
mixtureへ拡張します。クロックジッタ、非ガウス雑音、CDRは含みません。

OMAから光レベル、PD電流、各雑音成分へ至るシステム全体の説明は、
[OMA-to-BER／フォトダイオード受信系の物理モデル](MechanismExplanation.md)を
参照してください。本資料は、その中の単一しきい値最適化を詳しく導出します。

## 2. 受信電流の確率モデル

送信ビット0と1に対応するサンプル電流を、次の正規分布として扱います。

```text
bit 0: I ~ N(mu0, sigma0^2)
bit 1: I ~ N(mu1, sigma1^2)
```

各変数の単位は次のとおりです。

- `mu0_a`: ビット0の平均電流 [A]
- `mu1_a`: ビット1の平均電流 [A]
- `sigma0_a`: ビット0のRMS電流雑音 [A]
- `sigma1_a`: ビット1のRMS電流雑音 [A]
- `threshold_a`: 判定しきい値 [A]

実装では次の物理条件を要求します。

```text
mu1 > mu0
sigma0 > 0
sigma1 > 0
```

`sigma0` と `sigma1` が異なるのは、ショット雑音やRIN雑音が光レベルに依存する
ためです。例えば、本コードベースの基本雑音モデルはレベル `i` ごとに次式です。

```text
sigma_shot_i = sqrt(2 q (I_i + F_dark I_dark) Bn)
sigma_tia    = i_n sqrt(Bn)
sigma_rin_i  = Rpd P_i sqrt(RIN_linear Bn)

sigma_i = sqrt(
    sigma_shot_i^2
    + sigma_tia^2
    + sigma_rin_i^2
    + sigma_thermal^2
)
```

ここで、`q` は電気素量、`Bn` は等価雑音帯域 [Hz]、`Rpd` はPDの受光感度
[A/W]、`F_dark` は暗電流Fano係数です。既定値 `F_dark=1` は理想Poisson雑音を
表します。信号電流が大きいビット1側では、一般にショット雑音とRIN雑音も大きく
なるため、`sigma1 > sigma0` となる場合があります。

## 3. 単一しきい値による判定

判定規則は次のとおりです。

```text
I < gamma  -> bit 0
I >= gamma -> bit 1
```

`gamma` は判定しきい値 [A]です。

ビット0送信時には、電流がしきい値を上回ると誤りです。

```text
P(error | 0) = Q((gamma - mu0) / sigma0)
```

ビット1送信時には、電流がしきい値を下回ると誤りです。

```text
P(error | 1) = Q((mu1 - gamma) / sigma1)
```

ガウスQ関数は次式です。

```text
Q(x) = 1 / sqrt(2 pi) * integral[x, infinity] exp(-u^2 / 2) du
     = 0.5 * erfc(x / sqrt(2))
```

0と1の出現確率が等しいと仮定すると、BERは次式になります。

```text
BER(gamma) = 0.5 * [
    Q((gamma - mu0) / sigma0)
    + Q((mu1 - gamma) / sigma1)
]
```

この式の符号方向は重要です。`gamma` を高くするとビット0の誤りは減少しますが、
ビット1の誤りは増加します。最適しきい値は、このトレードオフを最小化します。

## 4. 等分散の場合

`sigma0 = sigma1` の場合、2つの正規分布は同じ形状です。等確率ビットに対する
最適しきい値は平均電流の中点になります。

```text
gamma_opt = (mu0 + mu1) / 2
```

この場合、両側の正規化距離は等しくなります。

```text
(gamma_opt - mu0) / sigma0
=
(mu1 - gamma_opt) / sigma1
```

## 5. 不等分散の場合の停留条件

`sigma0 != sigma1` の場合、中点は一般にBER最小点ではありません。

標準正規分布の確率密度関数を次のように定義します。

```text
phi(x) = 1 / sqrt(2 pi) * exp(-x^2 / 2)
```

Q関数の微分は次式です。

```text
dQ(x) / dx = -phi(x)
```

BERを `gamma` で微分すると、次式が得られます。

```text
dBER / dgamma = 0.5 * [
    -phi((gamma - mu0) / sigma0) / sigma0
    +phi((mu1 - gamma) / sigma1) / sigma1
]
```

停留点では微分が0なので、次の密度交点条件が成立します。

```text
phi((gamma - mu0) / sigma0) / sigma0
=
phi((mu1 - gamma) / sigma1) / sigma1
```

指数関数を展開して対数を取ると、次式になります。

```text
((gamma - mu0) / sigma0)^2
-
((mu1 - gamma) / sigma1)^2
=
2 ln(sigma1 / sigma0)
```

これは `gamma` に関する二次方程式です。

## 6. A単位で直接最適化する問題

以前の実装では、`gamma` をA単位のまま有界数値最適化していました。この方法で
絶対許容誤差が例えば `1e-5` なら、許容誤差は10 µAに相当します。

一方、Siフォトニクス受信器では次のような電流スケールが普通に現れます。

```text
mu0 = 0 A
mu1 = 1e-6 A
```

このとき探索区間全体が1 µAしかないため、10 µAの絶対許容誤差は探索区間より
大きくなります。その結果、最適化ルーチンがBER最小点へ十分に収束する前に終了
する可能性があります。

この問題は物理式の誤りではなく、次元を持つ値に固定の絶対許容誤差を適用する
ことによる数値スケールの問題です。

## 7. しきい値の無次元化

電流レベル間隔を次のように定義します。

```text
Delta_mu = mu1 - mu0
```

しきい値を無次元変数 `t` で表します。

```text
t = (gamma - mu0) / Delta_mu
gamma = mu0 + t Delta_mu
```

`gamma` を2つの平均電流の間に制限すると、探索範囲は電流スケールによらず
常に次の範囲になります。

```text
0 <= t <= 1
```

さらに、次の無次元量を定義します。

```text
d0 = Delta_mu / sigma0
d1 = Delta_mu / sigma1
```

密度交点条件へ代入すると、次の二次方程式になります。

```text
A t^2 + B t + C = 0
```

係数は次のとおりです。

```text
A = d0^2 - d1^2
B = 2 d1^2
C = -d1^2 - 2 ln(sigma1 / sigma0)
```

この形式では、すべての係数が無次元です。`mu` と `sigma` を同じ倍率で変更しても
`d0`、`d1`、`t` は変化しません。そのため、nA、µA、mAのどの電流スケールでも
同じ正規化しきい値が得られます。

## 8. 二次方程式の数値安定化

### 8.1 係数の正規化

二次方程式の係数を次の代表値で割ります。

```text
scale = max(abs(A), abs(B), abs(C))
```

```text
An = A / scale
Bn = B / scale
Cn = C / scale
```

これにより、二次係数が実質的に0かどうかを、電流単位や係数の絶対的な大きさに
依存せず判定できます。`An` が浮動小数点精度に対して十分小さい場合は、退化した
一次方程式として扱います。

```text
t = -Cn / Bn
```

### 8.2 安定した二次方程式の解法

通常の解の公式は次式です。

```text
t = (-Bn +/- sqrt(D)) / (2 An)
D = Bn^2 - 4 An Cn
```

`-Bn` と `sqrt(D)` が近い場合、減算による桁落ちが起きる可能性があります。
そのため実装では、次の安定化された形を使用します。

```text
q = -0.5 * (Bn + sign(Bn) * sqrt(D))

t1 = q / An
t2 = Cn / q
```

丸め誤差で判別式がごく小さい負値になる場合は、浮動小数点精度の範囲内で0へ
丸めます。明確に負の場合は実数の密度交点がないものとして扱います。

## 9. 最小BER候補の選択

二次方程式の実根は、必ずしも単一しきい値BERの最小点とは限りません。
また、分散比が極端な場合には区間端が最良となる可能性もあります。

そのため実装では次を候補にします。

```text
t = 0
t = 1
0 <= t <= 1 を満たす二次方程式の全実根
```

各候補をA単位へ戻します。

```text
gamma_candidate = mu0 + t Delta_mu
```

その後、元のBER式で各候補を評価し、BERが最小の候補を返します。

```text
gamma_opt = argmin BER(gamma_candidate)
```

密度交点を求める式と、最終的な目的関数であるBER式を分けているため、複数根や
端点を含む場合でも結果を明示的に検証できます。

## 10. 数値例

次の実電流条件を考えます。

```text
mu0    = 0 A
mu1    = 1.0e-6 A
sigma0 = 1.0e-7 A
sigma1 = 2.0e-7 A
```

レベル間隔は次のとおりです。

```text
Delta_mu = 1.0e-6 A
```

無次元量は次のようになります。

```text
d0 = Delta_mu / sigma0 = 10
d1 = Delta_mu / sigma1 = 5
```

したがって二次方程式の係数は次のとおりです。

```text
A = 10^2 - 5^2 = 75
B = 2 * 5^2 = 50
C = -5^2 - 2 ln(2)
```

区間内のBER最小根は次の値です。

```text
t_opt     = 0.3470550625549...
gamma_opt = 3.470550625549e-7 A
```

中点の `5.0e-7 A` ではなく、雑音が小さいビット0側へしきい値が移動します。
これは、広がりの大きいビット1分布による低電流側への裾を考慮した結果です。

## 11. スケール不変性

すべての電流を正の係数 `k` で拡大します。

```text
mu0'    = k mu0
mu1'    = k mu1
sigma0' = k sigma0
sigma1' = k sigma1
```

このとき次式が成立します。

```text
Delta_mu' / sigma0' = Delta_mu / sigma0
Delta_mu' / sigma1' = Delta_mu / sigma1
```

したがって二次方程式の係数と `t_opt` は変化しません。

```text
t_opt' = t_opt
gamma_opt' = k gamma_opt
```

また、Q関数へ入力される正規化距離も変わらないため、BERも不変です。
これが、電流の絶対スケールから最適化精度を分離できる理由です。

## 12. Q推定値との違い

本コードベースは便宜的な受信Qも計算します。

```text
Q_rx = (mu1 - mu0) / (sigma1 + sigma0)
BER_Q = 0.5 erfc(Q_rx / sqrt(2))
```

これは広く使われる簡便な品質指標ですが、不等分散ガウスに対する厳密な単一
しきい値BER最小化そのものではありません。

- `BER_Q` はQ推定値から得る近似BER
- `BER` は選択したしきい値で2つの誤り確率を直接計算した値

したがって、特に `sigma0` と `sigma1` の差が大きい場合、両者は一致しないことが
あります。コードは両方を結果として返し、区別できるようにしています。

## 13. モデル上の制約

### 13.1 単一しきい値への制約

分散が異なる2つのガウス分布に対する一般的な尤度比判定では、数学的には2つの
交点が生じ、複数の判定領域が最適になる場合があります。しかし、通常のNRZ/OOK
受信器は1つの電流比較しきい値を使用します。

この実装は、その物理的な受信器構成に合わせて次の範囲の単一しきい値を最適化
します。

```text
mu0 <= gamma <= mu1
```

したがって、任意の多領域ベイズ分類器の実装ではありません。

### 13.2 等確率ビット

BER式はビット0と1の出現確率をそれぞれ0.5と仮定しています。事前確率が異なる
場合、BERの重みと密度交点条件に事前確率比が追加されます。

### 13.3 ガウス近似

ショット雑音、TIA雑音、RIN雑音などを合成したサンプル電流をガウス分布として
近似しています。非常に小さい光子数、強い非線形性、クリッピング、パターン依存
ISIなどが支配的な場合には、より詳細なモデルが必要です。

### 13.4 スカラーBERと決定論的ISIの関係

本資料の `mu0` と `mu1` は、メモリレスな2レベルの平均値です。帯域制限によって
同じビットのサンプル値が過去パターンごとに分裂する場合、単一の `mu0`、`mu1`
だけでは決定論的ISIを表現できません。

その場合は、まず波形／ISIヘルパーで次の最悪サンプル開口を確認する必要があります。

```text
vertical_eye_opening = min_one_level - max_zero_level
```

この開口をどのように確率分布へ結び付けるかは追加のモデリング判断であり、現在の
コードは決定論的ISI指標を `optimum_threshold()` へ自動的に合成しません。

### 13.5 物理単位を分離した時間領域RX Eye

`oma_ber.time_domain` のPhase 1は、PD入力P0/P1をW、PD出力をA、TIA出力をVとして
別々に保持します。PDの静的変換と任意の飽和を先に適用し、その後にPD電気応答と
TIAトランスインピーダンスを因果フィルタとして適用します。

```text
P_PD(t) [W]
  -> i_signal(t) = g_PD(P_PD(t)) [A]
  -> i_total(t) = i_signal(t) + Idark [A]
  -> H_PD(f) [A/A]
  -> Z_TIA(f) [V/A]
  -> v_RX(t) [V]
```

位相 `phi` での決定論的開口は、正極性の場合、

```text
Veye(phi) = min(v_one(phi)) - max(v_zero(phi))
```

です。反転TIAでは0/1の上下関係を反転して計算します。Phase 1は `Veye(phi)` を最大
にする位相を求めますが、雑音とジッタを含むBER最適化ではありません。

Phase 2では、パターンごとの平均値と雑音分散を使うGaussian mixtureへ拡張し、
しきい値と位相をBER基準で最適化します。

### 13.6 Gaussian-mixtureしきい値最適化

位相 `phi` で得られるビットb、パターン成分jのTIA出力を、

```text
V | (b, j, phi) ~ N(mu_bj(phi), sigma_bj(phi)^2)
```

とします。`mu_bj` はPD/TIA帯域を通過した決定論的ISI波形、`sigma_bj` は各雑音PSDを
発生ノードからPD/TIA応答へ通したRMS電圧です。正極性TIAでは、0/1を等確率、各
ビット内の成分を等重みとしたBERは、

```text
BER(gamma, phi) = 0.5 * [
  mean_j Q((gamma - mu_0j(phi)) / sigma_0j(phi))
  + mean_j Q((mu_1j(phi) - gamma) / sigma_1j(phi))
]
```

です。反転TIAでは電圧の符号を反転して同じ問題へ変換します。単一ガウス2成分の
場合と異なり、Gaussian mixtureの密度は複数回交差し得るため、閉形式の候補だけに
依存しません。実装は次の順序で大域性と実電圧スケールでの精度を確保します。

1. 全成分平均値を少なくとも最大RMS雑音の8倍まで覆う電圧区間を作る。
2. Gaussian tailを `log_ndtr`、成分和を `logsumexp` で評価し、BERが浮動小数点で
   0へunderflowする領域でもlog-BERの順位を保つ。
3. 1025点の決定論的グリッドでlog-BERの最小近傍を選ぶ。
4. 選んだグリッド点の両隣を境界として、有界スカラー最適化する。
5. 最後に元のGaussian-mixture BER式で `threshold_v` とBERを再評価する。公開BER
   は表現可能範囲を下回ると0になり得るが、しきい値探索自体はlog領域で継続する。

各候補位相で上記しきい値最適化を実行し、BERが最小の位相を選びます。同じBERの
場合は有効Qが大きい位相、さらに0.5 UIに近い位相を選びます。

```text
phi_opt = argmin_phi min_gamma BER(gamma, phi)
```

したがってPhase 1の「決定論的最悪Eye開口を最大にする位相」とPhase 2の
「Gaussian-mixture BERを最小にする位相」は、一般には一致しません。

### 13.7 fractional samplingでの雑音分散

フィルタ後の隣接サンプル雑音は相関しています。fractional samplingの線形補間係数
を `alpha` とすると、

```text
V_phi = (1-alpha) V[n] + alpha V[n+1]

Var[V_phi]
  = (1-alpha)^2 Var[V[n]]
    + alpha^2 Var[V[n+1]]
    + 2 alpha (1-alpha) Cov[V[n], V[n+1]]
```

です。平均波形だけを補間し、RMS雑音を線形補間する方法は分散演算と一致しません。
Phase 2はPD/TIAインパルス応答からlag-1共分散を求め、この式で各パターン成分の
`sigma_bj(phi)` を計算します。

### 13.8 スカラーBERとの使い分け

- `calculate_ber_from_oma()` / `calculate_ber_from_optical_levels()`:
  メモリレス2レベルとユーザー指定 `noise_bandwidth_hz` による高速な感度解析。
- `analyze_statistical_tia_eye()`:
  決定論的ISI、実際の離散PD/TIA応答による等価雑音帯域、パターン依存雑音を結合
  した電圧しきい値・位相解析。

両経路は同じ近似ではありません。時間領域応答を使った後に、同じ帯域効果を
スカラーOMAペナルティとして重ねると二重計上になります。Phase 2はランダム雑音
波形を生成するMonte Carloではなく、有限パターンの解析的Gaussian mixtureです。

### 13.9 Phase 3複素応答とBER最適化の接続

Phase 3の測定応答は、信号平均と雑音分散へ同じFIRとして適用されます。PD応答
`h_PD[n]` とTIA応答 `z_TIA[n]` の合成応答を、

```text
h_RX[n] = h_PD[n] * z_TIA[n]
```

とすると、パターン成分の決定論的平均は、入力電流波形との畳み込みです。

```text
mu_bj(phi) = sample_phi {i_pattern[n] * h_RX[n]}
```

PD以前で発生する白色雑音PSD `S_PD` の一定レベルに対するTIA出力分散は、

```text
sigma_PD^2
  = S_PD * fs/2 * sum_n |h_RX[n]|^2
```

TIA入力換算雑音PSD `S_TIA` はTIA応答だけを使います。

```text
sigma_TIA^2
  = S_TIA * fs/2 * sum_n |z_TIA[n]|^2
```

この `mu_bj(phi)` と `sigma_bj(phi)` が、13.6のGaussian-mixture BERへ直接渡され
ます。したがって、1次LPFと実測応答でDC gainが同じでも、ringing、群遅延、
pre/post-cursor ISI、ENBWが違えば、最適しきい値、最適位相、BERは変化します。

FIRへ変換する前に既知の純遅延だけを除去できますが、最適位相を良く見せる目的で
impulse peakを自動移動してはいけません。また、FIR切捨て後のDC gain誤差とenergy
切捨て比を確認せず、測定応答のBERを確定値として扱ってはいけません。

### 13.10 Phase 4の重み付きmixtureと固定しきい値

Phase 4では、パターン成分jに加えてtiming quadrature node rを持つ分布を考えます。

```text
V | (b, j, r, phi0)
  ~ N(mu_bjr(phi0), sigma_bjr(phi0)^2)

mu_bjr(phi0)
  = mu_bj(phi0 + delta_t,r/T_UI)
```

timing nodeの正規化重みを `w_r`、bit bに属するpattern数を `N_b` とすると、各成分
の重みは `w_r/N_b` です。正極性TIAの固定しきい値BERは、

```text
BER_J(gamma, phi0) = 0.5 * [
  sum_r w_r/N_0 * sum_j Q(
    (gamma - mu_0jr(phi0))/sigma_0jr(phi0)
  )
  + sum_r w_r/N_1 * sum_j Q(
    (mu_1jr(phi0) - gamma)/sigma_1jr(phi0)
  )
]
```

です。`optimize_gaussian_mixture_threshold()` は、従来の等重み成分に加えて、正の
`zero_weights` と `one_weights` を受け取れるよう拡張されています。内部で各重みを
正規化し、

```text
log(sum_k w_k Q(x_k))
```

を `log_ndtr` と `logsumexp` で評価します。大量のtiming×pattern成分で一時配列が
過大にならないよう、しきい値軸をchunk分割します。

各nominal phaseの `gamma` は全timing nodeを結合した後に1回だけ最適化します。

```text
min_gamma sum_r w_r BER_r(gamma)
```

は、

```text
sum_r w_r min_gamma BER_r(gamma)
```

とは異なります。後者は瞬時ジッタごとに受信しきい値を変更する非現実的な下限値
なので採用しません。

`decision_delay_s` はPD/TIAの因果遅延と1 UI内の位相座標を分離します。BER比較では
ジッタなし／ありの双方へ同じ値を渡します。これは表示用の自動Eye中央寄せではなく、
実際のサンプル時刻 `k*T_UI + phi0*T_UI + decision_delay_s + delta_t` の一部です。

## 14. TX/linkから最適しきい値までの実装データフロー

### 14.1 RX入力経路

BERコアへは、OMA直接入力、汎用TX/link、WDMリンクの各経路から到達できます。

従来のOMA入力経路：

```text
calculate_ber_from_oma(oma_dbm, er_db, pd, rx)
  -> OMAをWへ変換
  -> OMA/ERからP0/P1を復元
  -> calculate receiver result
```

汎用TX/link経路：

```text
ExternalLaser
  -> OOKModulator（光1レベル基準IL）
  -> OpticalLossの列
  -> PD入力P0/P1
  -> calculate_ber_from_optical_levels(p0_w, p1_w, pd, rx)
```

WDMリンク経路：

```text
WdmChannel
  -> External Laser CW power
  -> laser-to-MOD coupling loss
  -> OOK MOD（光1レベル基準IL）
  -> AWG split loss
  -> TX fiber coupling loss
  -> OpticalFiberLink（length_km * attenuation_db_per_km）
  -> RX fiber coupling loss
  -> Design margin（任意）
  -> PD-input P0/P1
  -> calculate_ber_from_optical_levels()
```

`analyze_wdm_channel()` は1波長、`analyze_wdm_link()` は共通ファイバを通る複数
波長を解析します。各チャネルの最終P0/P1が同じRX境界コアへ渡されるため、WDM
経路でもガウス雑音と最適しきい値の計算式は変わりません。

`calculate_ber_from_oma()` は互換ラッパーであり、最終的な電流・雑音・しきい値計算
はP0/P1直接入力と共通です。したがって、同じP0/P1なら入力経路によらず同じBERに
なります。

### 14.2 WDMレベル点の意味

WDMレベルダイヤでは、MODより前と後で保持する物理量が異なります。

```text
External Laser / MOD input:
  cw_power_w, cw_power_dbm

MOD output / AWG / Fiber / RX / PD input:
  p0_w, p1_w, pavg_w, oma_w, oma_dbm, er_db
```

各 `WdmLevelPoint` の `incremental_loss_db` は直前要素の損失、
`cumulative_loss_db` はExternal Laserからの累積損失です。MOD出力以降の累積値は
光1レベル経路を基準とします。

```text
cumulative_one_level_loss_db
  = laser_to_modulator_coupling_loss_db
    + modulator_insertion_loss_db
    + awg_split_loss_db
    + tx_fiber_coupling_loss_db
    + fiber_propagation_loss_db
    + rx_fiber_coupling_loss_db
    + margin_db
```

ファイバ伝搬損失は次式です。

```text
fiber_propagation_loss_db
  = length_km * attenuation_db_per_km
    + additional_loss_db
```

各損失はP0とP1へ同じパワー伝送率として適用するため、現在のモデルではERが保存
されます。波長ごとのAWG損失や結合損失は各 `WdmChannel` で個別に設定できます。
ファイバのdB/kmも `fiber_attenuation_db_per_km` でチャネル別に上書きでき、未指定
なら共通 `OpticalFiberLink` の値を使用します。

### 14.3 PD入力から電流分布まで

`calculate_ber_from_optical_levels()` は、PD入力基準面のP0/P1を検証し、次の順序で
計算します。

```text
P0/P1
  -> photocurrent_a()
  -> I0/I1
  -> total_noise_rms_a()
  -> sigma0/sigma1
  -> ber_from_gaussian_levels()
  -> optimum_threshold()
```

返り値の `i0_a`、`i1_a`、`threshold_a` は暗電流を除いた光電流座標です。
DC動作点確認用には暗電流を加えた値も返します。

```text
i0_total_a = i0_a + dark_current_a
i1_total_a = i1_a + dark_current_a
threshold_total_a = threshold_a + dark_current_a
```

共通暗電流は電流差とBER判定距離を変えませんが、ショット雑音とDC動作点へ寄与
します。

### 14.4 TX/link結果とRX単体比較

`analyze_laser_to_receiver()` の `pd_input` は、最終光基準面のP0/P1を保持します。
この値を `compare_photodiodes_at_optical_levels()` へ渡すと、光リンク条件を固定した
ままPDだけを交換できます。

```python
comparison = compare_photodiodes_at_optical_levels(
    p0_w=system_result.pd_input.p0_w,
    p1_w=system_result.pd_input.p1_w,
    photodiodes={"baseline": pd_a, "candidate": pd_b},
    rx=rx,
)
```

各PDについて、同じP0/P1から `I0/I1`、`sigma0/sigma1`、最適しきい値、Q、BERを
再計算します。この分離により、受光感度、暗電流、Fano係数、熱雑音、静的飽和の
寄与をTX/link変動と混ぜずに比較できます。

一方、PDの `bandwidth_3db_hz` と `capacitance_f` は現在このスカラーBER経路へ
自動的には入らないため、帯域・ISI寄与まで含むPD比較では別の波形解析が必要です。

## 15. 検証方針

実装は次の観点でテストされています。

1. µAスケールの既知条件で解析的なしきい値が得られること。
2. nAスケールの条件で、無次元グリッドによる独立参照探索と一致すること。
3. `mu` と `sigma` を同倍率で変えても正規化しきい値が変わらないこと。
4. 等分散では平均電流の中点を返すこと。
5. 極端な分散比でも、選択したしきい値のBERが中点BER以下であること。
6. しきい値が常に平均電流間の候補から選ばれること。

この検証により、単に「しきい値が2レベルの間にある」だけでなく、実電流スケール
でBER最小化の精度が維持されることを確認しています。

さらに、TX/link/RX分離について次を検証しています。

7. 光1基準ILから計算したP1がレーザーパワーとILの式に一致すること。
8. 受動損失後もERが保存され、P0/P1/OMAが同じ伝送率で減衰すること。
9. レベルダイヤの累積損失とOMA低下量が一致すること。
10. 統合解析とP0/P1直接入力で同じBERが得られること。
11. 同一P0/P1で受光感度の高いPDが高いQを示すこと。
12. WDMレベルダイヤが指定されたTX、ファイバ、RX基準面の順序になること。
13. 距離とdB/kmから求めたファイバ損失が累積損失へ正しく加算されること。
14. 複数波長が一意なチャネル名と波長で個別に評価されること。

時間領域Phase 2については次も検証しています。

15. 理想離散応答の等価雑音帯域がNyquist帯域 `fs/2` と一致すること。
16. 高oversamplingの1次LPFが連続時間値 `pi*f3dB/2` へ収束すること。
17. 各雑音成分が一定レベルで `PSD * Bn * |gain|^2` と一致すること。
18. フィルタ後の隣接サンプル共分散がfractional sampling分散へ入ること。
19. メモリレス単一成分のGaussian mixtureがスカラーBERと一致すること。
20. 通常BERが0へunderflowする条件でもlog-BER探索が中点しきい値を保つこと。
21. 統計Eye密度の電圧積分が各位相で1となること。
22. 選択位相のBERが全候補位相の最小値であること。

## 16. まとめ

実装全体と、実電流スケールで安定したしきい値を得るための要点は次のとおりです。

- External Laser、光1基準ILの変調器、受動損失をP0/P1レベルで接続する。
- WDMチャネルごとにAWG、結合、距離依存ファイバ損失を追跡する。
- PD入力P0/P1をTX/linkとRXの明示的な境界にする。
- 既存OMA APIと新しいP0/P1 APIで同じBERコアを共有する。
- 同一P0/P1でPDだけを交換し、PDスペック寄与を分離する。
- A単位の固定絶対許容誤差に依存しない。
- しきい値をレベル間隔で無次元化する。
- ガウス密度交点を解析的に求める。
- 二次方程式を係数正規化と安定した解の公式で解く。
- 区間内の全交点と端点を元のBER式で比較する。
- Q推定BERと、最適しきい値による直接BERを区別する。
- 雑音源ごとの片側PSDを物理的な発生ノードからPD/TIA応答へ伝搬する。
- 離散応答から等価雑音帯域、パターン依存分散、隣接共分散を求める。
- Gaussian-mixture BERをlog領域で比較し、電圧しきい値と位相を最適化する。

この構造により、TX/linkのレベルダイヤとRX単体のPD検証を分離しながら、nA、µA、
mAの違いに左右されない一貫した正規化しきい値とBERを得られます。
