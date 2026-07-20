# OMA-to-BER／フォトダイオード受信系の物理モデル

## 1. この資料の目的と適用範囲

この資料は、本リポジトリが実装するOMA-to-BER解析の物理モデル、計算順序、
単位、仮定、適用限界を説明します。

中心となるモデルは、NRZ/OOK信号を対象にしたスカラー・ガウス雑音モデルです。
加えて、物理単位を分離した時間領域経路では、決定論的ISIとPSDベースのガウス
雑音モーメントをパターン依存Gaussian mixtureとして結合できます。光OMAと消光比
から光レベルを求め、フォトダイオード（PD）の電流へ変換し、レベル依存雑音、
判定しきい値、Q推定値、BERを計算します。

```text
OMA, ER
  -> 光0/1レベル P0, P1
  -> PD電流レベル I0, I1
  -> レベル別RMS雑音 sigma0, sigma1
  -> 単一判定しきい値 gamma
  -> Q推定値、BER
```

帯域ペナルティ、決定論的波形、ISI、アイダイアグラム、周波数応答、リンク
バジェット、PD飽和の補助モデルもあります。External Laserから受信器までの
受動リンクは `analyze_laser_to_receiver()` で明示的に接続できます。帯域、ISI、
周波数応答はスカラーBERへ自動的には反映されませんが、`oma_ber.time_domain` の
Phase 2ではPD/TIA応答から等価雑音帯域と統計Eyeを計算します。

本ツールは受信器のシステムレベル設計議論を目的とした簡略モデルです。TDECQなど
の規格適合性、ランダム雑音サンプルを生成するMonte Carlo波形シミュレーション、
デバイス物理シミュレーションを提供するものではありません。

最適しきい値の導出と数値安定化については、
[実電流スケールにおけるガウスBER最適しきい値](explanation.md)も参照してください。

## 2. 単位と対数表現

### 2.1 dBmとW

dBmは1 mWを基準とする絶対光パワーです。

```text
P_W = 1e-3 * 10 ** (P_dBm / 10)
P_dBm = 10 * log10(P_W / 1e-3)
```

例：

```text
  0 dBm = 1e-3 W = 1 mW
-10 dBm = 1e-4 W = 0.1 mW
```

`watt_to_dbm()` では、対数が定義できない0以下のパワーを拒否します。

### 2.2 dBと線形比

ERなどのパワー比には10を係数とする変換を使います。

```text
ratio_linear = 10 ** (ratio_dB / 10)
ratio_dB = 10 * log10(ratio_linear)
```

一方、電圧・電流の振幅比や電流OMA圧縮量には20を使います。

```text
amplitude_dB = 20 * log10(amplitude_ratio)
```

パワー比と振幅比を混同すると、対数値が係数2倍でずれるため、コードでは用途ごとに
式を分けています。

### 2.3 RIN

RINは通常dB/Hzで与えられ、線形値の単位は1/Hzです。

```text
RIN_linear_per_Hz = 10 ** (RIN_dB_per_Hz / 10)
```

例えば `-150 dB/Hz` は `1e-15 1/Hz` です。

## 3. NRZ/OOKの光レベル

### 3.1 OMA

NRZ/OOKの光変調振幅（Optical Modulation Amplitude）は、光1レベルと光0
レベルの差です。

```text
OMA = P1 - P0
```

OMAをdBmで表す場合、`P1` や平均光パワーのdBmではなく、差分パワーOMAをWへ
換算した値のdBmです。

### 3.2 消光比

消光比（Extinction Ratio; ER）は2レベルのパワー比です。

```text
ER_linear = P1 / P0
ER_linear = 10 ** (ER_dB / 10)
```

有限の正しい2レベルを得るには `ER_linear > 1`、すなわち `ER_dB > 0` が必要
です。

### 3.3 OMAとERからのレベル復元

```text
OMA = P1 - P0
ER  = P1 / P0
```

を連立すると、次式になります。

```text
P0   = OMA / (ER - 1)
P1   = ER * OMA / (ER - 1)
Pavg = (P0 + P1) / 2
```

`Pavg` は0と1が等確率で出現する場合の平均光パワーです。マーク率が0.5でない
信号では、この単純平均は実際の時間平均と一致しません。

実装は `OMA > 0` と `ER > 1` を検証し、ゼロ変調やレベル反転を拒否します。

## 4. フォトダイオード電流

### 4.1 線形受光感度モデル

PDの受光感度 `Rpd` [A/W]を使うと、線形領域の光電流は次式です。

```text
I0 = Rpd * P0
I1 = Rpd * P1
Delta_I = I1 - I0 = Rpd * OMA
```

`Delta_I` は電気的な信号振幅であり、受光感度が高いほど同じ光OMAから大きな
電流アイ開口が得られます。

`Photodiode` は次の物理パラメータを保持します。

- `responsivity_a_per_w`: 受光感度 [A/W]
- `dark_current_a`: 暗電流 [A]
- `bandwidth_3db_hz`: 3 dB帯域 [Hz]
- `capacitance_f`: 容量 [F]
- `saturation_power_w`: 飽和モデルの基準光パワー [W]
- `return_loss_db`: リターンロス [dB]
- `bias_v`: バイアス電圧 [V]
- `shunt_resistance_ohm`: シャント抵抗 [ohm]
- `temperature_k`: 温度 [K]
- `dark_current_fano_factor`: 暗電流ショット雑音電力の係数 [無次元]

現在の主スカラーBER計算へ直接入るのは、受光感度、暗電流、暗電流Fano係数、
任意の飽和パワー、任意のシャント抵抗と温度です。3 dB帯域、容量、リターンロス、
バイアス電圧は保存されますが、自動的にはBERへ反映されません。

主計算の `i0_a` と `i1_a` は暗電流を含まない光電流座標です。DC動作点確認用に、
共通暗電流を加えた総電流も返します。

```text
i0_total_a = i0_a + dark_current_a
i1_total_a = i1_a + dark_current_a
```

共通暗電流は0/1の電流差を変えません。

```text
i1_total_a - i0_total_a = i1_a - i0_a = Delta_I
```

同様に、`threshold_a` は光電流座標、`threshold_total_a` はDC総電流座標です。

```text
threshold_total_a = threshold_a + dark_current_a
```

この共通オフセットは理想的な線形受信器のBERを変えませんが、実回路ではTIAの
DCヘッドルームやオフセット補償範囲を評価するために重要です。現在の実装は総電流
を報告しますが、TIAのDC飽和判定までは行いません。

### 4.2 tanh静的飽和モデル

`saturation_power_w` が指定された場合、主BERフローは次の簡略化したtanhモデル
を使用します。

```text
I(P) = Rpd * Psat * tanh(P / Psat)
```

低入力では `tanh(x) ~= x` のため、線形モデルへ近づきます。

```text
P << Psat なら I(P) ~= Rpd * P
```

高入力では電流が次の値へ漸近します。

```text
P >> Psat なら I(P) -> Rpd * Psat
```

圧縮された電流OMAと線形予測値の差は次式で表します。

```text
compression_dB = 20 * log10(Delta_I_linear / Delta_I_compressed)
```

非常に高い入力では `I0` と `I1` が同じ飽和電流へ近づくため、`Delta_I` は単調に
増えません。途中で最大となった後に縮小し、BERが再び悪化する場合があります。
これはこの静的モデルの重要な性質です。

`Psat` はモデルの形を決めるパラメータであり、校正なしに実デバイスの規定限界や
1 dB圧縮点と同一視してはいけません。

## 5. レベル依存雑音モデル

各レベル `i` の総RMS電流雑音は、独立と仮定した雑音電力の二乗和平方根です。

```text
sigma_i = sqrt(
    sigma_shot_i**2
    + sigma_tia**2
    + sigma_rin_i**2
    + sigma_thermal**2
)
```

すべての `sigma` の単位はA RMSです。独立性を仮定しているため、相関雑音やパターン依存雑音はこの和には含まれません。

### 5.1 ショット雑音

```text
sigma_shot_i = sqrt(2 * q * (I_i + F_dark * I_dark) * Bn)
```

- `q = 1.602176634e-19 C`: 電気素量
- `I_i`: レベル `i` の光電流 [A]
- `I_dark`: 暗電流 [A]
- `F_dark`: 暗電流Fano係数、既定値1.0 [無次元]
- `Bn`: 等価雑音帯域 [Hz]

暗電流は信号ではありませんが、キャリアの離散到来に伴うショット雑音へ寄与
します。`F_dark=1` は理想Poissonショット雑音です。1より大きい値は、測定された
暗電流雑音が理想値を上回る場合の簡略な過剰雑音係数として使えます。信号光電流
項の係数は現在1に固定されています。飽和モデル有効時には、`I_i` として圧縮後の
光電流を使います。

### 5.2 TIA入力換算電流雑音

白色の入力換算電流雑音密度を、指定した雑音帯域で積分する簡略式です。

```text
sigma_tia = i_n * sqrt(Bn)
```

`i_n` の単位はA/sqrt(Hz)です。実際のTIA雑音密度が周波数依存する場合、本来は伝達関数を含むスペクトル積分が必要です。本モデルでは、ユーザーが適切な等価雑音帯域と等価白色雑音密度を与えるものとします。

### 5.3 RIN雑音

```text
sigma_rin_i = Rpd * P_i * sqrt(RIN_linear_per_Hz * Bn)
```

RINは光パワーに比例する相対雑音なので、通常は光1レベル側で大きくなります。
`rin_db_per_hz=None` の場合、RIN項は0です。

tanh飽和が有効でも、現在のRIN項は小信号受光感度 `Rpd` による上式を使用します。
より詳細には、動作点での微分受光感度を使うモデルが考えられます。

```text
dI/dP = Rpd * sech(P / Psat)^2
```

ただし、この微分飽和RINモデルは現在の主計算には実装されていません。

### 5.4 PDシャント抵抗の熱雑音

`shunt_resistance_ohm` が指定された場合、Johnson-Nyquist電流雑音を加えます。

```text
sigma_thermal = sqrt(4 * k_B * T * Bn / Rsh)
```

- `k_B = 1.380649e-23 J/K`: ボルツマン定数
- `T`: 絶対温度 [K]
- `Rsh`: シャント抵抗 [ohm]

シャント抵抗が未指定なら、この項は0です。

### 5.5 雑音帯域と3 dB帯域の違い

`noise_bandwidth_hz` と `bandwidth_3db_hz` は同じ量ではありません。

- 3 dB帯域は伝達関数振幅が低周波値の `1/sqrt(2)` になる周波数です。
- 等価雑音帯域は、実際の雑音伝達関数を同じ雑音電力の矩形帯域へ置き換えた幅です。

例えば1次ローパスの等価雑音帯域は理想条件で `pi/2` 倍の3 dB帯域ですが、実際の
PD/TIAでは伝達関数、ピーキング、後段フィルタで関係が変わります。本ツールは
`noise_bandwidth_hz` を自動導出せず、明示的な入力として扱います。

## 6. ガウスBERモデル

### 6.1 2つの電流分布

サンプリングされた電流を次の正規分布で近似します。

```text
bit 0: N(mu0, sigma0^2)
bit 1: N(mu1, sigma1^2)
```

主NRZ/OOKフローでは次の対応です。

```text
mu0 = I0
mu1 = I1
```

単一判定しきい値 `gamma` に対する判定規則は次のとおりです。

```text
I < gamma  -> bit 0
I >= gamma -> bit 1
```

条件付き誤り確率は次式です。

```text
P(error | 0) = Q((gamma - mu0) / sigma0)
P(error | 1) = Q((mu1 - gamma) / sigma1)
```

0と1が等確率なら、BERは次式です。

```text
BER = 0.5 * [
    Q((gamma - mu0) / sigma0)
    + Q((mu1 - gamma) / sigma1)
]
```

Q関数は次式です。

```text
Q(x) = 0.5 * erfc(x / sqrt(2))
```

### 6.2 最適しきい値

`optimize_threshold=True` の場合、`mu0 <= gamma <= mu1` の単一しきい値の中で
BERを最小にします。

等分散では中点です。

```text
sigma0 = sigma1 なら gamma = (mu0 + mu1) / 2
```

不等分散では、ガウス密度の交点条件を使います。

```text
phi((gamma - mu0) / sigma0) / sigma0
=
phi((mu1 - gamma) / sigma1) / sigma1
```

実装はA単位の固定絶対許容誤差を使う数値最適化ではなく、電流レベル間隔で無次元化した二次方程式を解析的に解きます。区間内の全停留点と両端でBERを評価し、最小の候補を返します。このため、nA、µA、mAの電流スケールを同倍率で変更しても、正規化しきい値とBERは不変です。

詳細な導出、係数、安定した二次方程式解法、数値例は[最適しきい値の詳細説明](explanation.md)を参照してください。

`optimize_threshold=False` の場合は、分散差にかかわらず中点を使用します。

### 6.3 単一しきい値モデルの意味

不等分散ガウスの一般的な尤度比判定では、数学的に複数の判定領域が最適になる
場合があります。本実装は通常のNRZ/OOK比較器を意図し、平均電流間にある単一
しきい値だけを最適化します。任意の多領域ベイズ分類器ではありません。

## 7. Q推定値

コードは一般的な受信Q推定値も返します。

```text
Q_rx = (mu1 - mu0) / (sigma1 + sigma0)
BER_Q = 0.5 * erfc(Q_rx / sqrt(2))
```

`BER_Q` は簡便な近似です。`calculate_ber_from_oma()` が返す `ber` は、選択した
しきい値をガウスBER式へ直接代入した値です。不等分散が大きい場合、`ber` と
`ber_from_q` は一致しないことがあります。

## 8. 目標BERに必要なOMA

`required_oma_dbm()` は次式の根を探索します。

```text
BER(OMA_required) = target_BER
```

探索範囲の下端ですでに目標を満たす場合は下端を返します。下端と上端で目標BERを
挟める場合は、スカラー根探索でOMAを求めます。

### 8.1 単調性の前提

線形PDモデルの代表的な範囲では、OMAを増やすと信号電流差が増え、BERは単調に
改善します。この条件では、両端で根を挟む方法が有効です。

一方、tanh飽和モデルでは高入力で `I0` と `I1` が同じ飽和電流へ近づくため、
BERは非単調になり得ます。途中で目標BERを満たしていても探索上端で再び悪化して
いると、現在の `required_oma_dbm()` は「到達不能」と判定する可能性があります。

したがって、飽和モデルを有効にしたrequired-OMA計算では次が必要です。

- 事前にOMA sweepでBERの単調区間を確認する。
- 探索範囲を最初のBER改善区間へ限定する。
- 強い飽和を横断する範囲で、単一のrequired OMA結果を無条件に信用しない。

これは現在の探索アルゴリズムの既知の制約です。

## 9. 1次ローパス帯域ペナルティ

簡易帯域モデルは1次ローパスの振幅応答を使います。

```text
|H(f)| = 1 / sqrt(1 + (f / f3dB)^2)
```

対応する振幅ペナルティは次式です。

```text
penalty_dB = -20 * log10(|H(f)|)
effective_OMA_dBm = OMA_dBm - penalty_dB
```

`f = f3dB` では `|H| = 1/sqrt(2)`、ペナルティは約3.0103 dBです。

これは単一周波数の振幅減衰をOMAペナルティへ換算する工学的補助式です。NRZ信号
の全スペクトル、立上り、ISI、等価雑音帯域を同時に扱うものではありません。

## 10. 決定論的波形、ISI、アイダイアグラム

本章の10.1から10.5は、既存の簡易波形／ISI互換経路を説明します。物理ドメインを
分離した新しい受信機Eye経路は10.6以降で説明します。

### 10.1 NRZ波形

波形ヘルパーは0/1ビット列を、指定した低レベル・高レベルの矩形NRZサンプル列へ
変換します。固定seedの疑似ランダム列は再現可能ですが、ここでのPRBSヘルパーは
通信規格で定義された特定多項式のPRBS生成器を意味しません。

### 10.2 ISI

ISI（Inter-Symbol Interference）は、あるシンボルの応答が後続シンボル時刻まで
残り、同じ論理値でも過去パターンによってサンプル値が変わる現象です。

理想的にはシンボル中央の値が2レベルへ一意に対応します。

```text
0 bit -> P0 または I0
1 bit -> P1 または I1
```

帯域が不足すると遷移が遅くなります。例えば長い0の後の1はサンプル時刻までに
高レベルへ到達せず、長い1の後の0は低レベルへ整定しない場合があります。

```text
...0001 の1サンプル != ...1111 の1サンプル
...1110 の0サンプル != ...0000 の0サンプル
```

このパターン依存性は、ランダム雑音がなくてもアイを閉じます。

### 10.3 インパルス応答と畳み込み

ISIヘルパーは入力波形とインパルス応答を線形畳み込みします。1次ローパスの
離散インパルス応答はDCゲインが1になるよう正規化されるため、定常レベルを保ち
ながら遷移だけを帯域制限できます。

畳み込み出力は入力波形と同じ長さへ切り出されます。初期過渡を含むため、評価時は
必要に応じて先頭シンボルを除外する必要があります。

### 10.4 サンプルド・アイ指標

主な決定論的指標は次のとおりです。

- 0サンプル平均値
- 1サンプル平均値
- 1サンプル最小値
- 0サンプル最大値
- 垂直アイ開口
- 簡易ISIペナルティ

最悪値による垂直アイ開口は次式です。

```text
vertical_eye_opening = min_one_level - max_zero_level
```

```text
vertical_eye_opening > 0 なら決定論的サンプルアイは開いている
vertical_eye_opening <= 0 なら単一しきい値で全サンプルを分離できない
```

簡易ISIペナルティは次式です。

```text
isi_penalty_db = 20 * log10(
    mean_level_separation / vertical_eye_opening
)
```

アイ開口が0以下なら、この指標は無限大です。これは工学的診断値であり、規格適合
指標ではありません。

### 10.5 アイダイアグラム

アイヘルパーは波形を通常2 UIの重なりトレースへ分割します。matplotlibによる
表示は帯域制限ISIの可視化に有用ですが、ランダム雑音、ジッタ、CDR、BER等高線を
自動追加しません。

光OMAは送信側の `P1 - P0` ですが、受信判定余裕はフィルタ後のサンプル開口です。

```text
min_one_level - max_zero_level < P1 - P0
```

となる場合、名目OMAから予想するより有効QとBERが悪化し得ます。現在の決定論的
ISI指標は、スカラー・ガウスBERへ自動結合されません。

### 10.6 物理ドメインを分離したRX Eye

`oma_ber.time_domain` は、簡易Eye経路とは別に、PD入力から判定点までを次の単位で
分離します。

```text
PD入力光パワー P(t) [W]
  -> PDの静的光電変換 i_raw(t) [A]
  -> PD電気応答 H_PD(f) [A/A]
  -> TIAトランスインピーダンス Z_TIA(f) [V/A]
  -> TIA出力 v(t) [V]
```

線形PDでは、信号光電流と暗電流を次式で分けます。

```text
i_signal(t) = Rpd * P(t)
i_total(t)  = i_signal(t) + Idark
```

`saturation_power_w` が設定されている場合は、既存の静的tanh圧縮をPD帯域より前に
適用します。

```text
i_signal(t) = Rpd * Psat * tanh(P(t) / Psat)
```

この順序により、飽和後の電流波形をPD/TIAの因果的電気応答へ入力できます。PD帯域
を光パワーWのまま処理しないため、TIA出力まで進んでも単位を追跡できます。

### 10.7 厳密な時間グリッドとPRBS

`TimeGrid` はシンボルレートと整数samples/UIから時間軸を作ります。

```text
T_UI = 1 / symbol_rate_baud
sample_rate_hz = symbol_rate_baud * samples_per_symbol
dt = 1 / sample_rate_hz
```

従来ヘルパーのように `sample_rate_hz / symbol_rate_baud` を丸めません。新経路の
`prbs_bits()` はPRBS7、PRBS9、PRBS15、PRBS31の多項式を使います。従来のseed付き
乱数ヘルパーは互換経路の再現可能な例示用として残します。

### 10.8 因果フィルタと初期過渡

PDとTIAは、分子・分母係数、サンプルレート、応答種別を保持する
`DiscreteTransferFunction` で表します。

```text
response_kind = "dimensionless"       # PD電流応答 A/A
response_kind = "transimpedance_ohm"  # TIA応答 V/A
```

1次LPFはprewarp付き双一次変換を使い、離散時間上の指定周波数で3 dBとなるように
します。複数段フィルタでは、各段の極から求めた整定サンプル数を加算し、評価ビット
列の周期的な直前パターンをwarm-upとして与えます。warm-up部分はEye評価前に除外
されるため、ゼロ初期状態を最悪レベルへ混入させません。

### 10.9 クロック位相と判定中心Eye

ビットkのサンプル時刻は次式です。

```text
t_sample(k) = k * T_UI + sampling_phase_ui * T_UI + decision_delay_s
```

非整数サンプル位置は線形補間します。`decision_delay_s` はチャネル全体遅延とビット
ラベルの対応、`sampling_phase_ui` は1 UI内の判定位相を表すため、両者を混同しません。

各位相で決定論的開口を計算します。正極性では、

```text
eye_opening(phi) = min(one(phi)) - max(zero(phi))
```

反転TIAでは、

```text
eye_opening(phi) = min(zero(phi)) - max(one(phi))
```

です。実装は平均レベルから極性を判定し、最大開口位相を選びます。同値位相が複数
ある場合は0.5 UIに最も近い位相を選びます。表示トレースは判定時刻を0 UIとして
`-1 UI` から `+1 UI` まで切り出し、2 UI終端を含みます。

これはISIだけを考えた最適位相です。ランダムジッタ、CDR追従、雑音を含むBER最適
位相ではありません。

### 10.10 Phase 1とPhase 2の境界

新しいRX EyeのPhase 1に含むものは次のとおりです。

- PD入力P0/P1からの矩形NRZ/OOK光パワー
- PD responsivity、暗電流DC、静的飽和
- 因果的なPD/TIA 1次応答または明示的な離散伝達関数
- PRBS7/9/15/31
- warm-up、fractional sampling、遅延、正負極性
- PD電流EyeとTIA電圧Eye

Phase 2で追加した項目は次のとおりです。

- 光電流ショット、暗電流ショット、TIA、PDシャント熱、RINの片側PSD
- PD/TIA応答からの等価雑音帯域
- パターン依存分散と隣接サンプル共分散
- Gaussian-mixture統計Eye密度
- BERを最小にする電圧しきい値とクロック位相

次の項目はまだ含みません。

- ランダム雑音サンプルを加えたMonte Carlo波形
- BER等高線、bathtub曲線
- ランダム／周期ジッタ、CDR
- 光電界、変調器chirp、AWG複素応答、ファイバ分散
- WDM beat/crosstalk

したがってPhase 2は、物理単位と処理順を保った解析的なガウス統計RX解析ですが、
ランダム波形シミュレーションまたは規格適合試験ではありません。

### 10.11 雑音源の発生ノード

`calculate_tia_output_noise()` は、すべての雑音を同じ帯域へ入れず、発生位置に応じた
伝達関数を使います。

```text
光電流shot、Idark shot、PDシャント熱、RIN
  -> H_PD(f) * Z_TIA(f)

TIA入力換算電流雑音
  -> Z_TIA(f)
```

各片側入力電流PSDは次式です。

```text
S_photo(t) = 2 q i_signal(t)                         [A^2/Hz]
S_dark     = 2 q F_dark I_dark                       [A^2/Hz]
S_TIA      = i_n^2                                   [A^2/Hz]
S_RIN(t)   = (Rpd P(t))^2 RIN_linear                 [A^2/Hz]
S_thermal  = 4 k_B T / R_shunt                       [A^2/Hz]
```

`I_dark` 自体は0/1の共通DC電流なので理想判定距離を変えませんが、`S_dark` を増やし、
両レベルのRMS雑音を増加させます。`F_dark=1` はPoisson暗電流、1より大きい値は過剰
雑音を表します。RINは光パワーに比例する振幅揺らぎなのでPSDは光パワーの2乗に
比例します。

静的tanh飽和PDでは、小信号RIN変換に局所傾きを使います。

```text
i(P) = Rpd Psat tanh(P / Psat)
di/dP = Rpd {1 - tanh^2(P / Psat)}
S_RIN(t) = {(di/dP) P(t)}^2 RIN_linear
```

これは飽和したDC電流をそのままRIN振幅とみなすより、入力光揺らぎに対する局所
小信号応答を表します。ショット雑音は実際の飽和後信号電流 `i_signal(t)` を使います。

### 10.12 離散応答と等価雑音帯域

サンプルレート `fs` の片側白色PSD `S1` を離散白色雑音へ対応させると、Nyquist
帯域 `0 ... fs/2` に含まれる入力分散は次式です。

```text
variance_input = S1 * fs / 2
```

因果インパルス応答 `h[n]` とDC利得 `H(0)` に対する片側等価雑音帯域は、

```text
Bn = (fs / 2) * sum_n |h[n]|^2 / |H(0)|^2
```

です。したがって、一定PSDの出力分散は `S1 * Bn * |H(0)|^2` になります。
Phase 2は実際に指定された離散PD/TIA応答から `pd_tia_noise_bandwidth_hz` と
`tia_noise_bandwidth_hz` を別々に求めます。スカラーBERのユーザー指定
`noise_bandwidth_hz` を再利用しません。

理想的な無帯域制限離散応答では `Bn=fs/2` です。連続時間1次LPFの
`Bn=pi*f3dB/2` は、サンプルレートが帯域より十分高い場合の極限です。有限の
サンプルレートでは、prewarp付き双一次変換後の離散応答から求めた値を採用します。

### 10.13 パターン依存分散と隣接共分散

ショット雑音とRINのPSDは `P(t)` または `i_signal(t)` に依存します。このため、
出力分散もビット履歴に依存します。入力白色雑音の時刻mにおける分散を `u[m]`、
伝達関数のインパルス応答を `h[k]` とすると、

```text
Var[y[n]] = sum_k h[k]^2 u[n-k]
Cov[y[n], y[n+1]] = sum_k h[k] h[k+1] u[n-k]
```

です。実装は評価PRBSを周期的に延長してPSDを畳み込むため、雑音分散の開始過渡を
Eyeへ混入させません。

非整数サンプル位置を線形補間すると、平均値だけでなく分散も変わります。

```text
y_phi = (1-alpha) y[n] + alpha y[n+1]

Var[y_phi]
  = (1-alpha)^2 Var[y[n]]
    + alpha^2 Var[y[n+1]]
    + 2 alpha (1-alpha) Cov[y[n], y[n+1]]
```

共分散項を省くと、ローパス後の強く相関した雑音をfractional samplingした際にRMS
雑音を過小評価するため、Phase 2ではlag-1共分散も保持します。

### 10.14 Gaussian mixtureとBER

位相 `phi` でサンプルしたビットbの各パターン成分jを、

```text
V | (b, j, phi) ~ N(mu_bj(phi), sigma_bj(phi)^2)
```

と近似します。`mu_bj` は決定論的TIA波形、`sigma_bj` はPSD伝搬結果から得ます。
0/1は各0.5の事前確率、各ビット内の有限パターン成分は等重みです。正極性TIAの
単一電圧しきい値 `gamma` に対するBERは、

```text
BER(gamma, phi) = 0.5 * [
  mean_j Q((gamma - mu_0j(phi)) / sigma_0j(phi))
  + mean_j Q((mu_1j(phi) - gamma) / sigma_1j(phi))
]
```

です。反転TIAは電圧座標の符号を内部で反転して同じ式を使います。実装は広い電圧
範囲の1025点を評価して大域的な候補を探し、その近傍だけを有界1次元最適化します。
各位相でしきい値を最適化した後、最小BERの位相を採用します。BERが同値の場合は
有効Qが大きい位相、さらに0.5 UIへ近い位相を選びます。

有効Q表示値は、各ビット内のガウス雑音とパターン平均値の広がりをまとめた診断値
です。

```text
sigma_eff,b^2
  = mean_j {sigma_bj^2 + (mu_bj - mean_j(mu_bj))^2}

Q_eff = |mean(mu_1j) - mean(mu_0j)| / (sigma_eff,0 + sigma_eff,1)
```

BERはこのQ近似からではなく、上記Gaussian-mixture BERを直接最小化して求めます。

### 10.15 統計Eye密度と解釈上の制約

`analyze_statistical_tia_eye()` は、各位相で0/1成分の正規密度を足し合わせた
`density_per_v` [1/V]を返します。`plot_statistical_eye_density()` は1 UIの密度を
2 UIへ周期表示し、BER最適位相が0 UIとなるよう中心合わせします。

この密度は有限PRBSに対する解析的Gaussian mixtureです。次の効果は含みません。

- 非ガウスな低光子数ショット統計
- 雑音源間の相関
- 周波数依存のTIA雑音密度やRIN
- ジッタと電圧雑音の時間・振幅結合
- CDR追従、判定器帯域、ヒステリシス、offset
- 光電界伝搬、分散、chirp、WDM beat/crosstalk

したがって結果は、白色PSD、線形電気応答、ガウス近似の範囲内でのシステム設計
値です。規格適合性や極低BERの保証にはMonte Carlo、importance sampling、回路／
光伝搬シミュレーション、実測との比較が別途必要です。

### 10.16 Phase 3の実測複素応答統合

Phase 3では、1次LPFの代わりにPD電流応答 `H_PD(f)` [A/A]またはTIA
トランスインピーダンス `Z_TIA(f)` [V/A]を複素振幅応答として入力できます。

```text
frequency_hz,magnitude_db,phase_deg
```

振幅は光パワーdBではなく、電流／電圧の振幅dBです。

```text
H(f) = 10^(magnitude_db/20) * exp(j*phase_deg*pi/180)
```

新しい `measured_frequency_response_from_csv()` は、応答の単位種別と物理基準面も
保持します。

```text
response_kind="dimensionless"        -> H_PD: A/A
response_kind="transimpedance_ohm"   -> Z_TIA: V/A
```

実波形へ `irfft` 変換するには、0 HzとNyquistの応答が実数で、負周波数側が
Hermitian共役となる必要があります。Phase 3経路はDC点を必須とし、入力上限が
`sample_rate_hz/2` まで届かない応答を拒否します。DC以下や測定上限以上を暗黙に
holdしません。

既知の純遅延 `tau_removed` だけを次式で除去できます。

```text
H_deembedded(f)
  = H_measured(f) * exp(+j*2*pi*f*tau_removed)
```

これは測定fixtureなどの既知遅延を取り除く操作です。最大インパルスtapを自動で
0へ移動すると、実際のpre-cursorや反射を純遅延と誤認するため、自動移動は行い
ません。

複素応答は、振幅とunwrap位相をrFFTグリッドへ補間してから `irfft` します。

```text
f_k = k*fs/Nfft,  k=0 ... Nfft/2
h_periodic[n] = irfft(H(f_k))
```

有限FFTでは周期末尾が負時間側に相当します。全energyに対する負時間energyと、
採用FIR長より後ろの正時間tail energyを別々に評価します。

```text
E_total = sum_n |h_periodic[n]|^2

pre_echo_energy_ratio
  = E_negative_time / E_total

discarded_tail_energy_ratio
  = E_discarded_positive_time / E_total
```

どちらかが指定許容値を超える場合、`fir_transfer_from_measured_response()` は応答を
自動修復せず `ValueError` を返します。測定のphase reference、除去遅延、FFT時間窓、
FIR長のいずれが不足しているかを入力側で確認します。

検証済みFIRは `DiscreteTransferFunction` となり、Phase 1/2と同じ経路へ入ります。

```text
deterministic signal mean  -> h[n]
white-noise variance       -> h[n]^2
lag-1 noise covariance     -> h[n]h[n+1]
ENBW                        -> fs/2 * sum(h[n]^2)/|H(0)|^2
```

したがって、実測応答を信号Eyeだけへ使い、雑音帯域には別の1次LPFを使う不整合を
避けられます。`MeasuredResponseDiagnostics` は、測定DC gain、FIR DC gain、その相対
誤差、ENBW、低周波群遅延、pre-echo比、tail切捨て比を返します。

Phase 3はCSVから直接A/AまたはV/Aが得られる場合を対象とします。S21から
transimpedanceへのreference impedance変換、fixture de-embedding、passivity
enforcement、非因果応答の修復は実装していません。

### 10.17 Phase 4の残留タイミングジッタ

Phase 4は、CDR後または外部で規定された、データEyeと判定クロックの間に残る時間
誤差を扱います。シンボルkのサンプル時刻は、

```text
t_sample(k)
  = k*T_UI + phi0*T_UI + decision_delay_s + delta_t
```

です。`phi0` はnominal sampling phase、`delta_t` は残留ジッタ [s]です。

初期モデルはGaussian random jitterと、位相が一様な単一sinusoidal jitterです。

```text
delta_t,RJ ~ N(0, sigma_RJ^2)
delta_t,SJ = A_SJ*sin(theta),  theta ~ Uniform(0, 2*pi)
delta_t = delta_t,RJ + delta_t,SJ
```

独立なRJとSJの総RMS表示値は、

```text
sigma_total
  = sqrt(sigma_RJ^2 + A_SJ^2/2)
```

です。秒からUIへの変換は `jitter_ui = jitter_s/T_UI` です。

nominal phase `phi0` と固定電圧しきい値 `gamma` に対するBERは、Phase 2のパターン
依存Gaussian-mixture BERを時間誤差分布について平均します。

```text
BER_J(gamma, phi0)
  = integral BER(gamma, phi0 + delta_t/T_UI)
             p(delta_t) d(delta_t)
```

実装ではGaussian RJをGauss-Hermite求積、SJの一様位相を周期中点求積で決定論的に
積分します。RJとSJを同時指定した場合は積求積を使います。

重要なのは、瞬時ジッタ値ごとにしきい値を再最適化しないことです。実受信器の電圧
しきい値は各クロック揺らぎを知って追従できないため、各nominal phaseに対して
全timing nodeに共通する1つのしきい値を選びます。

```text
gamma_opt(phi0) = argmin_gamma BER_J(gamma, phi0)
phi_opt = argmin_phi0 BER_J(gamma_opt(phi0), phi0)
```

timing nodeとPRBS pattern componentは重み付きGaussian mixtureへまとめます。0/1の
事前確率は各0.5のまま、timing nodeの重みは求積重み、同じbit内のpatternは等重み
です。しきい値探索はPhase 2と同じlog-BER計算を使うため、極低BERでもtiming node
間の順位を通常確率のunderflowで失いません。

ジッタによって `phi < 0` または `phi >= 1` になっても、単純にbit labelを隣接bitへ
付け替えません。

```text
decision target = bit[k]
sample time may enter bit[k-1] or bit[k+1]
```

信号波形、雑音分散、lag-1共分散は有限PRBSの周期延長として補間しますが、判定対象
は元の `bit[k]` に固定します。これにより、UI境界を越えたサンプルを誤って正解bit
として数えることを防ぎます。

`plot_jittered_ber_bathtub()` は、横軸nominal phase [UI]、縦軸jitter-averaged BERの
bathtubを表示します。この曲線は各nominal phaseで固定電圧しきい値を最適化した
結果です。

因果的なPD/TIA応答には群遅延があるため、応答生成時刻を基準にしたままではEye中心
が0.5 UIに現れるとは限りません。`decision_delay_s` は、この回路遅延と1 UI内の
nominal phaseを明示的に分離する座標です。正値は物理サンプル時刻を後ろへ移動し、
波形自体を自動で中央寄せしません。ジッタ有無を比較するときは同じ値を使い、既知
回路遅延またはジッタなし解析で求めた位相整列値として別途報告します。

Phase 4への入力はsamplerから見た残留ジッタです。CDR入力ジッタ、CDR tracking
response、dual-Dirac DJ、DCD、data-dependent jitter、TJ@BER外挿は含みません。
CDR前後を接続するには、次PhaseでジッタPSDとCDR誤差伝達関数を導入します。

## 11. 周波数応答／Sパラメータ補助機能

11.1と11.2は従来の簡易互換経路です。物理単位、因果性、pre-echo、FIR切捨てを
検証してPhase 2へ接続する場合は10.16の `oma_ber.time_domain` Phase 3経路を使います。

### 11.1 CSV形式

`frequency_response_from_csv()` は次の列を読みます。

```text
frequency_hz,magnitude_db,phase_deg
```

`magnitude_db` は光パワー比ではなく、電圧または電流の振幅利得dBです。

```text
H(f) = 10 ** (magnitude_db / 20) * exp(j * phase_rad)
```

周波数は正、厳密な昇順である必要があります。

### 11.2 実インパルス応答への変換

`impulse_response_from_frequency_response()` は複素応答の実部・虚部をreal FFTの
周波数グリッドへ線形補間し、`irfft` で実インパルス応答へ変換します。

現在の境界処理は次のとおりです。

- DCから最小入力周波数までは、最初の複素応答値を保持する。
- 最大入力周波数からNyquistまでは、最後の複素応答値を保持する。
- `sample_rate_hz` は最大入力周波数の2倍より大きくする。

この外挿は便宜的な仮定です。特にDC利得、Nyquist近傍、未測定帯域の応答が重要な
場合は、物理的に妥当な周波数点をCSVへ含める必要があります。

また、任意の測定応答を有限長 `irfft` へ直接変換するため、因果性補正、時間原点
調整、窓関数、パッシビティ保証、デエンベディングは自動では行いません。

## 12. 飽和モデルの比較用ヘルパー

飽和モジュールには、主フローのtanh以外にも静的モデルがあります。

### tanhモデル

```text
I = R * Psat * tanh(P / Psat)
```

### ソフトクリップモデル

```text
I = R * P / (1 + (P / Psat)^order)^(1 / order)
```

### 指数モデル

```text
I = R * Psat * (1 - exp(-P / Psat))
```

### 有理関数モデル

```text
R_eff = R / (1 + (P / Psat)^order)
I = P * R_eff
```

これらは圧縮曲線の定性的比較を目的とし、空間電荷、電界遮蔽、熱、キャリア輸送、
帯域変化などを解くデバイス物理モデルではありません。測定データへ校正しない限り、
パラメータを実デバイスの保証値として扱うことはできません。

## 13. リンクバジェット

### 13.1 WDMチャネルの基準面

波長多重リンクは、波長チャネルごとに次の基準面を追跡します。

```text
External Laser
  -> Laser-to-MOD coupling loss
  -> MOD（光1レベル基準IL）
  -> AWG split/insertion loss
  -> TX fiber coupling loss
  -> OpticalFiberLink propagation loss
  -> RX fiber coupling loss
  -> Design margin（任意）
  -> PD input
```

`WdmChannel` は波長 [nm]、External Laser、MOD、AWG損失、TX/RX結合損失を保持
します。`OpticalFiberLink` は複数チャネルが通過する共通ファイバ条件を保持します。

```text
fiber_propagation_loss_db
  = length_km * attenuation_db_per_km
    + additional_loss_db
```

既定では全波長が `OpticalFiberLink.attenuation_db_per_km` を使います。波長依存の
ファイバ減衰を与える場合は、各 `WdmChannel.fiber_attenuation_db_per_km` で上書き
します。固定の `additional_loss_db` は全チャネルに共通です。

External LaserからMOD入力までのCWパワーは次式です。

```text
P_mod_input
  = P_laser * 10 ** (-laser_to_modulator_coupling_loss_db / 10)
```

MOD出力以降はP0/P1を個別に伝搬します。

```text
P0_out = P0_in * 10 ** (-loss_db / 10)
P1_out = P1_in * 10 ** (-loss_db / 10)
```

`analyze_wdm_channel()` は1波長を解析し、`analyze_wdm_link()` は複数の
`WdmChannel` を共通ファイバ、共通PD、共通Receiver条件で一括解析します。

現在は各波長を独立なスカラー信号として扱います。AWGクロストーク、チャネル間
干渉、四光波混合、ラマン効果、ファイバ分散、合計WDM光パワーによる非線形性は
含みません。また、共通PDモデルを使う一括APIでは受光感度の波長依存性を扱わない
ため、波長ごとにPD特性が異なる場合は `analyze_wdm_channel()` を個別のPDで実行
する必要があります。

### 13.2 External Laserと光1レベル基準の変調器

External Laserは、変調器入力基準面の連続光パワーで指定します。

```text
P_laser_W = 1e-3 * 10 ** (laser_power_dbm / 10)
```

OOK変調器の挿入損失 `IL` は光1レベル基準です。

```text
P1_tx = P_laser_W * 10 ** (-IL_db / 10)
ER_linear = 10 ** (ER_db / 10)
P0_tx = P1_tx / ER_linear
OMA_tx = P1_tx - P0_tx
Pavg_tx = (P0_tx + P1_tx) / 2
```

この定義ではレーザーパワー、IL、ERからOMAが決まります。`specified_oma_dbm` は
独立した測定値または仕様値との整合性確認にだけ使用し、計算レベルを上書きしません。

```text
OMA_error_db = calculated_OMA_dbm - specified_OMA_dbm
```

### 13.3 受動リンクのレベルダイヤ

各 `OpticalLoss` は名前と非負のパワー損失 [dB]を持ちます。同じ伝送率をP0とP1へ
適用します。

```text
transmission = 10 ** (-loss_db / 10)
P0_out = P0_in * transmission
P1_out = P1_in * transmission
```

したがって、レベル非依存な受動損失ではERが保存され、OMAと平均パワーが同じdB
だけ低下します。

```text
ER_out = ER_in
OMA_out_dBm = OMA_in_dBm - loss_db
```

レベルダイヤは各基準面について、要素名、単独損失、累積損失、P0、P1、Pavg、
OMA、ERを保持します。設計マージンは必要に応じて最後の仮想損失点として追加
されます。

### 13.4 従来のOMA専用ヘルパー

従来のリンクバジェットヘルパーは、送信OMAから損失とマージンを引いてPD入力
OMAを計算します。

```text
PD_input_OMA_dBm = TX_OMA_dBm - sum(losses_dB) - margin_dB
```

挿入損失とマージンは非負のdB量として入力します。この計算はdB領域の帳尻計算で
あり、反射、多重干渉、波長依存損失、モード依存損失を自動的には扱いません。

この関数は互換性のため維持されています。P0/P1/Pavgを含むレベルダイヤが必要な
場合は `build_level_diagram()` または `analyze_laser_to_receiver()` を使用します。

### 13.5 RX境界とPD比較

`calculate_ber_from_optical_levels()` は、PD入力基準面の `p0_w` と `p1_w` を直接
受け取ります。これによりTXやリンクを通さず、同一光入力条件でPD仕様だけを変更
できます。

```text
fixed PD-input P0/P1
  -> PD responsivity, dark current, Fano factor, saturation, shunt resistance
  -> I0/I1, sigma0/sigma1, Q, BER
```

`compare_photodiodes_at_optical_levels()` は複数の名前付きPDを同じP0/P1で評価し、
`required_oma_by_photodiode()` は目標BERに必要なPD入力OMAを比較します。

これらの比較で、受光感度、暗電流、過剰暗電流雑音、熱雑音、飽和がシステム感度へ
与える寄与を分離できます。ただし、`bandwidth_3db_hz` と `capacitance_f` はまだ
スカラーBERへ自動結合されないため、帯域・ISI寄与は別途評価する必要があります。

### 13.6 実装クラスと関数の対応

システム構成とPython APIの対応は次のとおりです。

| 物理要素・境界 | 実装 | 主な値 |
|---|---|---|
| External Laser | `ExternalLaser` | `output_power_dbm`, `output_power_w` |
| OOK変調器 | `OOKModulator` | 光1基準`insertion_loss_db`, `er_db` |
| 変調器出力 | `ModulatorOutput` | `levels`, calculated/spec OMA差 |
| 受動損失 | `OpticalLoss` | `name`, `loss_db` |
| 光基準面 | `OpticalLevelPoint` | P0/P1/Pavg/OMA/ER、累積損失 |
| TX/link/RX統合 | `analyze_laser_to_receiver()` | レベルダイヤとRX結果 |
| WDMチャネル | `WdmChannel` | 波長、Laser、MOD、AWG、TX/RX結合損失 |
| 共通ファイバ | `OpticalFiberLink` | 距離、dB/km、追加損失 |
| WDM基準面 | `WdmLevelPoint` | CWまたはP0/P1、区間、累積損失 |
| WDM単一波長 | `analyze_wdm_channel()` | 1波長のレベルダイヤとRX結果 |
| WDM複数波長 | `analyze_wdm_link()` | 名前付きチャネルごとの解析結果 |
| RX直接入力 | `calculate_ber_from_optical_levels()` | PD入力P0/P1からBER |
| PD比較 | `compare_photodiodes_at_optical_levels()` | 同一P0/P1でPDごとの結果 |
| 必要OMA比較 | `required_oma_by_photodiode()` | PDごとの必要入力OMA [dBm] |

統合解析の最小コードは次のとおりです。

```python
from oma_ber import (
    ExternalLaser,
    OOKModulator,
    OpticalLoss,
    Photodiode,
    Receiver,
    analyze_laser_to_receiver,
)

result = analyze_laser_to_receiver(
    laser=ExternalLaser(output_power_dbm=-10.0),
    modulator=OOKModulator(
        insertion_loss_db=3.0,
        er_db=6.0,
        specified_oma_dbm=-14.5,
    ),
    losses=[
        OpticalLoss(name="Fiber coupling", loss_db=1.0),
        OpticalLoss(name="Waveguide", loss_db=1.5),
    ],
    margin_db=0.5,
    pd=Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9),
    rx=Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    ),
)
```

`result.modulator_output` はTX出力、`result.level_diagram` は各光基準面、
`result.pd_input` は最終PD入力、`result.receiver_result` は電流・雑音・Q・BERを
保持します。

```python
for point in result.level_diagram:
    print(point.name, point.cumulative_loss_db, point.p0_w, point.p1_w, point.oma_dbm)

print(result.receiver_result["q_rx"])
print(result.receiver_result["ber"])
```

### 13.7 WDM実装例

```python
from oma_ber import (
    ExternalLaser,
    OOKModulator,
    OpticalFiberLink,
    Photodiode,
    Receiver,
    WdmChannel,
    analyze_wdm_link,
)

channels = [
    WdmChannel(
        name="lambda_1310",
        wavelength_nm=1310.0,
        laser=ExternalLaser(output_power_dbm=-4.0),
        modulator=OOKModulator(insertion_loss_db=3.0, er_db=6.0),
        laser_to_modulator_coupling_loss_db=1.0,
        awg_split_loss_db=2.0,
        tx_fiber_coupling_loss_db=1.0,
        rx_fiber_coupling_loss_db=1.0,
        fiber_attenuation_db_per_km=None,
    ),
]

fiber_link = OpticalFiberLink(
    length_km=10.0,
    attenuation_db_per_km=0.35,
    additional_loss_db=0.2,
)

pd = Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9)
rx = Receiver(
    noise_bandwidth_hz=25e9,
    input_current_noise_density_a_per_sqrt_hz=10e-12,
    rin_db_per_hz=-150,
)

results = analyze_wdm_link(
    channels=channels,
    fiber_link=fiber_link,
    pd=pd,
    rx=rx,
    margin_db=0.5,
)
```

各結果の `level_diagram` には、External LaserとMOD入力までは `cw_power_w`／
`cw_power_dbm`、MOD出力以降は `p0_w`、`p1_w`、`pavg_w`、`oma_dbm` が格納
されます。`cumulative_loss_db` はExternal Laserからの光1経路累積損失です。

```python
channel_result = results["lambda_1310"]
for point in channel_result.level_diagram:
    print(point.section, point.name, point.cumulative_loss_db)

print(channel_result.pd_input.oma_dbm)
print(channel_result.receiver_result["ber"])
```

### 13.8 既存OMA APIとの関係

既存の `calculate_ber_from_oma()` は公開APIとして維持されています。内部では
OMA/ERからP0/P1を生成した後、RX境界コアと同じ計算を使用します。

```text
calculate_ber_from_oma(OMA, ER)
  -> nrz_levels_from_oma_er()
  -> P0/P1
  -> receiver BER core
```

新しい統合フローは、レベルダイヤで得たP0/P1を直接渡します。

```text
analyze_laser_to_receiver()
  -> ExternalLaser + OOKModulator
  -> build_level_diagram()
  -> PD-input P0/P1
  -> calculate_ber_from_optical_levels()
```

両方が同じPD入力P0/P1を持つ場合、電流、雑音、Q、BERは一致します。これにより、
従来のOMAベース解析を壊さず、TX/linkとRX単体検証を分離できます。

## 14. 各モデルを組み合わせる際の注意

補助機能は独立しており、同じ物理効果を重複適用しないよう注意が必要です。

例えば、周波数応答を使って波形の振幅低下とISIを計算した後、同じ3 dB帯域に
由来するスカラー帯域ペナルティをさらにOMAから引くと、帯域効果を二重計上する
可能性があります。

同様に、実測TIA出力雑音を入力換算して使う場合、PDシャント熱雑音などがすでに
測定値へ含まれていないか確認する必要があります。

推奨する使い分けは次のとおりです。

- 初期感度見積り: スカラーOMA-to-BERモデル
- 帯域の一次見積り: 1次ローパス・ペナルティ
- パターン依存性の確認: 波形／ISI／アイ補助機能
- ISIとガウス雑音の結合: `oma_ber.time_domain` のPhase 2統計Eye
- 実測または回路応答の確認: 周波数応答CSVからのフィルタリング
- 高入力検討: 飽和sweepと、BER非単調性の確認

## 15. 主な仮定と非対応項目

- 主BERモデルはスカラーかつメモリレスである。
- 0/1サンプル電流をガウス分布で近似する。
- BER式では0と1を等確率とする。
- 判定器は平均電流間の単一しきい値を使う。
- 等価雑音帯域はユーザー入力である。
- TIA雑音密度とRINは指定帯域内で白色と近似する。
- 雑音源間の相関は扱わない。
- PD容量、バイアス、リターンロスは主BERへ自動反映されない。
- 飽和は静的な経験式であり、動的帯域変化を扱わない。
- スカラーBER経路では決定論的ISIと雑音を自動結合しない。Phase 2統計Eyeは別経路
  としてISIと解析的ガウス雑音モーメントを結合する。
- PAM4 BERは実装されていない。
- ランダム雑音サンプルを生成する時間領域Monte Carloは実装されていない。
- CDR、ジッタ、FEC、符号化利得は扱わない。
- 規格固有のマスク、TDECQ、コンプライアンス判定は行わない。

モデルの結果は、入力した仮定の範囲内でのシステムレベル推定として解釈して
ください。設計保証や規格適合性の主張には、測定、回路シミュレーション、デバイス
モデル、対象規格に沿った追加検証が必要です。
