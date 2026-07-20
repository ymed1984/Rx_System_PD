# Eye Diagram Phase 3以降の実装計画

## 1. 目的

この文書は、Phase 2まで実装済みのNRZ/OOK受信器Eye解析を基準に、次回以降の
実装をそのまま開始できる粒度で定義する実行計画です。

この計画では、Phase 3以降を次の順序で進めます。

```text
Phase 2（現在の確定境界）
  -> Phase 3: 実測／回路由来の複素PD・TIA応答を物理単位付きで統合
  -> Phase 4: 残留タイミングジッタとbathtub解析
  -> Phase 5: ジッタPSDとCDR追従／残差伝達
  -> Phase 6: Monte Carloによる解析モデルの限定的な相互検証
  -> Phase 7: 必要な場合のみ光電界・分散・chirp・WDM干渉へ拡張
```

Phase 3とPhase 4は実装済みです。次回はPhase 4の残留ジッタ出力を境界条件として、
「6. Phase 5」のジッタPSDとCDR追従／残差伝達から進めます。

## 2. Phase 2の確定済みベースライン

Phase 1とPhase 2で、次の経路が実装済みです。

```text
PD-input P0/P1 [W]
  -> NRZ/OOK optical power waveform [W]
  -> static PD conversion / optional saturation [A]
  -> causal H_PD(z) [A/A]
  -> causal Z_TIA(z) [V/A]
  -> deterministic TIA waveform [V]
  -> source-specific one-sided noise PSD propagation
  -> pattern-dependent variance and lag-1 covariance
  -> Gaussian-mixture statistical eye
  -> BER-optimum voltage threshold and sampling phase
```

Phase 2が扱う雑音源は次のとおりです。

```text
S_photo(t) = 2 q i_signal(t)                         [A^2/Hz]
S_dark     = 2 q F_dark I_dark                       [A^2/Hz]
S_TIA      = i_n^2                                   [A^2/Hz]
S_RIN(t)   = (Rpd P(t))^2 RIN_linear                 [A^2/Hz]
S_thermal  = 4 k_B T / R_shunt                       [A^2/Hz]
```

発生ノードに応じて、PD起因雑音は `H_PD*Z_TIA`、TIA入力換算雑音は `Z_TIA` のみを
通過します。現在のベースラインチェックは次のとおりです。

```text
uv run pytest       -> 285 passed
uv run ruff check . -> All checks passed
```

次回作業開始時にも、変更前にこの2コマンドを実行して基準状態を確認します。

## 3. 再計画で採用する優先順位

### 3.1 Phase 3は実測複素応答統合を優先する

Phase 3では、ジッタやCDRより先に、PDとTIAの複素周波数応答を新しい
`oma_ber.time_domain` 経路へ統合します。

理由は次のとおりです。

- ジッタBERは判定点付近の波形勾配 `dv/dt` に強く依存する。
- `dv/dt`、ringing、pre/post-cursor ISI、群遅延は、1次LPFでは十分に表現できない。
- Phase 2の信号平均、雑音分散、隣接共分散は同じインパルス応答から一貫して計算
  する必要がある。
- 実応答統合前にジッタを追加すると、仮の1次LPFに対するジッタ感度を精密化する
  ことになり、デバイス設計への寄与が限定的になる。

したがって、Phase 3の完了条件は「実測応答を表示できること」ではなく、同じ応答が
決定論的波形と雑音PSD伝搬の両方へ使われ、1次LPFとの比較が可能になることです。

### 3.2 CDRはスカラーRJ値から推定しない

CDRが追従できるかどうかは、総RMSジッタだけでなくジッタの周波数成分に依存します。
このため、Phase 4ではCDR後の残留ジッタを入力し、Phase 5で初めてジッタPSDとCDR
ループ応答を結合します。

```text
禁止する近似:
  residual_jitter_rms = arbitrary_factor * input_jitter_rms

採用する方向:
  S_t,residual(f) = |H_error(f)|^2 S_t,input(f)
```

CDR帯域だけから説明なしの係数を作らないことを実装上のガードレールとします。

## 4. Phase 3: 実測／回路複素応答の統合

### 4.1 スコープ

Phase 3で実装するものは次のとおりです。

- CSVによる複素周波数応答入力。
- PD電流応答 `H_PD(f)` [A/A] とTIAトランスインピーダンス
  `Z_TIA(f)` [V/A = ohm] の明示的な区別。
- 周波数グリッド、DC、Nyquist、位相、遅延、因果性の診断。
- 複素応答から因果的FIR `DiscreteTransferFunction` への変換。
- 変換したFIRをPhase 1の信号経路とPhase 2のPSD経路で共用。
- 1次LPFと実測応答の決定論的Eye、ENBW、Gaussian-mixture BER比較。

Phase 3で実装しないものは次のとおりです。

- Touchstone全仕様対応。
- 自動de-embedding、fixture除去、reference impedance変換。
- 非因果応答の自動修復。
- passivity enforcement。
- 光Sパラメータと電気Sパラメータの自動判別。
- ジッタ、CDR、Monte Carlo雑音波形。

### 4.2 既存実装を直接流用しない理由

既存の `src/oma_ber/sparameters.py` は互換用の簡易ヘルパーです。現在は、

- 周波数が0より大きいことを要求する。
- DC以下を最初の値で保持する。
- 測定上限からNyquistまで最後の複素値を保持する。
- 実部と虚部を別々に線形補間する。
- `irfft` 後の因果性、時間原点、pre-echoを検証しない。
- 応答単位と物理基準面を保持しない。

という境界があります。Phase 3では既存APIを破壊せず、新しい物理経路用APIを
`src/oma_ber/time_domain/measured.py` に追加します。

### 4.3 提案データモデル

```python
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

MeasuredResponseKind = Literal["dimensionless", "transimpedance_ohm"]


@dataclass(frozen=True)
class MeasuredFrequencyResponse:
    frequency_hz: np.ndarray
    response_complex: np.ndarray
    response_kind: MeasuredResponseKind
    reference_plane: str
    source_path: Path | None = None


@dataclass(frozen=True)
class MeasuredResponseDiagnostics:
    dc_gain: float
    equivalent_noise_bandwidth_hz: float
    group_delay_s: float
    pre_echo_energy_ratio: float
    discarded_tail_energy_ratio: float
    fft_size: int
    fir_length_samples: int
```

`reference_plane` は空文字を禁止します。例は次のように明示します。

```text
PD response:  "PD optical conversion output -> TIA input"
TIA response: "TIA input current -> TIA output voltage"
```

### 4.4 提案API

```python
def measured_frequency_response_from_csv(
    path: str | Path,
    response_kind: MeasuredResponseKind,
    reference_plane: str,
) -> MeasuredFrequencyResponse: ...


def fir_transfer_from_measured_response(
    response: MeasuredFrequencyResponse,
    sample_rate_hz: float,
    fft_size: int,
    fir_length_samples: int,
    removed_reference_delay_s: float = 0.0,
    causality_tolerance: float = 1e-6,
    tail_energy_tolerance: float = 1e-6,
) -> tuple[DiscreteTransferFunction, MeasuredResponseDiagnostics]: ...
```

CSV形式は既存形式を踏襲します。

```text
frequency_hz,magnitude_db,phase_deg
0.0,63.52,0.0
...
```

ただし、Phase 3経路では次を必須とします。

- `frequency_hz` は0を含み、厳密な昇順。
- 測定上限は `sample_rate_hz/2` 以上。
- `magnitude_db` は電流／電圧の振幅比なので `20*log10` 基準。
- `response_kind="dimensionless"` では線形振幅利得。
- `response_kind="transimpedance_ohm"` では線形値の単位がohm。
- DCとNyquistの複素値は、実波形のHermitian条件に整合する実数値。
- NaN、Inf、重複周波数を拒否。

DCまたはNyquistを暗黙外挿しません。不足時は `ValueError` とし、ユーザーが物理的な
境界値をCSVへ追加します。

### 4.5 複素応答からFIRへの変換

ターゲット周波数グリッドは次式です。

```text
f_k = k * fs / Nfft
k = 0 ... Nfft/2
```

変換手順は次のとおりです。

1. 振幅を線形値、位相をradへ変換する。
2. 位相をunwrapする。
3. 指定された既知の基準遅延だけを除去する。

```text
H_deembedded(f) = H_measured(f) * exp(+j 2 pi f tau_removed)
```

4. 振幅とunwrap位相をターゲットrFFTグリッドへ補間する。
5. `irfft` で実インパルス応答を得る。
6. 負時間側に相当する周期末尾のenergyを `pre_echo_energy_ratio` として診断する。
7. 先頭から `fir_length_samples` で切り出し、捨てた正時間tail energyを診断する。
8. 許容値を超えるpre-echoまたはtailがあれば自動補正せず `ValueError` にする。
9. FIR係数を `DiscreteTransferFunction(numerator=h, denominator=[1])` へ変換する。

```text
E_total = sum_n |h[n]|^2
pre_echo_energy_ratio = E_negative / E_total
discarded_tail_energy_ratio = E_discarded_positive / E_total
```

`removed_reference_delay_s` は測定fixtureや伝送路の既知の純遅延を除去するための入力
です。最大tapを自動で0へ移動する処理は禁止します。自動移動はpre-cursorや反射を
純遅延と誤認するためです。

### 4.6 信号応答と雑音応答の共用

作成したFIRは既存の次の引数へそのまま渡します。

```python
waveforms = simulate_pd_tia_waveform(
    ...,
    pd_current_response=measured_pd_transfer,
    tia_transimpedance_response=measured_tia_transfer,
)

noise = calculate_tia_output_noise(
    waveforms,
    noise_model,
)
```

これにより、信号平均は `h[n]`、白色入力雑音分散は `h[n]^2`、隣接共分散は
`h[n]h[n+1]` から計算されます。同じ応答を手作業で二重入力しません。

### 4.7 Phase 3のファイル構成

```text
src/oma_ber/time_domain/measured.py
src/oma_ber/time_domain/__init__.py
src/oma_ber/plotting.py
tests/test_time_domain_measured.py
tests/data/synthetic_pd_response.csv
tests/data/synthetic_tia_response.csv
examples/16_measured_response_statistical_eye.py
README.md
docs/MechanismExplanation.md
docs/explanation.md
```

実測データをテストfixtureとして登録する場合は、出所、測定条件、再配布可否を
確認します。確認できないデータは入れず、既知の解析式から作ったsynthetic response
を使います。

### 4.8 Phase 3のテスト

最低限、次を実装します。

1. 0 Hzを含まないCSVを拒否する。
2. Nyquistまで測定範囲がない応答を拒否する。
3. 応答単位と `response_kind` が保持される。
4. flatなdimensionless応答がidentity FIRへ近づく。
5. flatなtransimpedance応答のDC gainが指定ohm値と一致する。
6. synthetic 1次LPFの決定論的波形が解析1次LPFと一致する。
7. synthetic 1次LPFのENBWが解析1次LPFと一致する。
8. 純遅延を既知値で除去した結果が遅延なし応答と一致する。
9. 許容値を超えるpre-echoを拒否する。
10. 許容値を超えるFIR tail切捨てを拒否する。
11. PD応答とTIA応答の `response_kind` 取り違えを拒否する。
12. 実測FIR使用時もPhase 2のGaussian-mixture BERが有限値を返す。
13. plotting helperが1次LPFと実測応答のEye／BER比較を表示できる。
14. 全既存テストが回帰しない。

### 4.9 Phase 3の完了条件

```text
[x] 新APIが既存sparameters.py APIを破壊しない
[x] PDとTIAの単位が明示されている
[x] DCとNyquistに暗黙外挿がない
[x] pre-echoと切捨てenergyが数値化される
[x] 同じFIRが信号と雑音PSDへ適用される
[x] 1次LPFとの決定論的Eye、ENBW、BER比較が可能
[x] example 16がrepository rootから実行可能
[x] READMEと物理説明資料が更新されている
[x] uv run pytest が全件合格（Phase 3完了時点: 298 tests）
[x] uv run ruff check . が合格
```

## 5. Phase 4: 残留タイミングジッタとbathtub解析

### 5.1 前提

Phase 4はPhase 3完了後に開始します。ただし実測応答が入手できない場合でも、
synthetic複素応答でPhase 3の変換・診断・統合テストが完了していれば開始できます。

Phase 4が受け取るジッタは、CDR後または外部から与えられた「判定器に対する残留
タイミング誤差」です。CDR入力ジッタではありません。

### 5.2 最初に扱うジッタ成分

```python
@dataclass(frozen=True)
class ResidualTimingJitter:
    random_jitter_rms_s: float = 0.0
    sinusoidal_jitter_peak_s: float = 0.0
```

- `random_jitter_rms_s`: 平均0のGaussian RJ標準偏差 [s]。
- `sinusoidal_jitter_peak_s`: 位相が一様な正弦ジッタのpeak時間偏差 [s]。

Dual-DiracのDJ、DCD、data-dependent jitterは同時に追加しません。まずRJと単一SJを
解析的に検証し、必要性が確認された後に別パラメータとして追加します。

### 5.3 ジッタ平均BER

nominal sampling phaseを `phi0`、残留時間誤差を `delta_t` とすると、実際の位相は
次式です。

```text
phi = phi0 + delta_t / T_UI
```

固定電圧しきい値 `gamma` に対するジッタ平均BERは、

```text
BER_J(gamma, phi0)
  = integral BER(gamma, phi0 + delta_t/T_UI)
             p(delta_t) d(delta_t)
```

です。重要な制約として、`gamma` をジッタの瞬時値ごとに再最適化してはいけません。
各nominal phaseに対して1つの固定しきい値を最適化します。

```text
gamma_opt(phi0) = argmin_gamma BER_J(gamma, phi0)
phi_opt = argmin_phi0 BER_J(gamma_opt(phi0), phi0)
```

RJ積分はGauss-Hermite求積、SJの一様位相積分は周期台形求積を使います。RJとSJを
同時指定した場合は積求積にします。乱数seedや有限ランダム波形に依存させません。

### 5.4 UI境界を跨ぐサンプル

ジッタで `phi < 0` または `phi >= 1` となる場合、単純に `phi % 1` として同じ
bit labelへ付け直してはいけません。判定対象bit kは固定したまま、サンプル時刻だけ
が隣接UIへ移動します。

```text
t_sample(k) = k*T_UI + phi0*T_UI + delta_t
decision label = bit[k]
```

有限PRBSは周期波形として延長し、波形と雑音分散の両方を周期補間します。これを
Phase 4実装前の必須リファクタリングとします。

### 5.5 Phase 4の提案API

```python
def analyze_jittered_tia_eye(
    waveforms: ReceiverWaveformResult,
    noise: TiaOutputNoiseResult,
    jitter: ResidualTimingJitter,
    nominal_phases_ui: np.ndarray | None = None,
    decision_delay_s: float = 0.0,
    random_quadrature_points: int = 17,
    sinusoidal_phase_points: int = 32,
) -> JitteredEyeAnalysis: ...


def plot_bathtub_curve(
    analysis: JitteredEyeAnalysis,
    ax: Axes | None = None,
) -> Axes: ...
```

結果には最低限、次を保持します。

```text
nominal_phases_ui
optimized_threshold_v_per_phase
jitter_averaged_ber_per_phase
optimum_nominal_phase_ui
optimum_threshold_v
optimum_ber
random_jitter_rms_ui
sinusoidal_jitter_peak_ui
```

### 5.6 Phase 4のテストと完了条件

1. RJ=0、SJ=0でPhase 2 BERとしきい値へ一致する。
2. 完全にflatなメモリレス波形ではタイミングジッタがBERを変えない。
3. 遷移を持つ帯域制限波形ではRJ増加により最良BERが悪化する。
4. `random_jitter_rms_s` と `random_jitter_rms_ui*T_UI` が一致する。
5. UI境界を跨いでも判定対象bit labelが変わらない。
6. Gauss-Hermite点数を増やしたとき結果が収束する。
7. SJのみの積分が高分解能の独立参照積分と一致する。
8. 反転TIAでもBERとしきい値極性が正しい。
9. thresholdは瞬時ジッタごとでなくnominal phaseごとに1つである。
10. bathtubの左右端でBERが悪化する既知条件を確認する。

Phase 4はbathtub曲線を提供しますが、規格固有のTJ@BER外挿やdual-Dirac準拠を主張
しません。

### 5.7 Phase 4実装状況

次を実装済みです。

- `ResidualTimingJitter` による残留RJ RMS [s]とSJ peak [s]入力。
- Gauss-Hermite RJ求積と周期SJ位相求積。
- UI境界を越える周期fractional samplingと固定bit label。
- timing×patternの重み付きGaussian mixture。
- nominal phaseごとに1つの固定電圧しきい値を使うBER最適化。
- `JitteredEyeAnalysis` とBER bathtub。
- `plot_jittered_ber_bathtub()`。
- `examples/17_residual_jitter_bathtub.py`。

Phase 4の完了チェックは次のとおりです。

```text
[x] RJ=0、SJ=0でPhase 2と一致
[x] UI境界を越えても判定対象bit labelを保持
[x] Gauss-Hermite点数とSJ位相点数の収束テスト
[x] 反転TIAに対応
[x] 瞬時ジッタごとにthresholdを変更しない
[x] bathtub plottingとExample 17
[x] bathtub中央より左右UI端でBERが悪化する回帰テスト
[x] 全pytest合格（Phase 4完了時点: 310 tests）
[x] ruff check合格
```

## 6. Phase 5: ジッタPSDとCDR

### 6.1 入力モデル

CDR解析には時間誤差PSDを明示的に入力します。

```python
@dataclass(frozen=True)
class TimingJitterPsd:
    frequency_hz: np.ndarray
    one_sided_time_error_psd_s2_per_hz: np.ndarray


@dataclass(frozen=True)
class CdrJitterTransfer:
    frequency_hz: np.ndarray
    tracking_response_complex: np.ndarray
```

`tracking_response_complex = H_track(f)` とすると、判定器から見た残差は、

```text
H_error(f) = 1 - H_track(f)
S_t,residual(f) = |H_error(f)|^2 S_t,input(f)
sigma_t,residual^2 = integral_0^infinity S_t,residual(f) df
```

です。正弦ジッタ周波数 `f_sj` では、

```text
A_residual = A_input * |H_error(f_sj)|
```

を使います。求めた残留RJ/SJをPhase 4へ渡します。

### 6.2 CDRの段階的実装

1. 任意の測定／回路由来 `H_track(f)` の入力と残差積分。
2. 明示的な1次tracking LPFモデル。
3. 必要な場合のみ2次Type-IIモデルへ拡張。
4. VCO／DCO固有位相雑音は別PSD入力として加算。

1次モデルを追加する場合も、近似であることをクラス名とdocstringで明示します。

```text
H_track(f) = 1 / (1 + j f/f_cdr)
H_error(f) = j f/f_cdr / (1 + j f/f_cdr)
```

テストでは `f << f_cdr` で入力ジッタを追従し、`f >> f_cdr` で追従しない極限を
確認します。

### 6.3 Phase 5で扱わないもの

- bang-bang CDRの非線形limit cycle。
- cycle slip。
- data-dependent phase detector gain。
- SSCを含むクロックアーキテクチャ全体。
- 規格固有のjitter tolerance mask自動判定。

これらが必要な場合は、位相領域の線形モデルと分離したイベント駆動モデルとして
追加します。

## 7. Phase 6: Monte Carlo相互検証

### 7.1 目的

Monte CarloはPhase 2～5の置換ではなく、解析的分散、Gaussian-mixture BER、
ジッタ求積を検証する補助経路とします。

### 7.2 最小スコープ

- seedを必須にした再現可能なGaussian白色雑音生成。
- 雑音源ごとの発生ノードからのフィルタリング。
- RJとSJを加えたfractional resampling。
- サンプル平均、分散、lag-1共分散の解析値との比較。
- BERが十分観測できる `1e-2 ... 1e-5` 程度での解析BER比較。

極低BERを通常Monte Carloで検証しません。必要な場合はimportance samplingを独立
Milestoneとして設計します。

### 7.3 合格基準

```text
mean error                  < statistical tolerance
variance relative error     < tolerance based on sample count
lag-1 covariance error      < tolerance based on sample count
BER confidence interval     contains analytic BER
same seed                   -> identical result
different seed              -> statistically independent result
```

許容差はサンプル数から導出し、固定の都合のよい絶対値にしません。

## 8. Phase 7: 光電界・分散・WDM効果

次の要求が具体化した場合のみ実装します。

- External Laser linewidth／phase noise。
- modulator chirp。
- AWGの複素光振幅応答とcrosstalk。
- ファイバ色分散。
- WDM beat noise／隣接チャネル漏れ。
- PD総入力パワーによる共有飽和。

この段階では光パワー `P(t)` だけでなく複素光電界 `E(t)` を導入します。

```text
P(t) = |E(t)|^2
E_RX(f) = E_TX(f) H_AWG(f) H_fiber(f)
i_PD(t) = Rpd |sum_channels E_k(t)|^2
```

現在の独立WDMチャネル損失モデルと混在させず、光電界経路を別モジュールとして
実装します。偏波、コヒーレンス、レーザー間周波数差を定義しないままbeat項を
追加してはいけません。

## 9. 次回の実行手順

次回「Phase 5を開始」と指示された場合は、次の順序で作業します。

### Step 0: ベースライン固定

```bash
git status --short
uv run pytest
uv run ruff check .
```

- Phase 4までの未コミット変更を破壊しない。
- 既存の310テストが通ることを確認する。
- Phase 5と無関係なdirty fileを変更しない。

### Step 1: ジッタPSDとCDR応答のデータ契約

- `TimingJitterPsd` と `CdrJitterTransfer` を実装する。
- 周波数 [Hz]、片側時間誤差PSD [s^2/Hz]、複素追従応答を明示する。
- DC、厳密昇順、有限値、周波数範囲、補間範囲を先に検証する。

### Step 2: CDR誤差伝達と残差積分

- `H_error(f) = 1 - H_track(f)` を同じ周波数グリッドで評価する。
- `|H_error(f)|^2*S_t,input(f)` を積分し、残留RJ RMS [s]を返す。
- 低周波追従、高周波非追従、全追従／無追従の極限を検証する。

### Step 3: Phase 4への明示的な接続

- PSD積分で得た残留RMSを `ResidualTimingJitter` へ渡す補助APIを追加する。
- 入力PSD、残差PSD、積分帯域、残留RMS [s/UI]を結果へ保持する。
- CDR入力ジッタとPhase 4入力の残留ジッタを同じ値として二重計上しない。

### Step 4: 可視化とexample

- 入力ジッタPSD、CDR追従／誤差伝達、残差PSDを比較するplot helperを追加する。
- `examples/18_cdr_residual_jitter.py` を追加する。
- 1次tracking LPFは近似モデルであり規格準拠CDRではないことを明記する。

### Step 5: 文書と全検証

- READMEへ入力PSDと実行方法を追記する。
- `docs/MechanismExplanation.md` へCDR誤差伝達とPSD積分の式を追記する。
- `docs/explanation.md` へPhase 4残留ジッタ入力との接続を追記する。
- 次を実行する。

```bash
uv run pytest
uv run ruff check .
uv run python examples/18_cdr_residual_jitter.py
```

## 10. 実装時のレビュー項目

Phase 3以降の各PRまたは作業単位で、次を確認します。

### 単位と基準面

- `dimensionless` と `transimpedance_ohm` を混同していないか。
- magnitude dBに20、光パワーdBに10を使っているか。
- input/output reference planeが文書化されているか。
- 秒、UI、Hz、baudを暗黙変換していないか。

### 周波数応答

- DCとNyquistを暗黙外挿していないか。
- phase unwrapと純遅延の符号が正しいか。
- `irfft` の周期末尾energyを無視していないか。
- truncationでDC gainまたはENBWが変化していないか。
- 信号用応答と雑音用応答が同一か。

### BERとタイミング

- 0/1のQ関数符号が正しいか。
- 反転TIAの判定方向が正しいか。
- thresholdを瞬時ジッタごとに再最適化していないか。
- UI境界を跨いだときbit labelを誤ってwrapしていないか。
- 極低BERの比較を通常確率値だけで行っていないか。

### モデル境界

- 実測応答を規格準拠モデルと過大表現していないか。
- CDR入力とCDR後残留ジッタを混同していないか。
- Gaussian analytic eyeとMonte Carlo波形を混同していないか。
- 同じ帯域／損失／雑音を二重計上していないか。

## 11. 保留する設計判断

次の項目は、入力データまたは用途が決まるまで実装しません。

1. Touchstone parserを自作するか `scikit-rf` を追加するか。
2. 実測TIA応答をZ21、S21＋reference impedance、直接transimpedanceのどれで受けるか。
3. fixture／probe／packageのde-embedding責任範囲。
4. 2次CDRのloop order、damping、phase-detector gain。
5. dual-Dirac DJと規格固有TJ外挿。
6. optical field modelの偏波・coherence・laser spacing条件。

これらは推測で決めず、具体的な測定形式または検証目的が提示された時点で、この
ExecPlanへ入力契約と合格基準を追記してから実装します。

## 12. 最終的な推奨

Phase 3の実測複素応答統合とPhase 4のCDR後残留ジッタbathtubは完了しました。
次はPhase 5でジッタPSDとCDR追従を接続します。入力ジッタPSDから残差PSDと残留RMS
を求め、その値だけをPhase 4へ渡すことで、CDR効果と残留ジッタを二重計上せず、
1次LPF、実測応答、残留ジッタ、CDRの寄与を個別に比較できます。
