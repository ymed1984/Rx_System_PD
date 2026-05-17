# CODEX_VISUALIZATION_TASK_PROMPTS.md

# OMA-to-BER Visualization: Codex Execution Prompts

このファイルは、既存の NRZ/OOK scalar Gaussian OMA-to-BER MVP に対して、グラフ・可視化スクリプトを実装するための Codex 実行プロンプト集です。

想定ファイル:
- `AGENTS.md`
- `docs/execplans/OMA_BER_MVP_EXECPLAN.md`
- `docs/execplans/OMA_BER_VISUALIZATION_EXECPLAN.md`

推奨は、まず inspection を行い、その後 small slice で `plotting.py`、`tests/test_plotting.py`、`examples/05_plot_oma_sweep.py` を実装することです。

---

## Prompt 0: Repository inspection only

最初にこれを投げて、Codexにリポジトリ構成・既存API・実行コマンドを確認させます。
この段階ではファイル変更させません。

```text
Read AGENTS.md and inspect this repository.

Also read:
- docs/execplans/OMA_BER_MVP_EXECPLAN.md
- docs/execplans/OMA_BER_VISUALIZATION_EXECPLAN.md

Do not modify files yet.

Summarize:
1. Current Python project structure.
2. Existing oma_ber public APIs relevant to OMA sweeps, BER calculation, required OMA, Photodiode, Receiver, and noise.
3. Available test and lint commands.
4. Whether matplotlib is already available in the project dependencies.
5. The smallest safe implementation slice for visualization.

Keep the proposed first implementation slice limited to:
- src/oma_ber/plotting.py
- tests/test_plotting.py
- examples/05_plot_oma_sweep.py

Do not propose PAM4, time-domain simulation, S-parameter import, eye diagrams, TDECQ, interactive dashboards, or standards-compliance plots.
```

---

## Prompt 1: Recommended first implementation slice

最初の実装用プロンプトです。基本的にはこれを使えばよいです。

```text
Proceed with the first visualization implementation slice according to:
- AGENTS.md
- docs/execplans/OMA_BER_MVP_EXECPLAN.md
- docs/execplans/OMA_BER_VISUALIZATION_EXECPLAN.md

Scope this task strictly to:
1. src/oma_ber/plotting.py
2. tests/test_plotting.py
3. examples/05_plot_oma_sweep.py

Do not implement:
- PAM4 visualization
- time-domain waveform simulation
- eye diagrams
- S-parameter, ISI, or TDECQ visualization
- interactive dashboards
- saturation visualization
- Monte Carlo or process variation plots
- standards-compliance claims

Implementation requirements:

A. Add `src/oma_ber/plotting.py` with the following helpers:

```python
def plot_oma_sweep(
    sweep_results: list[dict],
    ax: matplotlib.axes.Axes | None = None,
) -> matplotlib.axes.Axes: ...

def plot_q_vs_oma(
    sweep_results: list[dict],
    ax: matplotlib.axes.Axes | None = None,
) -> matplotlib.axes.Axes: ...

def plot_required_oma_vs_responsivity(
    responsivity_a_per_w_values: np.ndarray,
    required_oma_dbm_values: np.ndarray,
    ax: matplotlib.axes.Axes | None = None,
) -> matplotlib.axes.Axes: ...
```

Behavior:
- `plot_oma_sweep` plots BER versus OMA in dBm.
- `plot_q_vs_oma` plots `q_rx` versus OMA in dBm.
- `plot_required_oma_vs_responsivity` plots required OMA in dBm versus responsivity in A/W.
- Use matplotlib only.
- Return the `Axes` object.
- Do not call `plt.show()` in library helpers.
- Do not write files in library helpers.
- Use logarithmic y-axis for BER.
- Set explicit axis labels with units:
  - OMA [dBm]
  - BER
  - Q
  - Responsivity [A/W]
  - Required OMA [dBm]
- Add grid lines.
- Validate required dictionary keys and raise `ValueError` with a clear message when keys are missing.
- Do not silently clip BER values. If display-only handling is needed, document it clearly.

B. Add `tests/test_plotting.py`.

Test requirements:
- Use `matplotlib.use("Agg")` before importing pyplot.
- Verify each plotting helper returns a matplotlib `Axes`.
- Verify axis labels are non-empty and include units where relevant.
- Verify `plot_oma_sweep` uses a logarithmic y-axis.
- Verify missing required keys in sweep result dictionaries raise `ValueError`.
- Keep tests deterministic and independent of GUI backend.

C. Add `examples/05_plot_oma_sweep.py`.

Example behavior:
- Must run from repository root.
- Use existing MVP APIs and example Photodiode/Receiver values.
- Sweep OMA over a practical range, for example -24 dBm to -8 dBm.
- Produce a two-panel figure:
  1. BER versus OMA [dBm]
  2. Q versus OMA [dBm]
- Support CLI options:
  - `--save path/to/file.png`
  - `--show`
- If neither `--save` nor `--show` is provided, print a short summary and exit without writing files.
- The default command must not create or modify output files.
- When `--save` is provided, create parent directories if needed.
- Keep the example illustrative; do not make standards-compliance claims.

Validation:
Run the following commands if available:

```bash
uv run pytest
uv run ruff check .
uv run python examples/05_plot_oma_sweep.py
uv run python examples/05_plot_oma_sweep.py --save /tmp/oma_sweep.png
```

If the project does not use uv, use the repository's existing Python/test commands instead.

After saving `/tmp/oma_sweep.png`, confirm that the file exists. Remove the temporary PNG if appropriate.

At the end, summarize:
- files changed
- functions added
- examples added
- tests added
- commands run
- whether the PNG save check passed
- any limitations or assumptions
```

---

## Prompt 2: Add responsivity sweep visualization

`required_oma_dbm()` が既に実装されていて、最初の可視化が安定した後に使います。

```text
Add the responsivity sweep visualization according to docs/execplans/OMA_BER_VISUALIZATION_EXECPLAN.md.

Scope this task to:
- examples/06_plot_responsivity_sweep.py
- small tests only if needed
- no changes to core numerical models unless required for API compatibility

Requirements:
- Sweep `responsivity_a_per_w` over a reasonable range, for example 0.3 to 1.1 A/W.
- For each responsivity, compute `required_oma_dbm` for a target BER using existing MVP APIs.
- Plot required OMA [dBm] versus Responsivity [A/W].
- Use the existing `plot_required_oma_vs_responsivity()` helper.
- Support CLI options:
  - `--save path/to/file.png`
  - `--show`
  - `--target-ber`, default 1e-12
- If neither `--save` nor `--show` is provided, print a short summary and exit without writing files.
- Do not implement PAM4, waveform simulation, S-parameter import, saturation, Monte Carlo, or standards-compliance plots.

Validation:
Run:
```bash
uv run pytest
uv run ruff check .
uv run python examples/06_plot_responsivity_sweep.py
uv run python examples/06_plot_responsivity_sweep.py --save /tmp/responsivity_sweep.png
```

If the project does not use uv, use the repository's existing commands.

Summarize changed files, commands run, and assumptions.
```

---

## Prompt 3: Add noise breakdown visualization

これは、shot/TIA/RINの寄与を見たい段階で使います。
必要なら `noise_components_rms_a()` を追加しますが、coreはside-effect freeに保ちます。

```text
Add the noise breakdown visualization according to docs/execplans/OMA_BER_VISUALIZATION_EXECPLAN.md.

Scope this task to:
- a small deterministic noise component helper if needed
- examples/07_plot_noise_breakdown.py
- focused tests for any new helper

Requirements:
- For an OMA sweep, compute RMS current noise components:
  - shot noise at level 0
  - shot noise at level 1
  - TIA noise
  - RIN noise at level 0
  - RIN noise at level 1
  - total noise at level 0
  - total noise at level 1
- Plot RMS current noise [A] versus OMA [dBm].
- Use matplotlib only.
- Support CLI options:
  - `--save path/to/file.png`
  - `--show`
- If neither `--save` nor `--show` is provided, print a short summary and exit without writing files.

If adding a helper, prefer this API:
```python
def noise_components_rms_a(
    optical_power_w: float,
    photocurrent_a: float,
    pd: Photodiode,
    rx: Receiver,
) -> dict[str, float]: ...
```

Rules:
- Keep core `noise.py` side-effect free.
- Do not call matplotlib from `noise.py`.
- Do not change existing BER results.
- Do not introduce new numerical models beyond the existing shot, TIA, and RIN components.
- Do not implement PAM4, time-domain simulation, S-parameters, eye diagrams, or TDECQ.

Validation:
Run:
```bash
uv run pytest
uv run ruff check .
uv run python examples/07_plot_noise_breakdown.py
uv run python examples/07_plot_noise_breakdown.py --save /tmp/noise_breakdown.png
```

If the project does not use uv, use the repository's existing commands.

Summarize changed files, commands run, and any assumptions.
```

---

## Prompt 4: Add link-budget visualization

Tx OMAからPD入力OMAへの損失・マージンを見せる段階で使います。

```text
Add the illustrative link-budget visualization according to docs/execplans/OMA_BER_VISUALIZATION_EXECPLAN.md.

Scope this task to:
- examples/08_plot_link_budget.py
- small plotting/helper tests only if needed

Requirements:
- Demonstrate transmitter OMA [dBm] minus path losses and margin.
- Plot PD-input OMA [dBm] and resulting BER for a small set of illustrative link-loss cases.
- Use existing link-budget and BER APIs if available.
- If a link-budget helper already exists, reuse it.
- If no helper exists, keep any new helper minimal, deterministic, and well tested.
- Support CLI options:
  - `--save path/to/file.png`
  - `--show`
- If neither `--save` nor `--show` is provided, print a short summary and exit without writing files.
- Keep the plot explicitly illustrative.
- Do not make standards-compliance claims.
- Do not implement PAM4, time-domain waveforms, S-parameters, TDECQ, or interactive dashboards.

Validation:
Run:
```bash
uv run pytest
uv run ruff check .
uv run python examples/08_plot_link_budget.py
uv run python examples/08_plot_link_budget.py --save /tmp/link_budget.png
```

If the project does not use uv, use the repository's existing commands.

Summarize changed files, commands run, and assumptions.
```

---

## Prompt 5: Review visualization implementation

実装後のレビュー用です。

```text
Review the visualization implementation against:
- AGENTS.md
- docs/execplans/OMA_BER_VISUALIZATION_EXECPLAN.md

Focus on:
1. Whether plotting helpers are side-effect free.
2. Whether library functions avoid `plt.show()` and file writes.
3. Whether examples only write files when `--save` is provided.
4. Whether default example execution produces no files.
5. Whether BER plots use log y-axis.
6. Whether axis labels include units.
7. Whether missing dictionary keys raise clear `ValueError`.
8. Whether tests use the Agg backend.
9. Whether no PAM4, time-domain simulation, S-parameter, eye diagram, TDECQ, dashboard, saturation, or Monte Carlo features were introduced.
10. Whether the implementation changed numerical BER results unintentionally.

Run:
```bash
uv run pytest
uv run ruff check .
uv run python examples/05_plot_oma_sweep.py
uv run python examples/05_plot_oma_sweep.py --save /tmp/oma_sweep.png
```

If available, also run the other visualization examples.

Report:
- issues found
- fixes applied
- commands run
- remaining limitations
```

---

## One-shot master prompt

慣れてきた後、1回でまとめて実装したい場合のプロンプトです。
ただし、初回は Prompt 0 → Prompt 1 の2段階を推奨します。

```text
Read AGENTS.md and the following ExecPlans:
- docs/execplans/OMA_BER_MVP_EXECPLAN.md
- docs/execplans/OMA_BER_VISUALIZATION_EXECPLAN.md

First inspect the repository structure, existing Python APIs, test command, lint command, and dependency setup.

Then implement the first visualization slice only:
1. src/oma_ber/plotting.py
2. tests/test_plotting.py
3. examples/05_plot_oma_sweep.py

Required plotting helpers:
- plot_oma_sweep()
- plot_q_vs_oma()
- plot_required_oma_vs_responsivity()

Requirements:
- Use matplotlib only.
- Return Axes from plotting helpers.
- Do not call plt.show() in library helpers.
- Do not write files in library helpers.
- Use log y-axis for BER.
- Use explicit units in axis labels.
- Validate required dictionary keys and raise clear ValueError on missing keys.
- Use matplotlib Agg backend in tests.
- Example script supports --save and --show.
- Default example execution prints a short summary and writes no files.
- Save mode creates a PNG only when --save is provided.

Do not implement:
- PAM4 visualization
- time-domain waveform simulation
- eye diagrams
- S-parameter, ISI, or TDECQ visualization
- interactive dashboards
- saturation visualization
- Monte Carlo or process variation plots
- standards-compliance claims

Validation:
Run:
```bash
uv run pytest
uv run ruff check .
uv run python examples/05_plot_oma_sweep.py
uv run python examples/05_plot_oma_sweep.py --save /tmp/oma_sweep.png
```

If uv is not used in this repo, use the existing test/lint/Python commands.

Fix failures when safe. If a failure cannot be resolved safely, stop and explain the blocker.

At the end, summarize:
- files changed
- public functions added
- example scripts added
- tests added
- validation commands run
- whether the PNG save check passed
- limitations and next recommended visualization milestone
```
