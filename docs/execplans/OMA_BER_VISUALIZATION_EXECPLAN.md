# OMA-to-BER Visualization ExecPlan

## Purpose

Add plotting and visualization scripts for the existing NRZ/OOK scalar Gaussian OMA-to-BER MVP.

The goal is to make receiver behavior easier to inspect without changing the numerical core:

```text
OMA sweep
  -> BER and Q curves
Responsivity sweep
  -> BER and required OMA trends
Noise/source breakdown
  -> shot, TIA, and RIN contribution view
Link-budget helper
  -> transmitter OMA to photodiode-input OMA visualization
```

Plots must use the already implemented deterministic APIs and must not introduce PAM4, time-domain simulation, waveform eye diagrams, or standards-compliance claims.

## Non-goals

Do not implement:

- PAM4 visualization.
- Time-domain waveform or eye simulation.
- S-parameter, ISI, or TDECQ visualizations.
- Interactive dashboard dependencies.
- New numerical models that are not already covered by tests.

## Design Rules

- Core numerical modules should remain deterministic and side-effect free.
- Plotting should live in `src/oma_ber/plotting.py` and/or `examples/`.
- Example scripts may show figures or save optional PNGs when explicitly requested by CLI option.
- Default example behavior should run from repository root.
- Use `matplotlib` only.
- Keep units explicit in axis labels: dBm, W, A, Hz, A/W, BER.
- Use logarithmic y-axis for BER.
- Avoid silently clipping BER; if a plot needs a floor for display, document that it is display-only.
- Do not make standards-compliance claims.

## Proposed Files

```text
src/
  oma_ber/
    plotting.py
examples/
  05_plot_oma_sweep.py
  06_plot_responsivity_sweep.py
  07_plot_noise_breakdown.py
  08_plot_link_budget.py
tests/
  test_plotting.py
```

If the first pass should be smaller, implement only `plotting.py` and `examples/05_plot_oma_sweep.py`.

## Milestone V1: Plotting Helpers

### File

- `src/oma_ber/plotting.py`

### Required functions

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

### Behavior

- `plot_oma_sweep` plots `ber` versus `oma_dbm`.
- `plot_q_vs_oma` plots `q_rx` versus `oma_dbm`.
- `plot_required_oma_vs_responsivity` plots required OMA in dBm versus responsivity in A/W.
- Helpers return the `Axes` object and do not call `plt.show()`.
- Helpers should set clear axis labels and grid.
- Helpers should validate required dictionary keys and raise `ValueError` with a clear message when missing.

### Tests

- Verify each helper returns an `Axes`.
- Verify expected x/y labels are non-empty and include units.
- Verify missing keys raise `ValueError`.
- Use `matplotlib.use("Agg")` in tests to avoid GUI requirements.

## Milestone V2: OMA Sweep Plot Example

### File

- `examples/05_plot_oma_sweep.py`

### Behavior

- Runs from repository root.
- Uses example `Photodiode` and `Receiver` values from the MVP ExecPlan.
- Sweeps OMA over a practical dBm range, for example `-24` to `-8` dBm.
- Produces a two-panel figure:
  - BER versus OMA in dBm.
  - Q versus OMA in dBm.
- Accepts optional CLI arguments:
  - `--save path/to/file.png`
  - `--show`
- If neither option is provided, print a short summary and exit without writing files.

### Checks

```bash
uv run python examples/05_plot_oma_sweep.py
uv run python examples/05_plot_oma_sweep.py --save /tmp/oma_sweep.png
```

The default command must not write files.

## Milestone V3: Responsivity Sweep Plot Example

### File

- `examples/06_plot_responsivity_sweep.py`

### Behavior

- Sweeps `responsivity_a_per_w`.
- For each responsivity, computes `required_oma_dbm` for a target BER.
- Plots required OMA dBm versus responsivity A/W.
- Optional `--save` and `--show` flags follow the same pattern as V2.

### Notes

This is already in the future extension direction, but it uses only existing scalar MVP APIs.

## Milestone V4: Noise Breakdown Plot Example

### File

- `examples/07_plot_noise_breakdown.py`

### Behavior

- For an OMA sweep, compute per-level noise components:
  - shot noise at level 0 and level 1.
  - TIA noise.
  - RIN noise at level 0 and level 1.
  - total noise at level 0 and level 1.
- Plot RMS current noise in A versus OMA in dBm.
- This may require adding a small helper in `noise.py` or `plotting.py`.

### Preferred implementation

Keep core `noise.py` side-effect free. If adding a helper, make it return a plain dictionary:

```python
def noise_components_rms_a(
    optical_power_w: float,
    photocurrent_a: float,
    pd: Photodiode,
    rx: Receiver,
) -> dict[str, float]: ...
```

Add focused tests if this helper is introduced.

## Milestone V5: Link Budget Plot Example

### File

- `examples/08_plot_link_budget.py`

### Behavior

- Demonstrate transmitter OMA dBm minus path losses and margin.
- Plot PD-input OMA and resulting BER for a small set of link-loss cases.
- Keep it explicitly illustrative, not a standards-compliance plot.

## Validation Checklist

Before finishing visualization work, review:

- BER axis uses log scale.
- OMA axis is labeled dBm.
- Current/noise axes are labeled A.
- Responsivity axis is labeled A/W.
- No core numerical function calls `plt.show()`, writes files, or depends on GUI backend.
- Examples run from repository root.
- Optional save paths are only used when `--save` is provided.
- No PAM4 or time-domain simulation slipped in.

## Required Checks

Run:

```bash
uv run pytest
uv run ruff check .
uv run python examples/05_plot_oma_sweep.py
```

If save support is implemented, also run:

```bash
uv run python examples/05_plot_oma_sweep.py --save /tmp/oma_sweep.png
```

Report whether the PNG was created successfully, then remove temporary output if appropriate.

## Recommended First Implementation Slice

Start with:

1. `src/oma_ber/plotting.py`
2. `tests/test_plotting.py`
3. `examples/05_plot_oma_sweep.py`

This gives immediate value while keeping the change small and low risk.
