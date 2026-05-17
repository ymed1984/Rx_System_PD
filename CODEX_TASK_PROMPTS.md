# Codex Task Prompts for OMA-to-BER / PD Analysis Tool

Use these prompts one milestone at a time. Smaller tasks are easier to review and less likely to drift.

## Prompt 0: Repository inspection

```text
Read AGENTS.md and inspect this repository. Summarize the current project structure, Python environment, and available test/lint commands. Do not modify files yet. Then propose the smallest implementation path for the OMA-to-BER MVP described in docs/execplans/OMA_BER_MVP_EXECPLAN.md.
```

## Prompt 1: Create package skeleton

```text
Read AGENTS.md and docs/execplans/OMA_BER_MVP_EXECPLAN.md. Implement Milestone 0 only: create the minimal Python package skeleton for the OMA-to-BER tool, including pyproject.toml if needed, src/oma_ber/__init__.py, README.md placeholder, and test directory. Keep changes minimal. Run the available import/test command and report results.
```

## Prompt 2: Implement units and NRZ levels

```text
Read AGENTS.md and docs/execplans/OMA_BER_MVP_EXECPLAN.md. Implement Milestone 1 and Milestone 2 only: units.py, modulation.py, and their tests. Use explicit physical units in names and docstrings. Run pytest and ruff if available. Report changed files, tests run, and any assumptions.
```

## Prompt 3: Implement PD, receiver, and noise model

```text
Read AGENTS.md and docs/execplans/OMA_BER_MVP_EXECPLAN.md. Implement Milestone 3 and Milestone 4 only: Photodiode, Receiver, and noise calculation functions with tests. Validate invalid physical inputs. Do not implement BER or sweeps yet. Run pytest and ruff if available.
```

## Prompt 4: Implement BER model

```text
Read AGENTS.md and docs/execplans/OMA_BER_MVP_EXECPLAN.md. Implement Milestone 5 only: qfunc, threshold BER, optimum threshold, Q estimate, and Gaussian-level BER dictionary. Pay special attention to threshold sign conventions. Add tests for monotonicity and equal-noise threshold behavior. Run pytest and ruff if available.
```

## Prompt 5: Implement top-level calculation and OMA sweep

```text
Read AGENTS.md and docs/execplans/OMA_BER_MVP_EXECPLAN.md. Implement Milestone 6 only: calculate_ber_from_oma, sweep_oma, and required_oma_dbm with tests. Ensure BER decreases with OMA and required_oma_dbm fails clearly when the target BER is out of range. Run pytest and ruff if available.
```

## Prompt 6: Add link-budget helper and examples

```text
Read AGENTS.md and docs/execplans/OMA_BER_MVP_EXECPLAN.md. Implement Milestone 7 and Milestone 8 only: link_budget.py, tests, and executable examples. Examples should run from the repository root. Do not add PAM4 or time-domain simulation. Run pytest, ruff, and at least one example script if possible.
```

## Prompt 7: Complete README and perform review

```text
Read AGENTS.md and docs/execplans/OMA_BER_MVP_EXECPLAN.md. Implement Milestone 9 and Milestone 10: complete README.md and review the implementation. Check for unit mistakes, dB/dBm mistakes, BER sign mistakes, hidden file I/O, and incomplete examples. Run pytest and ruff. Summarize changed files, commands run, pass/fail status, and remaining limitations.
```

## Prompt 8: Ask Codex to review uncommitted changes

```text
Review the current uncommitted changes against AGENTS.md and docs/execplans/OMA_BER_MVP_EXECPLAN.md. Focus on numerical correctness, physical units, BER sign conventions, input validation, tests, and examples. Do not modify files unless you find a clear bug; if you modify files, run pytest and ruff again.
```

## Prompt 9: Start Phase 2 after MVP is stable

```text
The NRZ/OOK MVP is complete and tested. Read AGENTS.md and propose a new ExecPlan for Phase 2: PD parameter sweeps for responsivity, dark current, TIA noise, and required OMA at target BER. Do not implement yet. Create docs/execplans/PD_PARAMETER_SWEEP_EXECPLAN.md with milestones, APIs, tests, and examples.
```
