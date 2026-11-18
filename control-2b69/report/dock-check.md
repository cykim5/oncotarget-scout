# dock-check: WEAK

**the weakest known binder is at the 88.4th percentile, short of the 95th percentile bar.**

> **Not size matched:** these known binders fall outside the decoy size window, so their percentiles compare molecules of different sizes and should not be read as a result: UDP (25 heavy atoms), UMP (21 heavy atoms).

> The verdict is set by the size matched binders only: UDP-glucuronate, UDP-xylose, UDP-galactose.

## Known binders against size matched decoys

| Ligand | Heavy atoms | Score | Decoys beating it | Percentile |
| --- | --- | --- | --- | --- |
| UMP | 21 | -8.927 | 93 of 199 | 53.3 |
| UDP | 25 | -10.021 | 30 of 199 | 84.9 |
| UDP-glucuronate | 37 | -10.213 | 23 of 199 | 88.4 |
| UDP-galactose | 36 | -10.996 | 2 of 199 | 99.0 |
| UDP-xylose | 34 | -11.126 | 1 of 199 | 99.5 |

## How this was run

- scoring function: `vina`
- exhaustiveness 8, seed 42
- box centre [50.48, 21.14, 23.36], size [28.0, 28.0, 28.0]
- 199 decoys scored, median -8.784
- 36.7 minutes

- score against heavy atom count: rho 0.11 (p=1.2e-01, n=199)

### 1 decoys could not be docked

Listed so the count is not silently smaller than intended.

- `CHEMBL5315122`: Atom type B is not a valid AutoDock type (atom types are cas

## What this does and does not tell you

Passing means the protocol can find a binder it was not told
about, on this target, in this box. It does not mean any
particular compound in a screen binds, and it does not
transfer to a different box or a different target.

Failing is the stronger result. If a known binder cannot beat
random approved drugs here, then nothing in a ranked screen
on this target has been demonstrated, however good the
scores look.

