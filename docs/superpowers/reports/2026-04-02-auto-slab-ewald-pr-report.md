# Auto Slab Ewald Regression Test

## Summary

This change adds a focused unit test for `kspace_modify slab auto` on the `quasi-2d` branch.

The test uses:

- `kspace_style ewald`
- a neutral 100-particle deterministic random `p p f` fixture
- a `40 x 40 x 40` simulation box
- `ewald 1e-8` with `kspace_modify slab 10.0` as the fixed force reference

For target accuracies `1e-3`, `1e-4`, and `1e-5`, the test compares the total force field from `kspace_modify slab auto` against that tight reference and asserts:

`L2(force error) < 10 * accuracy`

## Motivation

Earlier exploratory checks showed that `slab auto` can differ significantly from `slab 10.0` on quasi-2D systems. The final regression is intentionally simpler and enforces a pragmatic acceptance criterion against a much tighter Ewald reference, instead of encoding exploratory sweeps or solver-specific tuning behavior into the test suite.

## Implementation

The regression is implemented in:

- `unittest/commands/test_kspace_auto_slab.cpp`

The fixture:

- creates 100 particles at deterministic pseudo-random positions
- alternates `+1/-1` charges to keep the system neutral
- uses `pair_style coul/long 8.0`
- runs Ewald with `slab auto` for the tested accuracies
- compares against an `ewald 1e-8` and `slab 10.0` reference

## Validation

Focused verification used:

```bash
module load gcc/13.3.0 && cmake --build build --target test_kspace_auto_slab -j2
module load gcc/13.3.0 && ./build/test_kspace_auto_slab
module load gcc/13.3.0 && ctest --test-dir build --output-on-failure -R KSpaceAutoSlab
```

Observed pass condition:

- `1e-3`: error below `1e-2`
- `1e-4`: error below `1e-3`
- `1e-5`: error below `1e-4`

## PR Notes

This is a test-only change plus this report artifact. No production solver code was modified in this final version.
