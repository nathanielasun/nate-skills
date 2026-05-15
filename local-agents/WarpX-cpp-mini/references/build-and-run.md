# Build, run, test, PR

## Quick build

```bash
git clone --depth 1 https://github.com/BLAST-WarpX/warpx.git
cd warpx

# 3D Cartesian (most common):
cmake -S . -B build -DWarpX_DIMS="3"
cmake --build build -j 8

# Resulting binary: build/bin/warpx.3d
```

All dimensions:
```bash
cmake -S . -B build -DWarpX_DIMS="1;2;3;RZ"
cmake --build build -j 8
```

## Common cmake options

| Option | Default | Meaning |
|---|---|---|
| `WarpX_DIMS` | `"3"` | Dimensions; `"1;2;3;RZ"` for all |
| `WarpX_COMPUTE` | `OMP` | Backend: `OMP`, `CUDA`, `HIP`, `SYCL`, `NOACC` |
| `WarpX_MPI` | `ON` | MPI support |
| `WarpX_PYTHON` | `OFF` | Build pywarpx Python bindings |
| `WarpX_FFT` | `OFF` | FFT (for PSATD) |
| `WarpX_QED` | `ON` | QED module |
| `WarpX_PRECISION` | `DOUBLE` | Field precision: `DOUBLE` or `SINGLE` |
| `WarpX_PARTICLE_PRECISION` | `DOUBLE` | Particle precision |
| `AMReX_CUDA_ARCH` | autodetect | CUDA capability: `80` (A100), `90` (H100) |

For typical low-temp plasma development (no PSATD, no QED, with Python):

```bash
cmake -S . -B build \
    -DWarpX_PYTHON=ON \
    -DWarpX_DIMS="1;2;3;RZ" \
    -DWarpX_FFT=OFF \
    -DWarpX_QED=OFF
cmake --build build -j 8
```

For GPU (A100):
```bash
cmake -S . -B build \
    -DWarpX_COMPUTE=CUDA \
    -DAMReX_CUDA_ARCH=80 \
    -DWarpX_DIMS="3"
cmake --build build -j 8
```

## Faster rebuild

```bash
# After C++ changes:
cmake --build build -j 8

# For Python builds:
cmake --build build -j 8 --target pip_install_nodeps
```

`pip_install_nodeps` skips dependency resolution — fastest.

## Run

```bash
# Serial:
./build/bin/warpx.3d inputs

# MPI:
mpiexec -n 4 ./build/bin/warpx.3d inputs

# Hybrid MPI + OpenMP:
export OMP_NUM_THREADS=4
mpiexec -n 8 ./build/bin/warpx.3d inputs        # 32 cores total
```

## Run tests

```bash
cd build
ctest -j 4                                       # all tests
ctest -j 4 -R MCC                                 # tests matching regex
ctest --output-on-failure                         # verbose on failure
```

Single test manually:
```bash
cd build/bin/test_<name>_run_dir/
./warpx.3d inputs_test_3d_...
python analysis.py
```

## Add a new test

Pattern:

1. Create directory: `Examples/Tests/<my_feature>/`
2. Add input: `inputs_test_3d_my_feature`
3. Add analysis script: `analysis_default.py` (uses checksum) or `analysis.py` (custom)
4. Generate reference data:
   ```bash
   ./build/bin/warpx.3d inputs_test_3d_my_feature
   python Tools/Checksum/generate_checksum.py \
       --evaluated-file diags/diag1 \
       --output-file Regression/Checksum/benchmarks_json/test_3d_my_feature.json
   ```
5. Register in `Examples/Tests/<my_feature>/CMakeLists.txt`:
   ```cmake
   add_warpx_test(
       test_3d_my_feature
       3
       inputs_test_3d_my_feature
       ""
       analysis_default.py
       OFF
   )
   ```

Keep tests fast — under 1 minute on a modern CPU. Use small grids (16-32 cells).

## Linting

```bash
# pre-commit hook (formats + lints):
pre-commit install                                # one-time setup
pre-commit run --all-files                        # run manually

# clang-tidy (runs in CI):
cmake -S . -B build -DWarpX_DIMS="3" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
./Tools/Linter/runClangTidy.sh
```

## Build failure cheat sheet

| Symptom | Likely cause |
|---|---|
| `undefined reference to ...` | Missing `.cpp` in `CMakeLists.txt` OR `Make.package` |
| `error: 'amrex::Real' was not declared` | Missing AMReX include |
| `error: identifier "..." is undefined` (NVCC) | Function in device code not marked `AMREX_GPU_HOST_DEVICE` |
| `static assertion failed: ...` | Type mismatch — usually Real vs ParticleReal |
| Slow build | Try `-DWarpX_PYTHON_IPO=OFF` or fewer `WarpX_DIMS` |
| AMReX errors after update | Wipe `build/` and reconfigure |

## Runtime debugging

### Trap floating-point exceptions
```bash
./build/bin/warpx.3d inputs \
    amrex.fpe_trap_invalid=1 \
    amrex.fpe_trap_zero=1 \
    amrex.fpe_trap_overflow=1
```

### Verbose output
```bash
./build/bin/warpx.3d inputs warpx.verbose=2
```

### Common runtime errors

| Error | Fix |
|---|---|
| "Particle outside box" | Position was modified without `pc.Redistribute()` |
| "PoissonSolver did not converge" | Tighten `required_precision` or increase `maximum_iterations`; check BCs |
| "ParmParse: parameter not understood" | Misspelled `warpx_*` arg → bad C++ parameter name |
| Segfault in MCC | Check cross-section file paths and product species names |

When in doubt, search the WarpX issue tracker for the error message.

## PR workflow

Before opening a PR:

- [ ] Builds clean with default options
- [ ] `clang-tidy` and `clang-format` pass (run `pre-commit`)
- [ ] New code follows `references/patterns.md` rules
- [ ] At least one test added or existing test extended
- [ ] CI passes
- [ ] Documentation updated (parameters → `Docs/source/usage/parameters.rst`)
- [ ] Commits squashed / clean history
- [ ] PR description: what + why + how

PR rules:
- **One feature per PR.**
- **Refactors separate from new features.**
- **No unrelated formatting changes.**

PR description template:
```
## What
[Brief feature description]

## Why
[Motivation; link to discussion if exists]

## How
- [Changes]
- [New tests]

## Validation
[How verified]
```

## Where to ask

- Bugs: `https://github.com/BLAST-WarpX/warpx/issues`
- Usage / design: `https://github.com/BLAST-WarpX/warpx/discussions`
- WarpX Slack (invite from readthedocs)

Always include: WarpX version (`git rev-parse HEAD`), OS, compiler, cmake config, complete error.
