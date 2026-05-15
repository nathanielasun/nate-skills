# Build, test, and PR workflow

Combined operational reference: building WarpX from source, configuring tests, opening PRs. Read this for setup, debugging build failures, or before submitting code.

## Quick build

```bash
git clone --depth 1 https://github.com/BLAST-WarpX/warpx.git
cd warpx

# 3D Cartesian build (most common):
cmake -S . -B build -DWarpX_DIMS="3"
cmake --build build -j 8

# Resulting binary: build/bin/warpx.3d
```

For multi-dimension build:
```bash
cmake -S . -B build -DWarpX_DIMS="1;2;3;RZ"
cmake --build build -j 8
# Binaries: build/bin/warpx.{1d,2d,3d,rz}
```

## Common CMake options

| Option | Default | Meaning |
|---|---|---|
| `WarpX_DIMS` | "3" | Dimensions to build; "1;2;3;RZ" for all |
| `WarpX_COMPUTE` | OMP | Backend: `OMP`, `CUDA`, `HIP`, `SYCL`, `NOACC` |
| `WarpX_MPI` | ON | MPI support |
| `WarpX_PYTHON` | OFF | Build Python bindings (pywarpx) |
| `WarpX_OPENPMD` | ON | openPMD I/O |
| `WarpX_QED` | ON | QED module |
| `WarpX_FFT` | OFF | FFT (needed for PSATD) |
| `WarpX_PRECISION` | DOUBLE | Field precision: `DOUBLE` or `SINGLE` |
| `WarpX_PARTICLE_PRECISION` | DOUBLE | Particle precision: same options |
| `AMReX_CUDA_ARCH` | (autodetect) | CUDA compute capability, e.g., `80`, `90` |

For typical low-temp plasma work (no PSATD, no QED, with Python):

```bash
cmake -S . -B build \
    -DWarpX_PYTHON=ON \
    -DWarpX_DIMS="1;2;3;RZ" \
    -DWarpX_FFT=OFF \
    -DWarpX_QED=OFF
cmake --build build -j 8
```

For GPU (NVIDIA, A100):
```bash
cmake -S . -B build \
    -DWarpX_COMPUTE=CUDA \
    -DAMReX_CUDA_ARCH=80 \
    -DWarpX_DIMS="3"
cmake --build build -j 8
```

## Faster iteration during development

```bash
# Initial configure once:
cmake -S . -B build -DWarpX_DIMS="3" -DWarpX_PYTHON=ON

# Rebuild after C++ changes:
cmake --build build -j 8

# Or for the Python build specifically:
cmake --build build -j 8 --target pip_install_nodeps
```

`pip_install_nodeps` skips dependency resolution — much faster on rebuild iterations.

For very fast iteration on a single file: set `-DWarpX_PYTHON_IPO=OFF` to skip link-time optimization (slow LTO step). Acceptable for development; turn back on for production.

## Dependencies

Required:
- CMake ≥ 3.24
- C++17 compiler (GCC 9+, Clang 12+, Intel 2021+)
- MPI library (OpenMPI, MPICH, Intel MPI)

Optional, picked up if available:
- HDF5 (for openPMD HDF5 backend)
- ADIOS2 (for openPMD BP backend)
- FFTW or cuFFT (for PSATD)
- CUDA / HIP / SYCL (for GPU)
- Python 3.8+ + pybind11 + mpi4py (for pywarpx)

The Superbuild option (default) downloads dependencies AMReX, picsar, openPMD-api, etc. automatically. To use system-installed versions instead:

```bash
cmake -S . -B build -DWarpX_amrex_internal=OFF -DAMReX_DIR=/path/to/amrex-install/lib/cmake/AMReX
```

Same pattern for `picsar`, `openpmd-api`, etc.

## Running a built binary

```bash
# Serial:
./build/bin/warpx.3d inputs

# MPI:
mpiexec -n 4 ./build/bin/warpx.3d inputs

# With OpenMP threading:
export OMP_NUM_THREADS=4
mpiexec -n 8 ./build/bin/warpx.3d inputs    # 32 cores total
```

The `inputs` file is plain text, key=value format. For PICMI workflows, the script generates this file via `sim.write_input_file('inputs')`.

## Running tests

WarpX uses CTest for the test suite. After build:

```bash
cd build
ctest -j 4              # run all tests
ctest -j 4 -R MCC       # only tests matching regex
ctest --output-on-failure  # verbose for failing tests
```

Tests are defined in `Examples/Tests/<test_name>/` with a `CMakeLists.txt` registering them. Each test has:
- Input file (`inputs_test_*` or `inputs_*_picmi.py`)
- Analysis script (`analysis*.py`) — Python script that checks results
- Reference data in `Regression/Checksum/benchmarks_json/<test_name>.json`

To run a single test manually:

```bash
cd build/bin/test_<name>_run_dir/
./warpx.3d inputs_test_3d_...
python analysis.py
```

## Adding a new test

For a new physics module, you need a test. Pattern:

1. **Create test directory**: `Examples/Tests/<my_feature>/`
2. **Add input file**: `inputs_test_3d_my_feature` (or `_picmi.py` for PICMI)
3. **Add analysis script**: `analysis_default.py` (uses checksum comparison) or `analysis.py` (custom)
4. **Generate reference data**:
   ```bash
   # Run once, save the checksum to benchmarks:
   ./build/bin/warpx.3d inputs_test_3d_my_feature
   python Tools/Checksum/generate_checksum.py --evaluated-file diags/diag1 --output-file Regression/Checksum/benchmarks_json/test_3d_my_feature.json
   ```
5. **Register in CMakeLists.txt**:
   ```cmake
   add_warpx_test(
       test_3d_my_feature
       3
       inputs_test_3d_my_feature
       ""                            # analysis args
       analysis_default.py
       OFF                           # CI on/off (start OFF, turn ON after a green run)
   )
   ```

Keep tests fast — under 1 minute on a modern CPU. Use small grids (16-32 cells), few particles per cell, short runtimes.

## Debugging build failures

| Symptom | Likely cause |
|---|---|
| `undefined reference to ...` (link error) | Missing `.cpp` in `CMakeLists.txt` or `Make.package` |
| `error: 'amrex::Real' was not declared` | Missing AMReX include |
| `error: identifier "..." is undefined` (NVCC) | Function called from device not marked `AMREX_GPU_HOST_DEVICE` |
| `static assertion failed: ...` | Type mismatch, usually Real vs ParticleReal |
| `Cannot create directory ...` during pip_install | Permissions issue; check Python venv is writable |
| Slow build | Try `-DWarpX_PYTHON_IPO=OFF` or fewer `WarpX_DIMS` |

If the error is in AMReX itself, you likely have a version mismatch. Wipe `build/` and reconfigure.

## Linting and static analysis

```bash
# clang-tidy (runs in CI):
cmake -S . -B build -DWarpX_DIMS="3" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build -j 8 -t clang-tidy

# clang-format check:
./Tools/Linter/runClangTidy.sh
```

The pre-commit hook system formats and lints automatically. Install it:

```bash
pre-commit install
```

After this, `git commit` will run formatting checks. To run manually:

```bash
pre-commit run --all-files
```

## Profiling

WarpX has built-in timers:

```bash
# Enable profiling in input file:
amrex.the_arena_is_managed = 0
warpx.do_dynamic_scheduling = 1
warpx.verbose = 2
```

Output goes to `tiny_profile.txt`. For more detail use AMReX's profile output or compile with `-DAMReX_TINY_PROFILE=ON`.

For GPU profiling, use `nsys` (NVIDIA Nsight Systems) or `rocprof` (AMD ROCm).

## PR workflow

### Before opening a PR

Checklist:
- [ ] Code builds in a clean checkout with default options
- [ ] `clang-tidy` and `clang-format` clean (run pre-commit)
- [ ] New code follows the coding rules (`references/coding-rules.md`)
- [ ] At least one test added or existing test extended
- [ ] CI passes (or expected failures are documented)
- [ ] Documentation updated (parameters, theory, implementation as needed)
- [ ] Commits squashed / rebased to a clean history
- [ ] PR description: what + why + how, with link to discussion if exists

### PR scope

- **One feature per PR.** Splitting style changes from substance is critical.
- **Refactors separate from new features.** Reviewers will ask to split anyway.
- **Don't change unrelated formatting.** Adds review burden.

### PR description template (informal)

```markdown
## What
Adds a new MCC process for X-ray ionization of Argon.

## Why
Closes #1234. Required for accurate modeling of laser-induced breakdown.

## How
- Added `MCCProcessType::XrayIonization` enum value
- New `XrayIonizationFunc` implementing the kinematics
- New test in `Examples/Tests/xray_ionization/`
- Cross-section data lives in `BLAST-WarpX/warpx-data` PR #56

## Validation
Compared total ionization rate vs. analytical curve in
`Examples/Tests/xray_ionization/analysis.py`. Agreement within 5% over
the 10-100 keV range.
```

### Review cycle

- CI must be green before merge.
- Two maintainer approvals typically required.
- Expect requests for splitting, tests, documentation. Address them in additional commits; squash on merge.

## Where to ask for help

- **Bug reports** with reproducer: GitHub issues at `https://github.com/BLAST-WarpX/warpx/issues`
- **Usage / design questions**: GitHub Discussions at `https://github.com/BLAST-WarpX/warpx/discussions`
- **Real-time chat**: WarpX Slack (invitation linked from readthedocs landing page)

Include: WarpX version (`git rev-parse HEAD`), OS, compiler, CMake configuration, complete error output. The smaller and more reproducible your example, the faster the response.

## Cheat sheet: common commands

```bash
# Clean build from scratch:
rm -rf build && cmake -S . -B build -DWarpX_DIMS="3" && cmake --build build -j

# Run a small test:
mpiexec -n 2 ./build/bin/warpx.3d Examples/Tests/<test>/inputs_test_3d_*

# Run with FPE traps (helpful for debugging NaN):
./build/bin/warpx.3d inputs amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1

# Higher verbosity:
./build/bin/warpx.3d inputs warpx.verbose=2

# Check available WarpX C++ APIs (Doxygen browse):
xdg-open https://warpx.readthedocs.io/en/latest/_static/doxyhtml/index.html
```
