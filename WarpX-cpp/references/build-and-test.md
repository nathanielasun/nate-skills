# Build and test

This reference covers the build system, dependencies, local development workflows, testing, profiling, and debugging. Read this when setting up a new development environment, when a build is failing, or when adding tests for new code.

## Table of contents

1. [Dependencies and build prerequisites](#dependencies-and-build-prerequisites)
2. [The CMake-based build, in one page](#the-cmake-based-build-in-one-page)
3. [CMake options that matter day-to-day](#cmake-options-that-matter-day-to-day)
4. [Local development workflows: faster compilation](#local-development-workflows-faster-compilation)
5. [Building the Python interface](#building-the-python-interface)
6. [The legacy GNUmake build (still used)](#the-legacy-gnumake-build-still-used)
7. [Adding new source files](#adding-new-source-files)
8. [Testing: writing and running](#testing-writing-and-running)
9. [clang-tidy: the linter](#clang-tidy-the-linter)
10. [Profiling](#profiling)
11. [Debugging](#debugging)

## Dependencies and build prerequisites

WarpX requires:
- **CMake ≥ 3.24** (the project tracks recent CMake; check the top-level `CMakeLists.txt` for the exact minimum).
- A C++17 compiler. GCC 9+, Clang 11+, Intel oneAPI, AMD ROCm clang, MSVC 2019+ are supported.
- **MPI** (optional but typical): MPICH, OpenMPI, Cray MPI, Intel MPI.
- **CUDA / HIP / SYCL** if building for GPU: CUDA 11.7+ for Nvidia, ROCm 5+ for AMD, oneAPI for Intel.
- **Python 3.10+** if building the PICMI interface.
- **FFTW or cuFFT/rocFFT/oneMKL** if `WarpX_FFT=ON` (PSATD solver).
- **HDF5 + openPMD-api** typically (built automatically as a dependency by default).

The major dependencies WarpX itself builds (the "superbuild" mode) are AMReX, PICSAR (QED), openPMD-api, pyAMReX, and pybind11. You can opt into using pre-installed versions of these via package manager (conda-forge, Spack, Brew, APT) if you need a specific configuration.

### Quick start with conda-forge (recommended for laptops)

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict

conda create -n warpx-dev -c conda-forge \
    blaspp boost ccache cmake compilers git lapackpp \
    "openpmd-api=*=mpi_mpich*" openpmd-viewer packaging \
    pytest python python-build make numpy pandas scipy \
    setuptools yt "fftw=*=mpi_mpich*" pkg-config matplotlib \
    mamba mpich mpi4py ninja pip virtualenv wheel

conda activate warpx-dev
# Then build with -DWarpX_MPI=ON
```

For Linux with NVIDIA GPU: add `mamba install -c nvidia -c conda-forge cuda cuda-nvtx-dev cupy`.

For macOS or Windows OpenMP support: `mamba install -c conda-forge llvm-openmp`.

### HPC systems

WarpX maintains profile scripts for major DOE/NSF systems (Frontier, Perlmutter, Polaris, Summit, Crusher, Adastra, Lassen, Cori, etc.) under `Tools/machines/<system-name>/`. Each contains a `<system>_warpx.profile.example` to source for environment setup. See the HPC Systems page in the readthedocs documentation for system-specific instructions.

## The CMake-based build, in one page

From the WarpX source root:

```bash
# Configure (creates build/ directory with cached options)
cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_MPI=ON -DWarpX_COMPUTE=OMP

# Compile (parallel build with N threads)
cmake --build build -j 8

# The binary appears in build/bin/warpx.3d.MPI.OMP.DP.OPMD.QED (the suffix encodes the configuration)
```

For multiple dimensions in one build:

```bash
cmake -S . -B build -DWarpX_DIMS="1;2;3;RZ"
```

This produces separate binaries: `warpx.1d`, `warpx.2d`, `warpx.3d`, `warpx.RZ`. Each is a monolithic binary for that geometry.

To rebuild after source changes, just rerun the build step:

```bash
cmake --build build -j 8
```

To reconfigure (after changing options or upgrading dependencies):

```bash
rm -rf build
cmake -S . -B build -DWarpX_DIMS=3 [...other options]
cmake --build build -j 8
```

You can inspect and modify options with `ccmake build` (interactive) or by passing `-D<OPTION>=<VALUE>` flags.

### Build types

The CMake `CMAKE_BUILD_TYPE` controls optimization and debug symbols:
- **`Release`** (default): full optimization, no debug symbols. Use for production runs.
- **`RelWithDebInfo`**: optimization with debug symbols. Use for profiling.
- **`Debug`**: no optimization, full debug symbols, runtime assertions enabled. Use when diagnosing wrong-output bugs. Very slow.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
```

For *minimal* debug symbols (just enough for stack traces, no performance cost), use `-DWarpX_BACKTRACE_INFO=ON`.

## CMake options that matter day-to-day

These are the ones you'll touch most often. See `https://warpx.readthedocs.io/en/latest/install/cmake.html` for the complete list.

| Option | Values (default in **bold**) | What it does |
|---|---|---|
| `WarpX_DIMS` | **`3`** / `1` / `2` / `RZ` / `RCYLINDER` / `RSPHERE`; semicolon-separated for multiple | Which geometries to build binaries for |
| `WarpX_COMPUTE` | `NOACC` / **`OMP`** / `CUDA` / `SYCL` / `HIP` | Compute backend |
| `WarpX_MPI` | **`ON`** / `OFF` | MPI support |
| `WarpX_PRECISION` | `SINGLE` / **`DOUBLE`** | Field floating-point precision |
| `WarpX_PARTICLE_PRECISION` | `SINGLE` / **`DOUBLE`** (default = `WarpX_PRECISION`) | Particle precision (can differ from field precision) |
| `WarpX_EB` | **`ON`** / `OFF` | Embedded boundaries (not supported in RZ/RCYLINDER/RSPHERE) |
| `WarpX_FFT` | `ON` / **`OFF`** | FFT-based solvers (PSATD, Ohm's law magnetostatic, IGF for ES) |
| `WarpX_OPENPMD` | **`ON`** / `OFF` | openPMD I/O (HDF5/ADIOS) |
| `WarpX_PYTHON` | `ON` / **`OFF`** | Python (PICMI) bindings — needed for input scripts |
| `WarpX_QED` | **`ON`** / `OFF` | QED physics via PICSAR |
| `WarpX_PETSC` | `ON` / **`OFF`** | PETSc nonlinear/linear solvers (advanced; for implicit EM) |
| `WarpX_APP` | **`ON`** / `OFF` | Build the standalone executable |
| `WarpX_LIB` | `ON` / **`OFF`** | Build as a library (auto-enabled by `WarpX_PYTHON`) |
| `WarpX_IPO` | `ON` / **`OFF`** | Link-time optimization (slower compile, faster runtime) |
| `CMAKE_BUILD_TYPE` | `Release` / `RelWithDebInfo` / `Debug` | Optimization level |

For development cycles:

```bash
# Fast iteration: debug symbols but optimized, no LTO
cmake -S . -B build \
    -DWarpX_DIMS=2 \
    -DWarpX_COMPUTE=OMP \
    -DWarpX_MPI=OFF \
    -DWarpX_QED=OFF \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DWarpX_BACKTRACE_INFO=ON
```

Building only 2D and turning off MPI/QED dramatically reduces compile time during development.

## Local development workflows: faster compilation

WarpX compile times are non-trivial (5-30 minutes for full rebuilds). Several techniques help.

### CCache

Enabled by default. Install `ccache` and CMake will find it automatically. Subsequent recompiles of unchanged translation units use the cache (nearly instant). To verify it's working:

```bash
ccache -s   # show stats
```

For developers who switch between many configurations, increase the cache size (default ~5 GB):

```bash
ccache --set-config max_size=50G
```

### Ninja instead of Make

```bash
cmake -S . -B build -G Ninja
cmake --build build
```

Ninja is faster than Make for parallel builds and produces better build-progress output.

### Working against local AMReX/PICSAR/openPMD source

If you're modifying AMReX or another dependency in tandem with WarpX:

```bash
cmake -S . -B build \
    -DWarpX_amrex_src=$HOME/src/amrex \
    -DWarpX_openpmd_src=$HOME/src/openPMD-api \
    -DWarpX_picsar_src=$HOME/src/picsar \
    -DWarpX_pyamrex_src=$HOME/src/pyamrex \
    -DWarpX_pybind11_src=$HOME/src/pybind11
```

WarpX will use the local clones rather than downloading its preferred versions. Changes to those local repos take effect on the next WarpX rebuild.

### Pre-installed dependencies

If you've installed AMReX/PICSAR/etc. via conda or system packages and want WarpX to use them:

```bash
cmake -S . -B build \
    -DWarpX_amrex_internal=OFF \
    -DWarpX_openpmd_internal=OFF \
    -DWarpX_picsar_internal=OFF \
    -DWarpX_pyamrex_internal=OFF \
    -DWarpX_pybind11_internal=OFF
```

Make sure `CMAKE_PREFIX_PATH` includes the installation prefixes of those libraries.

### Switching dimensions: a separate build directory

Each `WarpX_DIMS` configuration is a separate set of binaries. Keep separate build directories to avoid full rebuilds when switching:

```bash
cmake -S . -B build_3d -DWarpX_DIMS=3 -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake -S . -B build_2d -DWarpX_DIMS=2 -DCMAKE_BUILD_TYPE=RelWithDebInfo
# ... etc
```

## Building the Python interface

The PICMI Python interface (`pywarpx`) requires `WarpX_PYTHON=ON`:

```bash
cmake -S . -B build_py -DWarpX_DIMS="1;2;3;RZ" -DWarpX_PYTHON=ON
cmake --build build_py --target pip_install -j 8
```

This installs `pywarpx` into the active Python environment.

For repeated development cycles, skip dependency checks for speed:

```bash
cmake --build build_py --target pip_install_nodeps -j 8
```

For *much* faster Python builds, skip link-time optimization (loses some runtime performance):

```bash
cmake -S . -B build_py \
    -DWarpX_PYTHON=ON \
    -DWarpX_PYTHON_IPO=OFF \
    -DpyAMReX_IPO=OFF
```

Note: the Python interface is the domain of the (anticipated) `warpx-python` skill — this section is here only because the C++ side may need to build the Python bindings to test some changes (e.g., a new C++ API exposed to Python).

## The legacy GNUmake build (still used)

WarpX also maintains a GNUmake build system (inherited from AMReX). It's used in CI for some configurations and by some HPC users. Files: `Source/Make.package` and per-subdirectory `Make.package`. The build is invoked via:

```bash
cd Tools/machines/.../  # find a working makefile setup
make -j 8 DIM=3 COMP=gcc USE_MPI=TRUE USE_OMP=TRUE
```

**When you add a new `.cpp` file, you must add it to both `CMakeLists.txt` AND `Make.package` in the same subdirectory.** Headers (`.H` files) are not listed in either — they're picked up by the compiler include search.

The two build systems should produce identical binaries for the same options. If they diverge, that's a bug; either system can be the source of truth depending on how the file lists are maintained, but in practice CMake is more actively developed.

## Adding new source files

To add `Source/Particles/Collision/MyNewProcess.cpp` and `Source/Particles/Collision/MyNewProcess.H`:

1. Create the files. Headers are not built directly, but follow the include conventions in `references/style-and-conventions.md`.

2. Edit `Source/Particles/Collision/CMakeLists.txt`:
   ```cmake
   target_sources(WarpX
       PRIVATE
           CollisionBase.cpp
           CollisionHandler.cpp
           MyNewProcess.cpp   # <-- add this line
   )
   ```

3. Edit `Source/Particles/Collision/Make.package`:
   ```makefile
   CEXE_sources += CollisionBase.cpp
   CEXE_sources += CollisionHandler.cpp
   CEXE_sources += MyNewProcess.cpp   # <-- add this line
   CEXE_headers += CollisionBase.H
   CEXE_headers += MyNewProcess.H     # <-- and this (headers ARE listed in Make.package, not CMakeLists)
   ```

Note: the convention is that **headers are listed in `Make.package`** (for the legacy build) **but not in `CMakeLists.txt`** (where CMake finds them transitively). The contributing docs explicitly call this out: "Do not list header files (`.H`) here" referring to CMakeLists. Headers in `Make.package` are listed under `CEXE_headers`.

4. Rebuild. CMake re-runs configuration automatically when `CMakeLists.txt` changes.

## Testing: writing and running

WarpX has an extensive automated test suite (over 200 tests as of late 2025). New features without tests will not be accepted.

### Test types

- **Regression tests** (the bulk): a full WarpX simulation with a known input file. The simulation produces a checksum (a hash of key output quantities) that is compared against a previously-recorded baseline.
- **Unit tests**: smaller, faster, exercise a specific C++ function or class.
- **Analysis scripts**: optional Python script that runs after the simulation and computes additional diagnostics or comparisons (e.g., compare to an analytic solution).

### Running tests locally

Tests are built as CTest targets when `WarpX_TEST_*` options are enabled. To enable test cleanup, FPE-trapping, etc.:

```bash
cmake -S . -B build \
    -DWarpX_DIMS=3 \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DWarpX_TEST_CLEANUP=ON \
    -DWarpX_TEST_FPETRAP=ON \
    -DWarpX_BACKTRACE_INFO=ON \
    -DBUILD_TESTING=ON

cmake --build build -j 8
ctest --test-dir build -j 4 --output-on-failure
```

To run a single test by name pattern:

```bash
ctest --test-dir build -R "name_of_test" --output-on-failure -V
```

### Writing a regression test

The pattern: copy an existing test that's structurally similar, modify the input file, and update the checksum baseline.

1. **Pick an example to mimic.** For a new MCC process, find an existing MCC test under `Tests/Tests/`. For a new push, find an existing push test.

2. **Create a test directory**: `Tests/Tests/test_<descriptive_name>/`. Contents typically include:
   - `inputs_test_<name>`: the WarpX input file.
   - `analysis_<name>.py`: optional Python analysis (compares to analytic solution, etc.).
   - `<name>.json`: the checksum baseline (will be created when you run the test the first time).

3. **Register the test** in `Tests/Tests/CMakeLists.txt` or a closer-scoped CMakeLists. Use `add_warpx_test()` (a project-provided helper):
   ```cmake
   add_warpx_test(
       test_<name>      # test name
       3                # dimensionality
       1                # number of MPI ranks
       inputs_test_<name>   # input file
       analysis_<name>.py   # analysis script (or OFF)
       diags/diag1/        # output directory pattern for checksum
       OFF                 # restart test
   )
   ```

4. **Run the test, capture the baseline**: the first run will fail because no baseline exists. Use the helper script (typically `Tools/DevUtils/checksum/reset_all_checksums.sh` or a per-test variant) to write the current output as the new baseline. Commit the resulting JSON.

5. **Iterate**: run `ctest` to verify the test passes deterministically on multiple runs.

### Numerical reproducibility

Tests need to be reproducible. Sources of nondeterminism to be aware of:
- **GPU vs CPU**: bit-exact reproducibility across hardware is not guaranteed. Tests typically have separate baselines or use looser tolerances.
- **MPI rank count**: with collision modules using random numbers, different rank counts produce different results. Tests usually fix the rank count.
- **OpenMP threading**: similar concerns. CI runs with controlled thread counts.

For statistical tests (e.g., verifying an ionization rate matches an analytic value), use enough macroparticles that the Monte Carlo noise is well below the test tolerance, and choose tolerances appropriately.

## clang-tidy: the linter

WarpX uses `clang-tidy` for static analysis. CI runs it on every PR; local runs catch issues before submission.

```bash
# Generate compile_commands.json (needed by clang-tidy)
cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

# Run clang-tidy on changed files (see Tools/Linter/runClangTidy.sh or similar)
./Tools/Linter/runClangTidy.sh
```

Configuration is in `.clang-tidy` at the repo root. Common warnings flagged:
- Missing `m_` prefix on a member variable.
- `using namespace` at file scope.
- Missing `const` on a parameter that should be const.
- Range-based for loop on an integer type (suggests using `iota_view` or similar).
- Implicit conversions between signed/unsigned types.

For details on running clang-tidy, see the `How to run the clang-tidy linter` page in the developer docs.

## Profiling

WarpX is heavily profiled. AMReX provides `BL_PROFILE` macros (used in WarpX as `WARPX_PROFILE` historically, now `ABLASTR_PROFILE`):

```cpp
void my_function ()
{
    ABLASTR_PROFILE("my_function");
    // ... function body ...
}
```

These produce tiny structured timing data at runtime. After a run, the AMReX tiny profiler outputs a per-function timing summary.

For deeper profiling:
- **NVIDIA Nsight Systems / Nsight Compute** for CUDA builds.
- **rocprof** for HIP builds.
- **VTune** for CPU builds.
- **AMReX's full profiler** (more overhead, more detail).

See the `How to profile the code` page in the developer documentation for current tooling guidance.

## Debugging

### Compiler symbols and assertions

Build with `-DCMAKE_BUILD_TYPE=Debug` to enable runtime assertions and unoptimized code. The build is much slower but bugs are easier to catch.

A middle ground: `-DCMAKE_BUILD_TYPE=RelWithDebInfo -DWarpX_BACKTRACE_INFO=ON` gives optimized code with enough debug info for stack traces.

### FPE trapping

Floating-point exceptions (division by zero, invalid operations like `0/0` or `sqrt(-1)`, overflow) are silent by default in C++. WarpX can trap them via input file:

```
amrex.fpe_trap_invalid=1
amrex.fpe_trap_zero=1
amrex.fpe_trap_overflow=1
```

When a trap fires, the simulation aborts with a stack trace pointing to the offending operation. Indispensable for diagnosing NaN propagation.

CMake option `-DWarpX_TEST_FPETRAP=ON` enables these traps for all CI tests.

### Stack traces

AMReX installs a SIGSEGV/SIGBUS handler that prints a stack trace on crash. With `-DWarpX_BACKTRACE_INFO=ON`, the trace includes function names and line numbers.

To disable AMReX's signal handlers (so you can attach a debugger):

```
amrex.signal_handling=0
```

Or set `-DWarpX_TEST_DEBUGGER=ON` at build time.

### gdb / lldb / cuda-gdb

For interactive debugging on CPU:

```bash
gdb --args ./build/bin/warpx.3d inputs_test_problem
# Then in gdb: run, breakpoint, etc.
```

For CUDA:

```bash
cuda-gdb --args ./build/bin/warpx.3d.MPI.CUDA inputs_test_problem
```

For MPI debugging, start one rank under a debugger:

```bash
mpirun -np 4 xterm -e gdb --args ./build/bin/warpx.3d inputs_test_problem
```

(Requires X forwarding or a terminal multiplexer setup.)

### Common runtime crashes

- **Out of memory**: increase `amrex.fpe_trap_*=0`, reduce particle count, decrease grid resolution, or use a larger node. For GPU runs, check `nvidia-smi`/`rocm-smi` for memory pressure.
- **Particle goes off-grid**: usually a particle pusher producing a wildly large velocity. Trace back to see what created the bad state (boundary condition? Improperly initialized field?). FPE trapping often catches this earlier than the symptoms appear.
- **MPI hang**: usually a missing `FillBoundary` or `SumBoundary` on one rank but not others — the ranks deadlock waiting for each other. Reproduce with a small case under a debugger.
- **Wrong result on GPU, right on CPU**: almost always implicit `this` capture in a `ParallelFor` lambda. Check for any unqualified member access in kernels.
