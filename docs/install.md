# Installation

cellAdmix-core is built from a shared C++ core plus language bindings. The C++
dependencies are required for both the R package and the Python package.

## System Dependencies

The core build uses CMake and a C++20 compiler. Install:

- CMake 3.16 or newer for direct CMake builds; CMake 3.20 or newer to use
  `CMakePresets.json`.
- A C++20 compiler toolchain.
- `pkg-config`.
- Apache Arrow C++ and Parquet C++.
- Eigen3.
- nlohmann JSON.
- zlib.
- libTIFF.
- OpenJPEG.

On Ubuntu-like systems the package names are typically:

```bash
sudo apt-get install build-essential cmake pkg-config \
  libarrow-dev libparquet-dev libeigen3-dev nlohmann-json3-dev zlib1g-dev \
  libtiff-dev libopenjp2-7-dev
```

Arrow package names and availability vary by distribution. On GitHub Actions
and newer Ubuntu systems, the CI installs Arrow from the Apache Arrow apt
repository before installing `libarrow-dev` and `libparquet-dev`.

The CMake build first looks for imported CMake targets, then falls back to
`pkg-config` where needed. If dependencies are installed in a custom prefix,
set `CMAKE_PREFIX_PATH` for CMake builds or `PKG_CONFIG_PATH` for R builds.

Clone the repository:

```bash
git clone https://github.com/kharchenkolab/cellAdmix-core.git
cd cellAdmix-core
```

## R Bindings

Install the R package from the repository checkout:

```bash
MAKEFLAGS=-j8 R CMD INSTALL r
```

The `MAKEFLAGS` setting is optional, but it speeds up local compilation.

Optional R packages used by examples and plotting include:

```r
install.packages(c("ggplot2", "cowplot", "ggrepel", "rmarkdown", "knitr"))
```

Seurat integration additionally requires `Seurat` and `SeuratObject`.

Verify the R installation:

```bash
Rscript -e 'library(cellAdmixCore); packageVersion("cellAdmixCore")'
```

Because the R package compiles C++ sources from the repository root, the
supported R install route is cloning the repository and running
`R CMD INSTALL r` from the checkout.

## Python Bindings

The Python bindings are built with `scikit-build-core` and `pybind11`. The core
package currently supports Python 3.8 or newer. SpatialData integration should
use Python 3.11 or newer because current SpatialData packages require it.

Create and activate an environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
```

Install the core Python package:

```bash
python -m pip install -e python
```

The default (isolated) build fetches the Python build requirements
(`scikit-build-core`, `pybind11`, and CMake/Ninja wheels) automatically; only
the C++ system dependencies above need to be present. For repeated
development rebuilds, `--no-build-isolation` is faster but skips that
automatic step, so install the build tools into the environment first:

```bash
python -m pip install scikit-build-core pybind11 cmake ninja
python -m pip install -e python --no-build-isolation
```

Install plotting and SpatialData extras when needed:

```bash
python -m pip install -e "python[plot,spatialdata]"
```

If C++ dependencies are installed in a non-standard location, pass the prefix to
CMake through the build environment:

```bash
CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/path/to/prefix" \
  python -m pip install -e python --no-build-isolation
```

Verify the Python installation:

```bash
python -c "import celladmix as ca; print(ca.__version__)"
```

Direct GitHub installation for Python is expected to work because the Python
subdirectory build can still access the repository root during the
`scikit-build-core` build:

```bash
python -m pip install \
  "git+https://github.com/kharchenkolab/cellAdmix-core.git#subdirectory=python"
```

With optional plotting and SpatialData dependencies:

```bash
python -m pip install \
  "celladmix[plot,spatialdata] @ git+https://github.com/kharchenkolab/cellAdmix-core.git#subdirectory=python"
```

## Build Notes

The repository-level CMake build can be used for C++ tests and development:

```bash
cmake --preset tests
cmake --build --preset tests -j8
ctest --preset tests
```

For a no-tests optimized library build:

```bash
cmake --preset user
cmake --build --preset user -j8
```

If you prefer user-space dependency installation through vcpkg, install vcpkg,
set `VCPKG_ROOT`, then use:

```bash
cmake --preset user-vcpkg
cmake --build --preset user-vcpkg -j8
```

The GitHub Actions workflow builds the C++ core on Ubuntu, macOS, and
Windows, and runs Ubuntu smoke jobs for C++ tests, R install/load, and Python
install/import.
