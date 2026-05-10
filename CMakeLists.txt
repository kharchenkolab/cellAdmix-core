cmake_minimum_required(VERSION 3.16)

project(celladmix_core LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CELLADMIX_BUILD_TESTS "Build C++ tests" ON)
option(CELLADMIX_BUILD_BENCHMARKS "Build local development benchmarks" OFF)
option(CELLADMIX_BUILD_PYTHON "Build Python extension module" OFF)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install" CACHE PATH "Default cellAdmix install prefix." FORCE)
endif()

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(APPLE)
  foreach(_homebrew_prefix /opt/homebrew /usr/local)
    if(EXISTS "${_homebrew_prefix}")
      list(APPEND CMAKE_PREFIX_PATH "${_homebrew_prefix}")
    endif()
  endforeach()
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CellAdmixDependencies)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

celladmix_find_dependency(nlohmann_json
  TARGETS nlohmann_json::nlohmann_json
  PKG_CONFIG_NAME nlohmann_json
  OUT_TARGET CELLADMIX_NLOHMANN_JSON_TARGET
  INSTALL_HINT "Ubuntu: sudo apt-get install nlohmann-json3-dev\nmacOS: brew install nlohmann-json\nvcpkg: vcpkg install nlohmann-json")
celladmix_find_dependency(Eigen3
  TARGETS Eigen3::Eigen
  PKG_CONFIG_NAME eigen3
  OUT_TARGET CELLADMIX_EIGEN_TARGET
  INSTALL_HINT "Ubuntu: sudo apt-get install libeigen3-dev\nmacOS: brew install eigen\nvcpkg: vcpkg install eigen3")
celladmix_find_dependency(Parquet
  TARGETS Parquet::parquet_shared Parquet::parquet_static Parquet::parquet
  PKG_CONFIG_NAME parquet
  OUT_TARGET CELLADMIX_PARQUET_TARGET
  INSTALL_HINT "Ubuntu: install libparquet-dev from the Apache Arrow apt repository\nmacOS: brew install apache-arrow\nvcpkg: vcpkg install 'arrow[parquet]'")
celladmix_find_dependency(Arrow
  TARGETS Arrow::arrow_shared Arrow::arrow_static Arrow::arrow
  PKG_CONFIG_NAME arrow
  OUT_TARGET CELLADMIX_ARROW_TARGET
  INSTALL_HINT "Ubuntu: install libarrow-dev from the Apache Arrow apt repository\nmacOS: brew install apache-arrow\nvcpkg: vcpkg install 'arrow[compute,csv,parquet]'")
celladmix_find_dependency(TIFF
  TARGETS TIFF::TIFF
  PKG_CONFIG_NAME libtiff-4
  OUT_TARGET CELLADMIX_TIFF_TARGET
  INSTALL_HINT "Ubuntu: sudo apt-get install libtiff-dev\nmacOS: brew install libtiff\nvcpkg: vcpkg install tiff")
celladmix_find_dependency(OpenJPEG
  TARGETS OpenJPEG::OpenJPEG OpenJPEG::openjp2 openjp2
  PKG_CONFIG_NAME libopenjp2
  PREFER_PKG_CONFIG
  OUT_TARGET CELLADMIX_OPENJPEG_TARGET
  INSTALL_HINT "Ubuntu: sudo apt-get install libopenjp2-7-dev\nmacOS: brew install openjpeg\nvcpkg: vcpkg install openjpeg")

add_library(celladmix_core
  src/clustering.cpp
  src/bridge.cpp
  src/coherence.cpp
  src/domains.cpp
  src/graph.cpp
  src/input_store.cpp
  src/membrane.cpp
  src/ncv.cpp
  src/nmf_euclidean.cpp
  src/nmf_kl.cpp
  src/pipeline.cpp
  src/pipeline_common.cpp
  src/pipeline_store.cpp
  src/report_embedding.cpp
  src/run_store.cpp
  src/simulation.cpp
  src/tabular.cpp
  src/tabular_mask.cpp
  src/workflow.cpp
  src/xenium.cpp
)

if (CELLADMIX_BUILD_PYTHON)
  set_target_properties(celladmix_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

target_include_directories(celladmix_core PUBLIC include)
target_link_libraries(celladmix_core PUBLIC
  ${CELLADMIX_ARROW_TARGET}
  ${CELLADMIX_PARQUET_TARGET}
  ${CELLADMIX_EIGEN_TARGET}
  ${CELLADMIX_NLOHMANN_JSON_TARGET}
  ${CELLADMIX_OPENJPEG_TARGET}
  ${CELLADMIX_TIFF_TARGET}
  ZLIB::ZLIB
  Threads::Threads)

if (CELLADMIX_BUILD_PYTHON)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_core python/src/celladmix/_core_bindings.cpp)
  target_link_libraries(_core PRIVATE celladmix_core)
  install(TARGETS _core LIBRARY DESTINATION celladmix)
endif()

if (CELLADMIX_BUILD_TESTS)
  add_executable(celladmix_core_tests
    tests/test_clustering.cpp
    tests/test_main.cpp
    tests/test_bridge.cpp
    tests/test_graph.cpp
    tests/test_membrane.cpp
    tests/test_ncv.cpp
    tests/test_nmf.cpp
    tests/test_run_store.cpp
    tests/test_simulation.cpp
    tests/test_tabular.cpp
    tests/test_workflow.cpp
    tests/test_xenium.cpp
  )

  target_link_libraries(celladmix_core_tests PRIVATE celladmix_core)

  enable_testing()
  add_test(NAME celladmix_core_tests COMMAND celladmix_core_tests)
endif()

if (CELLADMIX_BUILD_BENCHMARKS)
  add_executable(xenium_load_bench
    benchmarks/xenium_load_bench.cpp
  )
  target_link_libraries(xenium_load_bench PRIVATE celladmix_core)
endif()
