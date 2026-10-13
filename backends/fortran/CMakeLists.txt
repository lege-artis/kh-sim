# KH-SIM Fortran backend — Task KH-006
# Mixed-language: Fortran physics + C++ HTTP shim (cpp-httplib)
# Requirements: CMake 3.20+, gfortran 11+, g++ 11+, Ninja
# FetchContent: cpp-httplib, nlohmann/json (no manual install needed)

cmake_minimum_required(VERSION 3.20)
project(kh-sim-fortran VERSION 0.1.0 LANGUAGES Fortran CXX)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Compiler flags ─────────────────────────────────────────────────────────────
# Note: each flag must be a separate argument inside generator expressions
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options(
    $<$<COMPILE_LANGUAGE:Fortran>:-O2>
    $<$<COMPILE_LANGUAGE:Fortran>:-fcheck=bounds>
    $<$<COMPILE_LANGUAGE:Fortran>:-fbacktrace>
  )
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(
    $<$<COMPILE_LANGUAGE:CXX>:-O2>
    $<$<COMPILE_LANGUAGE:CXX>:-Wall>
  )
endif()
add_compile_definitions(WIN32_LEAN_AND_MEAN NOMINMAX _CRT_SECURE_NO_WARNINGS)

# ── FetchContent dependencies ──────────────────────────────────────────────────
include(FetchContent)

FetchContent_Declare(
  httplib
  GIT_REPOSITORY https://github.com/yhirose/cpp-httplib.git
  GIT_TAG        v0.18.1
  GIT_SHALLOW    TRUE
)
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.3
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(httplib nlohmann_json)

# ── Fortran physics module ─────────────────────────────────────────────────────
# Compiled first; C++ shim links against the resulting .mod + .o files
add_library(kh_physics_f STATIC
  src/kh_constants.f90
  src/kh_fft.f90
  src/kh_grid.f90
  src/kh_initial_conditions.f90
  src/kh_nonlinear.f90
  src/kh_diagnostics.f90
  src/kh_poisson.f90
  src/kh_velocity.f90
  src/kh_forcing.f90
  src/kh_physics.f90
)
target_compile_options(kh_physics_f PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-O2>
)

# ── C++ HTTP shim (calls Fortran via iso_c_binding exports) ───────────────────
add_library(kh_shim STATIC
  src/kh_shim.cpp
)
target_include_directories(kh_shim PRIVATE src)
target_link_libraries(kh_shim PRIVATE
  kh_physics_f
  httplib::httplib
  nlohmann_json::nlohmann_json
)

# ── Main server executable ─────────────────────────────────────────────────────
# main.cpp directly includes <httplib.h> and <nlohmann/json.hpp>, so those
# interface targets must be linked here too (not just transitively via kh_shim).
add_executable(kh-sim-fortran src/main.cpp)
target_include_directories(kh-sim-fortran PRIVATE src)
target_link_libraries(kh-sim-fortran PRIVATE
  kh_shim
  httplib::httplib
  nlohmann_json::nlohmann_json
)
if(WIN32)
  target_link_libraries(kh-sim-fortran PRIVATE ws2_32)
endif()

# ── Validation test ────────────────────────────────────────────────────────────
# Guarded: Docker build context only contains src/; tests/ is host-only.
# EXISTS check prevents CMake generate failure when tests/ is absent (Docker).
# CTest still runs normally on host where tests/validation_test.cpp is present.
enable_testing()
if(EXISTS "${CMAKE_SOURCE_DIR}/tests/validation_test.cpp")
  add_executable(kh-sim-fortran-validate tests/validation_test.cpp)
  target_include_directories(kh-sim-fortran-validate PRIVATE src)
  target_link_libraries(kh-sim-fortran-validate PRIVATE
    kh_shim
    nlohmann_json::nlohmann_json
  )
  if(WIN32)
    target_link_libraries(kh-sim-fortran-validate PRIVATE ws2_32)
  endif()

  add_test(
    NAME validation
    COMMAND kh-sim-fortran-validate
            "${CMAKE_SOURCE_DIR}/../../shared/golden-vectors/v1-canonical-reference.json"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()
