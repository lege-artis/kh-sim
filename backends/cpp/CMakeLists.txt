# KH-SIM C++ backend — Task KH-005
# Requirements: CMake 3.20+, C++20 compiler (MSVC 19.29+ / GCC 11+ / Clang 13+)
# All dependencies fetched automatically via FetchContent (no vcpkg required).

cmake_minimum_required(VERSION 3.20)
project(kh-sim-cpp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Compiler flags ────────────────────────────────────────────────────────────
if(MSVC)
  add_compile_options(/W3 /O2 /EHsc)
  # Suppress cpp-httplib MSVC warnings
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS NOMINMAX WIN32_LEAN_AND_MEAN)
else()
  add_compile_options(-Wall -O2)
endif()

# ── Dependencies via FetchContent ─────────────────────────────────────────────
include(FetchContent)

# cpp-httplib — header-only HTTP server/client
FetchContent_Declare(
  httplib
  GIT_REPOSITORY https://github.com/yhirose/cpp-httplib.git
  GIT_TAG        v0.18.1
  GIT_SHALLOW    TRUE
)

# nlohmann/json — header-only JSON
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.3
  GIT_SHALLOW    TRUE
)

FetchContent_MakeAvailable(httplib nlohmann_json)

# ── Main server executable ────────────────────────────────────────────────────
add_executable(kh-sim-cpp
  src/main.cpp
  src/physics.cpp
)

target_include_directories(kh-sim-cpp PRIVATE src)

target_link_libraries(kh-sim-cpp PRIVATE
  httplib::httplib
  nlohmann_json::nlohmann_json
)

# Windows: link WinSock (required by cpp-httplib)
if(WIN32)
  target_link_libraries(kh-sim-cpp PRIVATE ws2_32)
endif()

# ── Validation test executable ────────────────────────────────────────────────
# Guarded: Docker build context only contains src/; tests/ is host-only.
# EXISTS check prevents CMake generate failure when tests/ is absent (Docker).
# CTest still runs normally on host where tests/validation_test.cpp is present.
enable_testing()
if(EXISTS "${CMAKE_SOURCE_DIR}/tests/validation_test.cpp")
  add_executable(kh-sim-cpp-validate
    tests/validation_test.cpp
    src/physics.cpp
  )
  target_include_directories(kh-sim-cpp-validate PRIVATE src)
  target_link_libraries(kh-sim-cpp-validate PRIVATE nlohmann_json::nlohmann_json)

  add_test(
    NAME validation
    COMMAND kh-sim-cpp-validate
            "${CMAKE_SOURCE_DIR}/../../shared/golden-vectors/v1-canonical-reference.json"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

# ── Golden-vector test executable (phase 6c) ──────────────────────────────────
# Mode 1: HTTP POST /simulate to running server on port 8003.
# Guarded by same EXISTS pattern as validation test (host-only; absent in Docker).
if(EXISTS "${CMAKE_SOURCE_DIR}/tests/test_golden_vectors.cpp")
  add_executable(test_golden_vectors
    tests/test_golden_vectors.cpp
  )
  target_include_directories(test_golden_vectors PRIVATE src)
  target_link_libraries(test_golden_vectors PRIVATE
    httplib::httplib
    nlohmann_json::nlohmann_json
  )
  if(WIN32)
    target_link_libraries(test_golden_vectors PRIVATE ws2_32)
  endif()

  add_test(
    NAME golden_vectors
    COMMAND test_golden_vectors
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

# ── Physics testbed test executable (phase 7d) ────────────────────────────────
# Mode 1: HTTP POST /simulate to running server on port 8003.
# Guarded by same EXISTS pattern (host-only; absent in Docker).
if(EXISTS "${CMAKE_SOURCE_DIR}/tests/test_physics_testbeds.cpp")
  add_executable(test_physics_testbeds
    tests/test_physics_testbeds.cpp
  )
  target_include_directories(test_physics_testbeds PRIVATE src)
  target_link_libraries(test_physics_testbeds PRIVATE
    httplib::httplib
    nlohmann_json::nlohmann_json
  )
  if(WIN32)
    target_link_libraries(test_physics_testbeds PRIVATE ws2_32)
  endif()

  add_test(
    NAME physics_testbeds
    COMMAND test_physics_testbeds
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

# ── Canonical-physics test executable (phase 7f) ──────────────────────────────
# Mode 1: HTTP POST /simulate to running server on port 8003.
# Guarded by same EXISTS pattern (host-only; absent in Docker).
# Every CP-7F testbed runs from a shared harness in tests/integration
# (test_canonical_physics_cp7fa.py, ..._cp7fbd.py, ..._cp7fc.py); C++ takes no part
# in CP-7F-E. This executable only prints where each one runs; one CTest entry.
if(EXISTS "${CMAKE_SOURCE_DIR}/tests/canonical_physics.cpp")
  add_executable(canonical_physics
    tests/canonical_physics.cpp
  )
  target_include_directories(canonical_physics PRIVATE src)
  target_link_libraries(canonical_physics PRIVATE
    httplib::httplib
    nlohmann_json::nlohmann_json
  )
  if(WIN32)
    target_link_libraries(canonical_physics PRIVATE ws2_32)
  endif()

  add_test(
    NAME canonical_physics_d
    COMMAND canonical_physics D
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()
