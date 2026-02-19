cmake_minimum_required(VERSION 3.16)
project(mtt_tracker VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ─── Optimisation flags ────────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -march=native -DNDEBUG)
endif()

# ─── Required: Eigen3 ─────────────────────────────────────────────────────────
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# ─── Optional: GeographicLib ──────────────────────────────────────────────────
# If found, link against it for production-grade geodetic ops.
# If absent, the built-in Bowring/Heikkinen implementation is used instead.
find_package(GeographicLib QUIET)

# ─── Threads ──────────────────────────────────────────────────────────────────
find_package(Threads REQUIRED)

# ─── Source files ─────────────────────────────────────────────────────────────
set(MTT_SOURCES
    src/geodesy.cpp
    src/main.cpp
)

add_executable(mtt_tracker ${MTT_SOURCES})

target_include_directories(mtt_tracker
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(mtt_tracker
    PRIVATE
        Eigen3::Eigen
        Threads::Threads
)

if(GeographicLib_FOUND)
    target_link_libraries(mtt_tracker PRIVATE ${GeographicLib_LIBRARIES})
    target_include_directories(mtt_tracker PRIVATE ${GeographicLib_INCLUDE_DIRS})
    target_compile_definitions(mtt_tracker PRIVATE MTT_USE_GEOGRAPHICLIB)
    message(STATUS "GeographicLib found – using for geodetic transforms")
else()
    message(STATUS "GeographicLib NOT found – using built-in Bowring/Heikkinen geodesy")
endif()

# ─── Warning flags ────────────────────────────────────────────────────────────
target_compile_options(mtt_tracker PRIVATE
    -Wall -Wextra -Wpedantic
    -Wno-unused-parameter
)

# ─── Install ──────────────────────────────────────────────────────────────────
install(TARGETS mtt_tracker DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

# ─── Simple unit-test target (no external framework) ─────────────────────────
option(MTT_BUILD_TESTS "Build unit tests" ON)
if(MTT_BUILD_TESTS)
    add_executable(mtt_test tests/test_main.cpp src/geodesy.cpp)
    target_include_directories(mtt_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(mtt_test PRIVATE Eigen3::Eigen Threads::Threads)
    target_compile_options(mtt_test PRIVATE -Wall -Wextra)
    enable_testing()
    add_test(NAME mtt_unit_tests COMMAND mtt_test)
endif()
