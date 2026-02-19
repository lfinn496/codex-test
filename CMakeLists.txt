cmake_minimum_required(VERSION 3.20)
project(mtt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(mtt
    src/geo.cpp
    src/pmbm_tracker.cpp
)

target_include_directories(mtt PUBLIC include)

target_compile_options(mtt PRIVATE -Wall -Wextra -Wpedantic)

add_executable(mtt_cli src/main.cpp)
target_link_libraries(mtt_cli PRIVATE mtt)
