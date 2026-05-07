cmake_minimum_required(VERSION 3.16)
project(log_bench CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -ffast-math")

find_path(SLEEF_INCLUDE_DIR sleef.h)
find_library(SLEEF_LIB sleef)
find_path(CRLIBM_INCLUDE_DIR crlibm.h)
find_library(CRLIBM_LIB crlibm)

include_directories(${SLEEF_INCLUDE_DIR} ${CRLIBM_INCLUDE_DIR} /usr/include)

add_executable(bench_log bench_log.cpp tu_algo.cpp)
target_link_libraries(bench_log ${SLEEF_LIB} ${CRLIBM_LIB} m)

add_executable(bench_root bench_root.cpp tu_algo.cpp)
target_link_libraries(bench_root ${SLEEF_LIB} ${CRLIBM_LIB} m)
