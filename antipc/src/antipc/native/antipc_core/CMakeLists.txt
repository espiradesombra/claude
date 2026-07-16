cmake_minimum_required(VERSION 3.16)
project(antipc_native C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(K3HASH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../k3hash")

add_library(antipc_native SHARED
    src/antipc_native.c
    src/mdc_factor.c
    src/sieve.c
    src/criva_est.c
    src/k3_stream.c
    src/geo_convergence.c
    src/mdc_trains.c
    src/criba_hibrida.c
    src/criba_modular6k.c
    src/newton_rapido.c
    src/mdc_ksweep.c
    src/k3_geo_aleatorovix.c
    src/k3_suite.c
    ${K3HASH_DIR}/src/k3hash.c
)

target_include_directories(antipc_native PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${K3HASH_DIR}/include
)

target_compile_definitions(antipc_native PRIVATE
    ANTIPC_BUILD_SHARED
    K3HASH_BUILD_SHARED
)

if(MSVC)
    target_compile_options(antipc_native PRIVATE /O2 /W3)
else()
    target_compile_options(antipc_native PRIVATE -O3 -Wall)
    target_link_libraries(antipc_native m)
endif()

set_target_properties(antipc_native PROPERTIES
    OUTPUT_NAME "antipc_native"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Release"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Release"
)