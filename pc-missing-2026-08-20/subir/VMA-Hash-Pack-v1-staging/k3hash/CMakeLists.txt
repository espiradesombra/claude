cmake_minimum_required(VERSION 3.15)
project(k3hash C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

option(K3HASH_BUILD_SHARED "Compilar como librería dinámica (.dll/.dylib/.so)" ON)
option(K3HASH_BUILD_EXAMPLES "Compilar ejemplos y CLI" ON)

if(K3HASH_BUILD_SHARED)
    add_library(k3hash SHARED src/k3hash.c)
    target_compile_definitions(k3hash PRIVATE K3HASH_BUILD_SHARED)
else()
    add_library(k3hash STATIC src/k3hash.c)
endif()

target_include_directories(k3hash PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# En macOS, fija el nombre de instalación de la dylib para que los
# binarios que la enlacen la encuentren correctamente.
if(APPLE AND K3HASH_BUILD_SHARED)
    set_target_properties(k3hash PROPERTIES
        MACOSX_RPATH ON
        INSTALL_NAME_DIR "@rpath"
    )
endif()

# En Windows, exporta todos los símbolos automáticamente si se usa MinGW;
# con MSVC se usan los __declspec(dllexport/dllimport) del propio header.
if(MINGW)
    set_target_properties(k3hash PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

set_target_properties(k3hash PROPERTIES
    VERSION 1.0.0
    SOVERSION 1
    PUBLIC_HEADER include/k3hash.h
)

if(K3HASH_BUILD_EXAMPLES)
    add_library(textutil STATIC examples/textutil.c)
    target_include_directories(textutil PUBLIC examples)

    add_executable(k3hash_cli examples/k3hash_cli.c)
    target_link_libraries(k3hash_cli PRIVATE k3hash)

    add_executable(k3dedup examples/k3dedup.c)
    target_link_libraries(k3dedup PRIVATE k3hash)

    add_executable(k3simil examples/k3simil.c)
    target_link_libraries(k3simil PRIVATE k3hash textutil)

    add_executable(k3search examples/k3search.c)
    target_link_libraries(k3search PRIVATE k3hash textutil)

    if(K3HASH_BUILD_SHARED)
        target_compile_definitions(k3hash_cli PRIVATE K3HASH_SHARED)
        target_compile_definitions(k3dedup PRIVATE K3HASH_SHARED)
        target_compile_definitions(k3simil PRIVATE K3HASH_SHARED)
        target_compile_definitions(k3search PRIVATE K3HASH_SHARED)
    endif()
endif()


# --- Instalación (opcional: cmake --install build) ---
include(GNUInstallDirs)
install(TARGETS k3hash
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION include
)
