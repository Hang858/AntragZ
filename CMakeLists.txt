cmake_minimum_required(VERSION 3.10)
project(Zitaka_KeyGen C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(-O2 -Wall -Wextra -Werror=implicit-function-declaration)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DZITAKA_DEBUG)
    message(STATUS "Build Type: DEBUG (Verbose logging enabled)")
else()
    message(STATUS "Build Type: RELEASE (Logging disabled)")
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src               
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ntru
    ${CMAKE_CURRENT_SOURCE_DIR}/src/prematrix
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sampler
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utils
    ${CMAKE_CURRENT_SOURCE_DIR}/PQCDevKitAPI
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math/kiss_fft
)


file(GLOB PQC_SOFT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PQCDevKitAPI/soft/*.c)

set(SOURCES
    src/math/fft_clean.c
    src/utils/rng.c
    src/ntru/keygen_basis.c
    src/ntru/keygen_ntru.c
    src/prematrix/prematrix.c
    src/sampler/migd.c
    src/sampler/eigd.c
    src/sampler/sample.c
    src/math/cholesky.c
    src/keygen.c
    src/math/ntt.c
    src/math/zint.c
    src/math/poly.c
    src/math/poly_ntt.c
    src/math/kiss_fft/kiss_fft.c
    src/utils/file_utils.c
    src/utils/hash.c
    src/utils/encode.c
    src/sampler/zalcon_samp.c
    src/sign.c
    src/verify.c
    ${PQC_SOFT_SOURCES}
)

list(APPEND SOURCES test/test_sign_impl.c)

add_executable(zitaka_test ${SOURCES})


target_link_libraries(zitaka_test 
    m
)
