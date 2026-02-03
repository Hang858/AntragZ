cmake_minimum_required(VERSION 3.10)
project(Zitaka_KeyGen C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(-O2 -Wall -Wextra -Werror=implicit-function-declaration)


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include           # 新增：建议存放 common.h, api.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src               # 兼容旧引用
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ntru
    ${CMAKE_CURRENT_SOURCE_DIR}/src/prematrix
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sampler
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utils
    ${CMAKE_CURRENT_SOURCE_DIR}/PQCDevKitAPI
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math/kiss_fft
)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)

include_directories(${FFTW3_INCLUDE_DIRS})
link_directories(${FFTW3_LIBRARY_DIRS})

file(GLOB PQC_SOFT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PQCDevKitAPI/soft/*.c)

set(SOURCES
    src/math/fft_clean.c
    src/utils/rng.c
    src/ntru/keygen_basis.c
    src/ntru/keygen_ntru.c
    src/prematrix/prematrix.c
    src/sampler/migd.c
    src/sampler/eigd.c
    src/math/cholesky.c
    src/keygen.c
    src/math/ntt.c
    src/math/zint.c
    src/math/poly.c
    src/math/kiss_fft/kiss_fft.c
    ${PQC_SOFT_SOURCES}
)

list(APPEND SOURCES test/test_keygen.c)

# if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/main.c")
#     list(APPEND SOURCES main.c)
# elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/test_prematrix.c")
#     list(APPEND SOURCES test/test_prematrix.c)
# endif()

add_executable(zitaka_test ${SOURCES})


target_link_libraries(zitaka_test 
    # ${FFTW3_LIBRARIES} 
    m
)

# 输出构建成功信息
message(STATUS "FFTW3 Libraries: ${FFTW3_LIBRARIES}")
message(STATUS "Hardware Include Path set to: ${CMAKE_CURRENT_SOURCE_DIR}/../hardware_inc")