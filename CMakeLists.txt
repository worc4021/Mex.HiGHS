
set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake CACHE STRING "Vcpkg toolchain file")
cmake_minimum_required(VERSION 3.24)
project(highs_mex)

set(BUILD_SHARED_LIBS OFF)
find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM)
find_package(highs CONFIG REQUIRED)

SET(USE_EIGEN OFF)
SET(BUILD_TESTS OFF)
add_subdirectory(MexUtilities)

matlab_add_mex(
    NAME highs_lp
    SRC src/linprog.cpp
    OUTPUT_NAME highs_lp
    LINK_TO 
        MexUtilities::MexUtilities  
        highs::highs
    R2018a
)

set_property(TARGET highs_lp PROPERTY
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

matlab_add_mex(
    NAME highs_qp
    SRC src/quadprog.cpp
    OUTPUT_NAME highs_qp
    LINK_TO 
        MexUtilities::MexUtilities  
        highs::highs
    R2018a
)

set_property(TARGET highs_qp PROPERTY
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

include(CTest)
matlab_add_unit_test(
    NAME highs_test
    UNITTEST_FILE highstest.m
    ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/matlab $<TARGET_FILE_DIR:highs_lp>
)