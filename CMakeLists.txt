cmake_minimum_required(VERSION 3.24)
set(PROJECT_LANGUAGES "Fortran")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_FORTRAN_STANDARD 90)
#set(CMAKE_FORTRAN_STANDARD_REQUIRED ON)
set(FORTRAN_INTEGER_SIZE_4B "")
set(FORTRAN_INTEGER_SIZE_8B "-fdefault-integer-8")
# fixed form, C preprocessed, all REAL as double, 4-byte integers (int32_t, int)
set(FORTRAN_INTEGER_SIZE ${FORTRAN_INTEGER_SIZE_4B})
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-form -cpp -fdefault-real-8 ${FORTRAN_INTEGER_SIZE} -fmax-stack-var-size=1024")
#add_compile_definitions(foo=42)
set(PROJECT_FAMILY_NAME "CubeSphere")
set(PROJECT_VERSION "0.5.0")
set(PROJECT_VERSION_NAME "RoselandRiemann")
set(PROJECT_VERSION_MONTH "August")
set(PROJECT_VERSION_YEAR "2026")
set(PROJECT_WORKING_RDIRS "_T" "_I")

set(CUSTOM_FILELIST
        # include this if building library to keep CMake from emitting spurious error messages.
        _L/Lmain.f
#    _L/001addmul/Laddmul.f
#    _L/002vert_coords_patch/Lvert_coords_patch.f
    _L/200mesh/Lmesh.f
)


include($ENV{DGCALC_HOME}/src/dgcalc.1.cmake)
set(PROJECT_NAME "${PROJECT_FAMILY_NAME}${PROJECT_VERSION}${PROJECT_VERSION_NAME}_${PROJECT_VERSION_MONTH}${PROJECT_VERSION_YEAR}")
project(${PROJECT_NAME} VERSION ${PROJECT_VERSION} LANGUAGES ${PROJECT_LANGUAGES})
include($ENV{DGCALC_HOME}/src/dgcalc.2.cmake)


