cmake_minimum_required(VERSION 3.9)
cmake_policy(VERSION 3.9)

project(osprey-dpd VERSION "0.0.0")

include(CheckIPOSupported)

set(COMPILE_OPTIONS
    -g3
    -O3
    -fno-omit-frame-pointer
    -DNDEBUG=1

    -Wall
    -Wextra
    -pedantic
    #-Werror
    -Wno-unused-variable
    -Wno-unused-parameter
    -Wno-reorder
    -Wno-sign-compare
    )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# cheating until source file count goes down
file(GLOB SRC_FILES src/*.cpp src/*/*.cpp)

INCLUDE_DIRECTORIES(
  src
  src/engines
)


find_package(TBB)
if(TBB_FOUND)
set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -DOSPREY_DPD_ENABLE_TBB=1 )
message("TBB found and enabled")
link_libraries(TBB::tbb)
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -DOSPREY_DPD_ENABLE_OPENMP=1 )
message("OpenMP found and enabled")
endif()

if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -mavx2)
endif()


add_executable(dpd ${SRC_FILES})
target_compile_options(dpd
  PRIVATE ${COMPILE_OPTIONS}  -DOSPREY_DPD_ENABLE_STATE_LOGGING=1
)

add_executable(dpd-bd ${SRC_FILES})
target_compile_options(dpd-bd
  PRIVATE ${COMPILE_OPTIONS} -DSimIdentifier=BD
)

add_executable(dpd-md ${SRC_FILES})
target_compile_options(dpd-md
  PRIVATE ${COMPILE_OPTIONS} -DSimIdentifier=MD
)

add_executable(dpd-2d ${SRC_FILES})
target_compile_options(dpd-2d
  PRIVATE ${COMPILE_OPTIONS} -DSimDimension=2
)


#############################################################
## Debug build

add_executable(dpd-debug ${SRC_FILES})
target_compile_options(dpd-debug
  PRIVATE ${COMPILE_OPTIONS}
  -DOSPREY_DPD_ENABLE_STATE_LOGGING=1
  -DOSPREY_ENABLE_DEBUG_ASSERT=1
  -DTraceOn
  -UNDEBUG
)
target_compile_options(dpd-debug PRIVATE -fsanitize=undefined -fsanitize=address)
target_link_options(dpd-debug PRIVATE -fsanitize=undefined -fsanitize=address)

#############################################################
## Most optimised version.

add_executable(dpd-opt ${SRC_FILES})
target_compile_options(dpd-opt
  PRIVATE ${COMPILE_OPTIONS}
)
# These get a good boost, and should be portable. fast-math appears fine in the 
# context of commands I've tried, but there may be parts where it doesn't work.
# AFAICT nothing in osprey depends on errno, or has very sensitive ordering.
target_compile_options(dpd-opt PRIVATE -fomit-frame-pointer -ffast-math)
target_link_options(dpd-opt PRIVATE -fomit-frame-pointer -ffast-math)

# -march=native may cause problems when compiling on head nodes of a different type
# to cluster nodes (e.g. head node has avx512, cluster nodes don't), but that is
# usually not the case.
# Note that -march=native doesn't work on Apple-M1+clang till version 15, which
# is not (as of May 2023) the default apple clang, so we only enable on gcc.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
target_compile_options(dpd-opt PRIVATE -march=native)
target_link_options(dpd-opt PRIVATE -march=native)
endif()

# Enable link-time optimisation if possible. This may substantially increase
# compilation time
check_ipo_supported(RESULT check_ipo_supported_result OUTPUT output)
if(check_ipo_supported_result)
  set_property(TARGET dpd-opt PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
