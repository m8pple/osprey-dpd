cmake_minimum_required(VERSION 3.9)
cmake_policy(VERSION 3.9)

project(osprey-dpd VERSION "0.0.0")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

    -DOSPREY_DPD_ENABLE_STATE_LOGGING=1
    )

include(CheckIPOSupported)

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

set(COMPILE_OPTIONS_DEBUG
  ${COMPILE_OPTIONS}
  -DOSPREY_ENABLE_DEBUG_ASSERT=1
  -DTraceOn
  -UNDEBUG
  -fsanitize=undefined
  -fsanitize=address
  -O0
)
set(LINK_OPTIONS_DEBUG
  -fsanitize=undefined
  -fsanitize=address
)


set(COMPILE_OPTIONS_RELEASE
  ${COMPILE_OPTIONS}
  -UOSPREY_DPD_ENABLE_STATE_LOGGING
  -UOSPREY_ENABLE_DEBUG_ASSERT
  -UTraceOn
)

INCLUDE_DIRECTORIES(
  src
  src/engines
)

#####################################################################
## Create "virtual" object libraries
##
## These are used to reduce the recompilation time when there are both
## programs and unit tests that depend on the main body of the code.
## 
## They must be virtual libraries (OBJECT) rather than shared libraries
## due to the self-registration pattern used in osprey. It relies on
## static initialisers running before main, but if linking with a 
## static library the translation units containing initialisers will
## be stripped out.

# cheating until source file count goes down
file(GLOB SRC_FILES src/*.cpp src/*/*.cpp)

# Get the list with everything except source containing main
file(GLOB SRC_FILES_NO_MAIN src/*.cpp src/*/*.cpp)
list(FILTER SRC_FILES_NO_MAIN EXCLUDE REGEX ".*/src/dmpc.cpp")

# Create a library that contains the default-compiled version of everything _except_ dmpc.cpp)
add_library(libdpd OBJECT ${SRC_FILES_NO_MAIN})
target_compile_options(libdpd PRIVATE ${COMPILE_OPTIONS})

# Create a library that contains the debug default-compiled version of everything _except_ dmpc.cpp)
add_library(libdpd-debug OBJECT ${SRC_FILES_NO_MAIN})
target_compile_options(libdpd-debug PRIVATE ${COMPILE_OPTIONS_DEBUG})

# Create a library that contains the release-compiled version of everything _except_ dmpc.cpp)
add_library(libdpd-release OBJECT ${SRC_FILES_NO_MAIN})
target_compile_options(libdpd-release PRIVATE ${COMPILE_OPTIONS_RELEASE})

#############################################################
## Main build targets

# Bog standard version. Nothing special
add_executable(dpd src/dmpc.cpp)
target_link_libraries(dpd libdpd)
target_compile_options(dpd
  PRIVATE ${COMPILE_OPTIONS}
)

# Debug version of dpd with limited sanitizers turned on
add_executable(dpd-debug  src/dmpc.cpp)
target_link_libraries(dpd-debug libdpd-debug)
target_compile_options(dpd-debug
  PRIVATE ${COMPILE_OPTIONS_DEBUG} 
)
target_link_options(dpd-debug PRIVATE ${LINK_OPTIONS_DEBUG})

# Release version. Max optimisation, min error checking
add_executable(dpd-release src/dmpc.cpp)
target_link_libraries(dpd-release libdpd-release)
target_compile_options(dpd-release
  PRIVATE ${COMPILE_OPTIONS_RELEASE}
)


################################################################################
## These are mainly added to check less common options still compile

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
## Most optimised version.

add_executable(dpd-opt ${SRC_FILES})
target_compile_options(dpd-opt
  PRIVATE ${COMPILE_OPTIONS_RELEASE}
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


############################################################################
## Unit tests

function(add_unittest cpp_name)

add_executable(${cpp_name} testing/unittests/${cpp_name}.cpp)

target_link_libraries(${cpp_name} libdpd)
target_compile_options(${cpp_name}
  PRIVATE ${COMPILE_OPTIONS}  -DOSPREY_DPD_ENABLE_STATE_LOGGING=1
)
add_test(NAME ${cpp_name} COMMAND ${cpp_name})

add_executable(${cpp_name}-debug testing/unittests/${cpp_name}.cpp)
target_link_libraries(${cpp_name}-debug libdpd-debug)
target_compile_options(${cpp_name}-debug
  PRIVATE ${COMPILE_OPTIONS_DEBUG} 
)
target_link_options(${cpp_name}-debug
  PRIVATE  ${LINK_OPTIONS_DEBUG}
)

endfunction()


add_unittest(test_parallel_services)
