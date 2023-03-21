cmake_minimum_required(VERSION 3.3)
cmake_policy(VERSION 3.3)

project(osprey-dpd VERSION "0.0.0")

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

set(COMPILE_OPTIONS_DEBUG
  -g3
  -O0

  -Wall
  -Wextra
  -pedantic
  #-Werror
  -Wno-unused-variable
  -Wno-unused-parameter
  -Wno-reorder
  -Wno-sign-compare

  -DTraceOn
  -UNDEBUG
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# cheating until source file count goes down
file(GLOB SRC_FILES src/*.cpp src/*/*.cpp)

INCLUDE_DIRECTORIES(
  src
  src/aa
  src/ac
  src/acf
  src/ae
  src/aev 
  src/cc
  src/cgc
  src/cgt
  src/ct
  src/ep_protocols
  src/et_tools
  src/ev_events
  src/interfaces
  src/log
  src/mc_monitor_cmds
  src/mp_parallel 
  src/pc 
  src/pm
  src/pr_processes
  src/ta_decorators
  src/tgu
)


add_executable(dpd ${SRC_FILES})
target_compile_options(dpd
  PRIVATE ${COMPILE_OPTIONS}
)
# These get a good boost, and should be portable. fast-math appears fine in the 
# context of commands I've tried, but there may be parts where it doesn't work.
# AFAICT nothing in osprey depends on errno, or has very sensitive ordering.
target_compile_options(dpd PRIVATE -march=native -fomit-frame-pointer -ffast-math)
target_link_options(dpd PRIVATE -march=native -fomit-frame-pointer -ffast-math)


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

add_executable(dpd-debug ${SRC_FILES})
target_compile_options(dpd-debug
  PRIVATE ${COMPILE_OPTIONS_DEBUG}
)
target_compile_options(dpd-debug PRIVATE -fsanitize=address)
target_link_options(dpd-debug PRIVATE -fsanitize=address)