cmake_minimum_required(VERSION 3.3)
cmake_policy(VERSION 3.3)

project(osprey-dpd VERSION "0.0.0")

set(COMPILE_OPTIONS
    -g3
    -O3
    -fno-omit-frame-pointer

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
file(GLOB SRC_FILES src/*.cpp)
add_executable(dpd ${SRC_FILES})

target_compile_options(dpd
  PRIVATE ${COMPILE_OPTIONS}
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
