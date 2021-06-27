cmake_minimum_required(VERSION 3.3)
cmake_policy(VERSION 3.3)

project(osprey-dpd VERSION "0.0.0")

find_package(MPI)

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

if(${CMAKE_VERSION} VERSION_GREATER "3.15.0") 
  target_precompile_headers(dpd PRIVATE src/StdAfx.h)
endif()


if ( MPI_FOUND )
include_directories(${MPI_INCLUDE_PATH})

# cheating until source file count goes down
file(GLOB SRC_FILES src/*.cpp)
add_executable(dpd-mpi ${SRC_FILES}  ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(dpd-mpi PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS} -DSimpMPS=2")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(dpd-mpi PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_compile_options(dpd-mpi
  PRIVATE ${COMPILE_OPTIONS}
)

if(${CMAKE_VERSION} VERSION_GREATER "3.15.0") 
  target_precompile_headers(dpd-mpi PRIVATE src/StdAfx.h)
endif()

endif()