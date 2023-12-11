cmake_minimum_required(VERSION 3.10)
project(MyProject)

option(_MPI "option for MPI" ON)
if(_MPI)
    add_compile_definitions(_MPI)
endif()

# 添加头文件路径
include_directories(include)

# 添加源代码
add_executable(myapp src/main.cpp src/grid.cpp src/input.cpp src/chazhi.cpp src/timer.cpp src/scalapack_diago.cpp)

# 添加ScaLAPACK库的路径
#find_package(SCALAPACK REQUIRED)
#include_directories(${SCALAPACK_INCLUDE_DIRS})
#link_directories(${SCALAPACK_LIBRARY_DIRS})

# 链接ScaLAPACK库
#target_link_libraries(myapp PRIVATE ${SCALAPACK_LIBRARIES} -llapacke)


# 查找并链接LAPACK库
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)


# 如果 lapacke.h 头文件不在默认的搜索路径中，需要手动添加
include_directories(${LAPACK_INCLUDE_DIRS})
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})


# 链接所有依赖库
target_link_libraries(myapp PRIVATE ${LAPACK_LIBRARIES} ${MPI_CXX_LIBRARIES} OpenMP::OpenMP_CXX -llapacke libscalapack-openmpi.so)


