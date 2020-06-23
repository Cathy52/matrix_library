cmake_minimum_required(VERSION 3.6)
PROJECT(project)

set(OpenCV_INCLUDE_DIRS D:/CS_Tools/opencv/mingw_build/install/include)
FIND_PACKAGE(opencv REQUIRED)
include_directories(${PROJECT_SOURCE_DIR} "D:/CS_Tools/opencv/mingw_build/install/include/opencv2")

# 生成可执行文件
add_executable(test main.cpp matrix.hpp SpMatrix.hpp my_exception.cpp)
target_link_libraries(test ${OpenCV_LIBS})