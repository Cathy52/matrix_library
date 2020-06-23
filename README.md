#### 一 简介

简单的矩阵运算库，以及 OpenCV 调用运算库做简单的图片处理

#### 二 开发环境

CLion + mingw64

#### 三 Run Demo

1.修改 CMakeList.txt 文件中的 OpenCV 安装地址，且需要配置 OpenCV_DIR（大小写要写对） 的环境变量

参考：

```C++
cmake_minimum_required(VERSION 3.6)
PROJECT(project)

set(OpenCV_INCLUDE_DIRS D:/CS_Tools/opencv/mingw_build/install/include)
FIND_PACKAGE(opencv REQUIRED)
include_directories(${PROJECT_SOURCE_DIR} "D:/CS_Tools/opencv/mingw_build/install/include/opencv2")

# 生成可执行文件
add_executable(test main.cpp matrix.hpp SpMatrix.hpp my_exception.cpp)
target_link_libraries(test ${OpenCV_LIBS})
```

环境变量：

Open_CV: opencv/mingw_build/install

2.直接运行 main.cpp 中的方法

#### 四 文件说明

Test image: img1.jpg img2.jpg

Source file: matrix.hpp my_exception.cpp SpMatrix.hpp

Other: CMakeList.txt




