
cmake_minimum_required(VERSION 3.12)

project(homography)

set(CMAKE_CXX_STANDARD 14)

add_executable(homography hw2.cpp)

find_package(OpenCV 4.0.1 REQUIRED)

target_include_directories(homography PUBLIC ${OpenCV_INCLUDE_DIRS})

target_link_libraries(homography ${OpenCV_LIBS})

