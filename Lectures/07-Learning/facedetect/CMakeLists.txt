cmake_minimum_required (VERSION 2.8.12)

project (FACEDETECT)

set (FACEDETECT_VERSION_MAJOR 1)
set (FACEDETECT_VERSION_MINOR 0)

set (CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}
                       /usr/local/opencv-3.2.0/share/OpenCV)

add_definitions (-Wall -Wno-unknown-pragmas)
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g0")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g0")
# Remove -rdynamic option, to allow stripping of symbols
set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)

find_package (OpenCV 3.2.0 EXACT REQUIRED)

add_executable (facedetect ${CMAKE_CURRENT_SOURCE_DIR}/facedetect.cpp)

target_include_directories (facedetect SYSTEM PUBLIC ${OpenCV_INCLUDE_DIRS})

target_link_libraries (facedetect ${OpenCV_LIBS})

