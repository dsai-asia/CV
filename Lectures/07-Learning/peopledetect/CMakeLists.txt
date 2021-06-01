cmake_minimum_required (VERSION 2.8.12)

project (PEOPLEDETECT)

set (PEOPLEDETECT_VERSION_MAJOR 1)
set (PEOPLEDETECT_VERSION_MINOR 0)

set (CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}
                       /usr/local/opencv-3.1.1/share/OpenCV)

add_definitions (-Wall -Wno-unknown-pragmas)
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g0")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g0")
# Remove -rdynamic option, to allow stripping of symbols
set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)

find_package (OpenCV 3.4.1 EXACT REQUIRED)

add_executable (peopledetect ${CMAKE_CURRENT_SOURCE_DIR}/peopledetect.cpp)

target_include_directories (peopledetect SYSTEM PUBLIC ${OpenCV_INCLUDE_DIRS})

target_link_libraries (peopledetect ${OpenCV_LIBS})

