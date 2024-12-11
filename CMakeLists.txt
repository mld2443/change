cmake_minimum_required(VERSION 3.27)
project(makechange)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi -fcx-fortran-rules -Wall -Wextra -Wpedantic")

add_executable(${PROJECT_NAME} change.cpp)
