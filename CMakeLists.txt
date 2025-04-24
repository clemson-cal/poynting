cmake_minimum_required(VERSION 3.15)
project(poynting VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find Python and pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# Add compiler flags for different platforms
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Add optimization flags for Release builds
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        add_compile_options(/O2)
    else()
        add_compile_options(-O3)
    endif()
endif()

# Define the extension modules
pybind11_add_module(_core MODULE
    src/poynting/ext/core/_core.cpp
)

# Install the extension modules to the correct location
install(TARGETS _core
    DESTINATION ${SKBUILD_PROJECT_NAME}/ext/core
)

# Output configuration information
message(STATUS "Python version: ${Python_VERSION}")
message(STATUS "Python interpreter: ${Python_EXECUTABLE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
