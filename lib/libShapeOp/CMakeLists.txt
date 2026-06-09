cmake_minimum_required(VERSION 3.25)
project(ShapeOp)

set(CMAKE_CONFIGURATION_TYPES Release Debug CACHE TYPE INTERNAL FORCE)
option(OPENMP TRUE)
  
# Tell CMake he'll be able to include ".cmake" configurations 
# files in the folder where the current file is located
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

# Qt
if ((NOT DEFINED ENV{QT_VERSION_TO_USE}) OR ($ENV{QT_VERSION_TO_USE} STREQUAL 5))
    find_package(Qt5 5.15 REQUIRED COMPONENTS Core Widgets Gui Concurrent Multimedia)
    find_package(Qt5Multimedia)
elseif($ENV{QT_VERSION_TO_USE} STREQUAL 6)
    find_package(Qt6 REQUIRED COMPONENTS Core Widgets Gui Concurrent Multimedia)
    #find_package(Qt6Multimedia)
else()
    message( FATAL_ERROR "ERROR: Invalid QT_VERSION_TO_USE envar provided (must be 5 or 6)" )
endif()
qt_wrap_cpp(CXX_FILES Solver.h)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
if(NOT EIGEN3_FOUND)
    message(ERROR " EIGEN not found.")
else()
    message("EIGEN found.")
endif() 

# Where are the headers found?
include_directories("src/")
include_directories("api/")

# h files
file(GLOB SHAPEOP_HEADERS "src/*.h")
# cpp files
file(GLOB SHAPEOP_SOURCES "src/*.cpp")
# h files
file(GLOB API_HEADERS "api/*.h")
# cpp files
file(GLOB API_SOURCES "api/*.cpp")

# Shared library
set(HEADERS ${SHAPEOP_HEADERS} ${API_HEADERS})
set(SOURCES ${SHAPEOP_SOURCES} ${API_SOURCES})
add_library(ShapeOp SHARED ${HEADERS} ${SOURCES})
set_target_properties(ShapeOp PROPERTIES VERSION ${SHAPEOP_VERSION_STRING} SOVERSION ${SHAPEOP_VERSION_MAJOR})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC Qt::Core Qt::Widgets Qt::Gui OpenMP::OpenMP_CXX Eigen3::Eigen)
endif()

# Rename executable to contain version number on windows. On mac, a so-chain is generated from the SOVERSION specified above.
if(CMAKE_SYSTEM_NAME MATCHES "Windows")
	set_target_properties(ShapeOp PROPERTIES OUTPUT_NAME ShapeOp.${SHAPEOP_VERSION_STRING})


endif()
