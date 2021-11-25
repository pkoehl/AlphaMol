cmake_minimum_required(VERSION 3.1.0)
project(alphamol)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Try to find the GNU Multiple Precision Arithmetic Library (GMP)
find_path(GMP_INCLUDE_DIR gmp.h
    PATHS /usr/include /usr/local/include )

find_library(GMP_LIBRARY NAMES gmp
    PATHS /usr/lib /usr/local/lib)

if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
    get_filename_component(GMP_LIBRARY_DIR ${GMP_LIBRARY} PATH)
    set(GMP_FOUND TRUE)
endif()

if(GMP_FOUND)
   if(NOT GMP_FIND_QUIETLY)
      MESSAGE(STATUS "Found GMP: ${GMP_LIBRARY}")
   endif()
elseif(GMP_FOUND)
   if(GMP_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find GMP")
   endif()
endif()

# Include JLCxx
find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

# AlphaMol includes
list(APPEND ALPHAMOL_INCLUDE_DIRS "DataStruct/include")
list(APPEND ALPHAMOL_INCLUDE_DIRS "Delcx/include")
list(APPEND ALPHAMOL_INCLUDE_DIRS "Alphacx/include")
list(APPEND ALPHAMOL_INCLUDE_DIRS "InOut/include")
list(APPEND ALPHAMOL_INCLUDE_DIRS "Volumes/include")
list(APPEND ALPHAMOL_INCLUDE_DIRS "project/include")


# build app
file(GLOB ALPHAMOL_SOURCE "project/src/*.cpp")
include_directories(${ALPHAMOL_INCLUDE_DIRS} ${GMP_INCLUDE_DIR})
add_library(AlphaMol SHARED ${ALPHAMOL_SOURCE})
target_link_libraries(AlphaMol gmp JlCxx::cxxwrap_julia)
