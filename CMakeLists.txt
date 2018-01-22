cmake_minimum_required(VERSION 2.8.9)
project(relz)
set(CMAKE_BUILD_TYPE Release)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O3")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -O3")
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR}/include)
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")

add_executable(compress "${CMAKE_SOURCE_DIR}/programs/compress.cpp")
add_executable(decompress "${CMAKE_SOURCE_DIR}/programs/decompress.cpp")
add_executable(test_rmq "${CMAKE_SOURCE_DIR}/programs/test_rmq.cpp")
add_executable(prueba_search "${CMAKE_SOURCE_DIR}/programs/prueba_search.cpp")

add_library(relz SHARED ${SOURCES})
link_directories(${CMAKE_BINARY_DIR})

target_link_libraries(compress relz pthread)
target_link_libraries(decompress relz pthread)
target_link_libraries(test_rmq relz pthread)
target_link_libraries(prueba_search relz pthread)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "/usr/" CACHE PATH "default install path" FORCE )
endif()

install(TARGETS relz DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(TARGETS compress DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/relz FILES_MATCHING PATTERN "*.h")

