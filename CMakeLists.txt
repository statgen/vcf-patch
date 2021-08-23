cmake_minimum_required(VERSION 3.2)
project(vcf-patch)

set(CMAKE_CXX_STANDARD 11)

find_package(savvy REQUIRED)

#find_path(HTS_INCLUDE_DIRS htslib/bgzf.h HINTS ${CMAKE_SOURCE_DIR}/../htslib)
#if(NOT HTS_INCLUDE_DIRS )
#    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
#endif()
#message(STATUS "HTS include: ${HTS_INCLUDE_DIRS}")
#
#
#find_library(HTS_LIBRARIES hts HINTS ${CMAKE_SOURCE_DIR}/../htslib)
#if(NOT HTS_LIBRARIES)
#    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
#endif()
#message(STATUS "HTS lib: ${HTS_LIBRARIES}")
#
#find_library(ZLIB z)
#if(NOT ZLIB)
#    message(FATAL_ERROR "libz library not found")
#endif()

add_executable(vcf-patch main.cpp)
#target_include_directories(vcf-patch PRIVATE ${HTS_INCLUDE_DIRS})
target_link_libraries(vcf-patch savvy)
