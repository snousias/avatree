
# Minimum required version of CMake
cmake_minimum_required(VERSION 2.8)


set(TARGETWEBSRV EmbeddedWebServer CACHE TYPE STRING)
project(${TARGETWEBSRV})
file(GLOB MODELLEREXEC_SRC
	"*.cpp"
	"*.h"
	)

	
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

set(CIVETWEB_INCLUDE_DIR "" CACHE TYPE STRING)
set(CIVETWEB_LIB_DIR "" CACHE TYPE STRING)

include_directories(${CIVETWEB_INCLUDE_DIR})
link_directories(${CIVETWEB_LIB_DIR})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR} )
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR} )

add_executable(${TARGETWEBSRV} ${MODELLEREXEC_SRC})
target_link_libraries(${TARGETWEBSRV} civetweb civetweb-cpp)



