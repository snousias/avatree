cmake_minimum_required (VERSION 2.8.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
set(TARGETWEBAPP LungModellingWebInterface CACHE TYPE STRING)
add_definitions( -DSRC_PATH="${CMAKE_SOURCE_DIR}" )

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR} CACHE TYPE PATH)


find_package(Threads REQUIRED)
set(BOOST_COMPONENTS system thread filesystem date_time)
# Late 2017 TODO: remove the following checks and always use std::regex
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        set(BOOST_COMPONENTS ${BOOST_COMPONENTS} regex)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST_REGEX")
    endif()
endif()



include_directories(.)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/Modeller/Modeller/include)
include_directories(${EIGEN3_INCLUDE_DIR})
	include_directories(${Boost_INCLUDE_DIRS}) 
	include_directories(${CMAKE_SOURCE_DIR}/include/OpenGL)
	include_directories(${CMAKE_SOURCE_DIR}/include)
	include_directories(${CMAKE_SOURCE_DIR}/src/Modeller/Modeller)
	include_directories(${CMAKE_SOURCE_DIR}/src/Modeller/Modeller/include)
	include_directories(${CGAL_DIR}/auxiliary/gmp/include)
	include_directories(${CGAL_BUILD_DIR}/include)
	include_directories(${CGAL_DIR}/include)
	
	link_directories(${CGAL_BUILD_DIR}/lib)
	link_directories(${CGAL_BUILD_DIR}/lib/$(Configuration))




if(APPLE)
    set(OPENSSL_ROOT_DIR "/usr/local/opt/openssl")
endif()


link_directories(${Boost_LIBRARY_DIR})
link_directories(${CMAKE_BINARY_DIR})
link_directories(${CMAKE_SOURCE_DIR}/lib)



file(GLOB SOURCE_FILES 
*.cpp
modellerToViewer.h
)


add_executable(
${TARGETWEBAPP} 
${SOURCE_FILES}  
${CMAKE_CURRENT_SOURCE_DIR}/images/favicon.ico)


target_link_libraries(${TARGETWEBAPP} ${Boost_LIBRARIES})
target_link_libraries(${TARGETWEBAPP} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(${TARGETWEBAPP} ModellerLib )
set_target_properties(${TARGETWEBAPP}        
        PROPERTIES        
        PROJECT_LABEL "${TARGETWEBAPP}"  
        FOLDER "WebServer"      
)


#TODO: add requirement for version 1.0.1g (can it be done in one line?)
#find_package(OpenSSL)
#if(OPENSSL_FOUND)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_OPENSSL")
#    target_link_libraries(http_examples ${OPENSSL_LIBRARIES})
#    include_directories(SYSTEM ${OPENSSL_INCLUDE_DIR})
#    add_executable(https_examples https_examples.cpp)
#    target_link_libraries(https_examples ${Boost_LIBRARIES})
#    target_link_libraries(https_examples ${OPENSSL_LIBRARIES})
#    target_link_libraries(https_examples ${CMAKE_THREAD_LIBS_INIT})
#endif()
#if(MSYS) #TODO: Is MSYS true when MSVC is true?
#    target_link_libraries(http_examples ws2_32 wsock32)
#    if(OPENSSL_FOUND)
#        target_link_libraries(https_examples ws2_32 wsock32)
#    endif()
#endif()

enable_testing()
#add_subdirectory(tests)

install(FILES server_http.hpp client_http.hpp server_https.hpp client_https.hpp DESTINATION include/WebServer)
