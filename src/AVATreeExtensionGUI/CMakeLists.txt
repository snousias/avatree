cmake_minimum_required(VERSION 2.8.12)
cmake_policy(SET CMP0020 NEW) # for Qt: auto-link to qtmain.lib
set(TARGETAPP_AVATreeExtensionGUI AVATreeExtensionGUI CACHE TYPE STRING)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR} CACHE TYPE PATH)

#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
file(GLOB SOURCE_FILES *.cpp)
file(GLOB HEADER_FILES *.h)
file(GLOB UI_FILES *.ui)
file(GLOB SHADER_FILES Shaders/*.glsl)
set(Srcs ${SOURCE_FILES})
set(Hdrs ${HEADER_FILES})
set(MOC_Hdrs ${HEADER_FILES})
set(UIs ${UI_FILES})
set(SHs ${SHADER_FILES})
set(QRCs Icons/icons.qrc Shaders/shaders.qrc)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR} CACHE TYPE PATH)





####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

if(UNIX)
#Include
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/src/Modeller/Modeller)
include_directories(${CMAKE_SOURCE_DIR}/src/Modeller/Modeller/include)
include_directories(${EIGEN3_INCLUDE_DIR})
#Link Directories
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
# The rest should just work for Qt5 sure!
# We have ui files, this will bring in the macro: qt5_wrap_ui
qt5_wrap_ui(UI_Srcs ${UIs})
# Wrap resources
qt5_add_resources(QRC_Srcs ${QRCs})
# Group resources
source_group("Resources" FILES
  ${UIs}
  ${QRCs}
  ${EXE_ICON} # Not present
)
source_group("Generated" FILES
  ${UI_Srcs}
  ${QRC_Srcs}
)
source_group("Shaders" FILES
  ${SHs}
)
file(COPY ${CMAKE_SOURCE_DIR}/dependencies DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/images DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/models DESTINATION ${CMAKE_BINARY_DIR})

file(GLOB STATIC_LIBS_COPIED ${CMAKE_BINARY_DIR}/dependencies/*.a)
link_directories(/usr/lib/)
add_executable(${TARGETAPP_AVATreeExtensionGUI} MACOSX_BUNDLE ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs} ${QRC_Srcs} ${SHs})
target_link_libraries(${TARGETAPP_AVATreeExtensionGUI} ${STATIC_LIBS_COPIED} ModellerLib Qt5::Core Qt5::Gui Qt5::Widgets Qt5::OpenGL)
qt5_use_modules(${TARGETAPP_AVATreeExtensionGUI} Core Gui Widgets OpenGL)
endif()

####################################################################
####################################################################
####################################################################
####################################################################
####################################################################





if(MSVC)
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
	link_directories(${Boost_LIBRARY_DIR})
	link_directories(${Boost_LIBRARY_DIR})
	link_directories(${CMAKE_BINARY_DIR})
	link_directories(${CMAKE_SOURCE_DIR}/lib)



	include_directories(${PCL_INCLUDE_DIRS})
	link_directories(${PCL_LIBRARY_DIRS})
	add_definitions(${PCL_DEFINITIONS})


	set(MODELLER_LIBRARIES  ${CGAL_DIR}/auxiliary/gmp/lib/libmpfr-4.lib ${CGAL_DIR}/auxiliary/gmp/lib/libgmp-10.lib) 
	link_directories(${CMAKE_SOURCE_DIR}/lib)
	link_directories(${Boost_INCLUDE_DIR}/lib/win64)
	link_directories(${CGAL_BUILD_DIR}/lib)
	link_directories(${CGAL_DIR}/auxiliary/gmp/lib)



#Link Directories
if (DEFINED ENV{QTDIR})
  set(CMAKE_PREFIX_PATH $ENV{QTDIR})
else()
  set(CMAKE_PREFIX_PATH CACHE PATH "Location of Qt")
endif()
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
# The rest should just work for Qt5 sure!
# We have ui files, this will bring in the macro: qt5_wrap_ui
qt5_wrap_ui(UI_Srcs ${UIs})
# Wrap resources
qt5_add_resources(QRC_Srcs ${QRCs})
# Group resources
source_group("Resources" FILES
  ${UIs}
  ${QRCs}
  ${EXE_ICON} # Not present
)
source_group("Generated" FILES
  ${UI_Srcs}
  ${QRC_Srcs}
)
source_group("Shaders" FILES
  ${SHs}
)


add_executable(${TARGETAPP_AVATreeExtensionGUI} ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs} ${QRC_Srcs} ${SHs})
target_link_libraries(${TARGETAPP_AVATreeExtensionGUI} ${MODELLER_LIBRARIES})
target_link_libraries(${TARGETAPP_AVATreeExtensionGUI} Qt5::Core Qt5::Gui Qt5::Widgets Qt5::OpenGL)
#target_link_libraries(${TARGETAPP_AVATreeExtensionGUI} ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} ${PCL_VISUALIZATION_LIBRARIES} ${PCL_SURFACE_LIBRARIES} ${PCL_FILTERS_LIBRARIES} ${PCL_FEATURES_LIBRARIES} )
target_link_libraries(${TARGETAPP_AVATreeExtensionGUI} ModellerLib )
set_target_properties(${TARGETAPP_AVATreeExtensionGUI}        
        PROPERTIES        
        PROJECT_LABEL "${TARGETAPP_AVATreeExtensionGUI}"  
        FOLDER "UI"
		)
qt5_use_modules(${TARGETAPP_AVATreeExtensionGUI} Core Gui Widgets OpenGL)
endif()


if( ACTIVATE_CONCURRENCY )
  if( TBB_FOUND )
    include_directories(${TBB_INCLUDE_DIR})
    target_link_libraries(${TARGETAPP_AVATreeExtensionGUI} ${TBB_LIBRARY_RELEASE} ${TBB_MALLOC_LIBRARY_RELEASE} ${TBB_PROXY_LIBRARY_RELEASE})
  endif()
endif()
