# Minimum required version of CMake
cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR} CACHE TYPE PATH)

function(enable_unity_build UB_SUFFIX SOURCE_VARIABLE_NAME)
  set(files ${SOURCE_VARIABLE_NAME})
  # Generate a unique filename for the unity build translation unit
  set(unit_build_file ${CMAKE_CURRENT_BINARY_DIR}/ub_${UB_SUFFIX}.cpp)
  # Exclude all translation units from compilation
  set_source_files_properties(${files} PROPERTIES HEADER_FILE_ONLY true)
  # Open the ub file
  FILE(WRITE ${unit_build_file} "// Unity Build generated by CMake\n")
  # Add include statement for each translation unit
  foreach(source_file ${files} )
    FILE( APPEND ${unit_build_file} "#include <${CMAKE_CURRENT_SOURCE_DIR}/${source_file}>\n")
  endforeach(source_file)
  # Complement list of translation units with the name of ub
  set(${SOURCE_VARIABLE_NAME} ${${SOURCE_VARIABLE_NAME}} ${unit_build_file} PARENT_SCOPE)  
endfunction(enable_unity_build)





set(TARGETMODEXE ModellerCMD CACHE TYPE STRING)
project(${TARGETMODEXE})
file(GLOB MODELLEREXEC_SRC
    "../Modeller/include/*.h"
    "../Modeller/*.cpp"
	"*.cpp"
	"../PoissonRecon/*.h"
    "../PoissonRecon/*.cpp"
	"../PoissonRecon/PoissonRecon/*.cpp"
)


FOREACH(MODELLEREXEC_SRC ${MODELLEREXEC_SRC})
 MESSAGE(STATUS "Process file: ${MODELLEREXEC_SRC}")
ENDFOREACH(MODELLEREXEC_SRC)

if (MSVC)

	include(${ITK_USE_FILE})

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

	set(MODELLER_LIBRARIES  ${CGAL_DIR}/auxiliary/gmp/lib/libmpfr-4.lib ${ITK_LIBRARIES} ${CGAL_DIR}/auxiliary/gmp/lib/libgmp-10.lib) 
	link_directories(${CMAKE_SOURCE_DIR}/lib)
	link_directories(${Boost_INCLUDE_DIR}/lib/win64)
	link_directories(${CGAL_BUILD_DIR}/lib)
	link_directories(${CGAL_DIR}/auxiliary/gmp/lib)


	# Comment to disable unit build
	# enable_unity_build(UnityBuildDemo ${MODELLEREXEC_SRC})
 
	add_executable(${TARGETMODEXE} ${MODELLEREXEC_SRC} )
	set_target_properties(${TARGETMODEXE}        
        PROPERTIES        
        PROJECT_LABEL "${TARGETMODEXE}"  
        FOLDER "CMD"
		)
	
	
	if( ACTIVATE_CONCURRENCY )
		#find_package( TBB REQUIRED )
		if( TBB_FOUND )
			#include(${TBB_USE_FILE})
			target_link_libraries(${TARGETMODEXE} ${TBB_LIBRARIES})
		endif()
	endif()


	target_link_libraries(${TARGETMODEXE} ${MODELLER_LIBRARIES} )
	target_link_libraries(${TARGETMODEXE} ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} ${PCL_VISUALIZATION_LIBRARIES} ${PCL_SURFACE_LIBRARIES} ${PCL_FILTERS_LIBRARIES} ${PCL_FEATURES_LIBRARIES})
	
endif()
















