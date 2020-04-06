cmake_minimum_required(VERSION 2.8.11)


# Activate concurrency?
option(ACTIVATE_CONCURRENCY "Enable concurrency" ON)

if( ACTIVATE_CONCURRENCY )
  find_package( TBB REQUIRED )
endif()


find_package(Boost)
find_package(CGAL)
find_package(Eigen3 REQUIRED)
find_package(PCL COMPONENTS common io visualization surface filters features)

find_package(Qt5Xml REQUIRED)
find_package(Qt5Network REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)

set(CGAL_BUILD_DIR CACHE PATH "") 

add_subdirectory(Modeller)
add_subdirectory(App)
add_subdirectory(AppV2) 
add_subdirectory(AVATreeExtensionGUI)  
#add_subdirectory(InitSim)
#add_subdirectory(EmbeddedWebServer)
add_subdirectory(WebServer)
##########################################################################################


