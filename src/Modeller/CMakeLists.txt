cmake_minimum_required(VERSION 2.8.11)
add_subdirectory(Modeller)
if(MSVC)
add_subdirectory(ModellerExecutable)
#add_subdirectory(ModelFeatureGeneration)
endif (MSVC)
##########################################################################################


