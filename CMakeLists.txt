# use cmake 2.8 or later
cmake_minimum_required (VERSION 2.8)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "LArPhysicsContent requires an out-of-source build.")
endif()

project(LArPhysicsContent)

# project version
# ATTN This package supports two build systems; please ensure version is specified here *and* in ups/product_deps
set(${PROJECT_NAME}_VERSION_MAJOR 00)
set(${PROJECT_NAME}_VERSION_MINOR 01)
set(${PROJECT_NAME}_VERSION_PATCH 00)
set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

#-------------------------------------------------------------------------------------------------------------------------------------------
# Dependencies
include(PandoraCMakeSettings)

# Prefer local include directory to any paths to installed header files
include_directories(./)

find_package(PandoraSDK 03.02.00 REQUIRED)
find_package(LArContent 03.11.00 REQUIRED)

foreach(package PandoraSDK LArContent)
    if(${package}_FOUND)
        include_directories(${${package}_INCLUDE_DIRS})
        link_libraries(${${package}_LIBRARIES})
        add_definitions (${${package}_DEFINITIONS})
    endif()
endforeach()

if(PANDORA_MONITORING)
    find_package(PandoraMonitoring 03.04.00 REQUIRED)
    include_directories(${PandoraMonitoring_INCLUDE_DIRS})
    link_libraries(${PandoraMonitoring_LIBRARIES})
    add_definitions(${PandoraMonitoring_DEFINITIONS})
    add_definitions("-DMONITORING")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

find_package(ROOT 5.26.00 REQUIRED COMPONENTS Eve Geom RGL EG TMVA)
include_directories(${ROOT_INCLUDE_DIRS})
link_libraries(${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES})
add_definitions(${ROOT_DEFINITIONS})
include(${ROOT_USE_FILE})

#-------------------------------------------------------------------------------------------------------------------------------------------
# Low level settings - compiler etc
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -pedantic -Wno-long-long -Wno-sign-compare -Wshadow -fno-strict-aliasing -std=c++11 ${CMAKE_CXX_FLAGS}")

include(CheckCXXCompilerFlag)
unset(COMPILER_SUPPORTS_CXX_FLAGS CACHE)
CHECK_CXX_COMPILER_FLAG(${CMAKE_CXX_FLAGS} COMPILER_SUPPORTS_CXX_FLAGS)

if(NOT COMPILER_SUPPORTS_CXX_FLAGS)
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support cxx flags ${CMAKE_CXX_FLAGS}")
endif()

#-------------------------------------------------------------------------------------------------------------------------------------------
# Build products

# - Generate any custom ROOT classes required for the ntuple
ROOT_GENERATE_DICTIONARY(PandoraNtupleClassesDict vector LINKDEF larphysicscontent/LinkDef.h MODULE PandoraNtupleClasses)

# - Collect sources - not ideal because you have to keep running CMake to pick up changes
file(GLOB_RECURSE LAR_PHYSICS_CONTENT_SRCS RELATIVE ${PROJECT_SOURCE_DIR} "larphysicscontent/*.cc" "test/*.cc")

# - Add library and properties
add_library(${PROJECT_NAME} SHARED ${LAR_PHYSICS_CONTENT_SRCS} PandoraNtupleClassesDict.cxx)
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})

add_library(PandoraNtupleClasses SHARED PandoraNtupleClassesDict.cxx)

#-------------------------------------------------------------------------------------------------------------------------------------------
# Install products

#Install the ROOT 6 PCM and rootmap files.
if (${ROOT_VERSION} VERSION_GREATER "6.0")
   install(
      FILES
      ${CMAKE_CURRENT_BINARY_DIR}/libPandoraNtupleClasses_rdict.pcm
      ${CMAKE_CURRENT_BINARY_DIR}/libPandoraNtupleClasses.rootmap
    ${CMAKE_CURRENT_BINARY_DIR}/lib/libPandoraNtupleClasses.so
      DESTINATION lib)
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

# - library
install(TARGETS ${PROJECT_NAME} DESTINATION lib COMPONENT Runtime)

# - headers
install(DIRECTORY ./larphysicscontent DESTINATION include COMPONENT Development FILES_MATCHING PATTERN "*.h")

# - support files
PANDORA_GENERATE_PACKAGE_CONFIGURATION_FILES(${PROJECT_NAME}Config.cmake ${PROJECT_NAME}ConfigVersion.cmake ${PROJECT_NAME}LibDeps.cmake)

#-------------------------------------------------------------------------------------------------------------------------------------------
# display some variables and write them to cache
PANDORA_DISPLAY_STD_VARIABLES()
