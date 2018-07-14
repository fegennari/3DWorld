cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_policy(VERSION 2.6)

project(gli)
enable_testing()

add_definitions(-D_CRT_SECURE_NO_WARNINGS)

include_directories(".")

if(MSVC)
	add_definitions(/Za)
endif()

################################
# Add package

macro(addExternalPackageGTC NAME DIRECTORY)
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/external/${DIRECTORY}/util)
	set(${NAME}_PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/${DIRECTORY})

	find_package(${NAME})
	if(NOT ${NAME}_FOUND)
		message("ERROR: ${NAME} not found")
	endif(NOT ${NAME}_FOUND)
endmacro(addExternalPackageGTC)

################################
# Add custom target to copy all data

set(TARGET_DATA_COPY DATA_COPY)
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
	add_custom_target(
		${TARGET_DATA_COPY}
		COMMAND ${CMAKE_COMMAND} -E echo "In source build")
else(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
	add_custom_target(
		${TARGET_DATA_COPY}
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data)
endif(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})

set(OGL_DATA_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)

################################
# Add GLM 

set(GLM_DIRECTORY glm-0.9.5.1)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/${GLM_DIRECTORY})

################################
# Add subdirectory

add_subdirectory(gli)
add_subdirectory(test)
#add_subdirectory(doc)
