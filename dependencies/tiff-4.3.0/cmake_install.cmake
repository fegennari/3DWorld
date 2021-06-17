# Install script for directory: E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/tiff")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/tiff/lib/pkgconfig/libtiff-4.pc")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "C:/Program Files (x86)/tiff/lib/pkgconfig" TYPE FILE FILES "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/libtiff-4.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/port/cmake_install.cmake")
  include("E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/libtiff/cmake_install.cmake")
  include("E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/tools/cmake_install.cmake")
  include("E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/test/cmake_install.cmake")
  include("E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/contrib/cmake_install.cmake")
  include("E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/build/cmake_install.cmake")
  include("E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/cmake_install.cmake")
  include("E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/html/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
