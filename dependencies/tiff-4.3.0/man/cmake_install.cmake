# Install script for directory: E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man

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
   "C:/Program Files (x86)/tiff/share/man/man1/fax2ps.1;C:/Program Files (x86)/tiff/share/man/man1/fax2tiff.1;C:/Program Files (x86)/tiff/share/man/man1/pal2rgb.1;C:/Program Files (x86)/tiff/share/man/man1/ppm2tiff.1;C:/Program Files (x86)/tiff/share/man/man1/raw2tiff.1;C:/Program Files (x86)/tiff/share/man/man1/tiff2bw.1;C:/Program Files (x86)/tiff/share/man/man1/tiff2pdf.1;C:/Program Files (x86)/tiff/share/man/man1/tiff2ps.1;C:/Program Files (x86)/tiff/share/man/man1/tiff2rgba.1;C:/Program Files (x86)/tiff/share/man/man1/tiffcmp.1;C:/Program Files (x86)/tiff/share/man/man1/tiffcp.1;C:/Program Files (x86)/tiff/share/man/man1/tiffcrop.1;C:/Program Files (x86)/tiff/share/man/man1/tiffdither.1;C:/Program Files (x86)/tiff/share/man/man1/tiffdump.1;C:/Program Files (x86)/tiff/share/man/man1/tiffgt.1;C:/Program Files (x86)/tiff/share/man/man1/tiffinfo.1;C:/Program Files (x86)/tiff/share/man/man1/tiffmedian.1;C:/Program Files (x86)/tiff/share/man/man1/tiffset.1;C:/Program Files (x86)/tiff/share/man/man1/tiffsplit.1")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "C:/Program Files (x86)/tiff/share/man/man1" TYPE FILE FILES
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/fax2ps.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/fax2tiff.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/pal2rgb.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/ppm2tiff.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/raw2tiff.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiff2bw.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiff2pdf.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiff2ps.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiff2rgba.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffcmp.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffcp.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffcrop.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffdither.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffdump.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffgt.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffinfo.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffmedian.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffset.1"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/tiffsplit.1"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/tiff/share/man/man3/libtiff.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFbuffer.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFClose.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFcodec.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFcolor.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFDataWidth.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFError.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFFieldDataType.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFFieldName.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFFieldPassCount.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFFieldReadCount.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFFieldTag.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFFieldWriteCount.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFFlush.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFGetField.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFmemory.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFOpen.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFPrintDirectory.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFquery.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadDirectory.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadEncodedStrip.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadEncodedTile.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadRawStrip.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadRawTile.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadRGBAImage.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadRGBAStrip.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadRGBATile.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadScanline.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFReadTile.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFRGBAImage.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFSetDirectory.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFSetField.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFsize.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFstrip.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFswab.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFtile.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFWarning.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFWriteDirectory.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFWriteEncodedStrip.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFWriteEncodedTile.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFWriteRawStrip.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFWriteRawTile.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFWriteScanline.3tiff;C:/Program Files (x86)/tiff/share/man/man3/TIFFWriteTile.3tiff")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "C:/Program Files (x86)/tiff/share/man/man3" TYPE FILE FILES
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/libtiff.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFbuffer.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFClose.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFcodec.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFcolor.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFDataWidth.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFError.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFFieldDataType.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFFieldName.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFFieldPassCount.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFFieldReadCount.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFFieldTag.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFFieldWriteCount.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFFlush.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFGetField.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFmemory.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFOpen.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFPrintDirectory.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFquery.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadDirectory.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadEncodedStrip.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadEncodedTile.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadRawStrip.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadRawTile.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadRGBAImage.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadRGBAStrip.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadRGBATile.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadScanline.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFReadTile.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFRGBAImage.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFSetDirectory.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFSetField.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFsize.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFstrip.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFswab.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFtile.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFWarning.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFWriteDirectory.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFWriteEncodedStrip.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFWriteEncodedTile.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFWriteRawStrip.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFWriteRawTile.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFWriteScanline.3tiff"
    "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/tiff-4.3.0/man/TIFFWriteTile.3tiff"
    )
endif()

