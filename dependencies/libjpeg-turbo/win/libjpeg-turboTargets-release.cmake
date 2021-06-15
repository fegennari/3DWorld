#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libjpeg-turbo::jpeg" for configuration "Release"
set_property(TARGET libjpeg-turbo::jpeg APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::jpeg PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/jpeg.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/jpeg62.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS libjpeg-turbo::jpeg )
list(APPEND _IMPORT_CHECK_FILES_FOR_libjpeg-turbo::jpeg "${_IMPORT_PREFIX}/lib/jpeg.lib" "${_IMPORT_PREFIX}/bin/jpeg62.dll" )

# Import target "libjpeg-turbo::turbojpeg" for configuration "Release"
set_property(TARGET libjpeg-turbo::turbojpeg APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::turbojpeg PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/turbojpeg.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/turbojpeg.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS libjpeg-turbo::turbojpeg )
list(APPEND _IMPORT_CHECK_FILES_FOR_libjpeg-turbo::turbojpeg "${_IMPORT_PREFIX}/lib/turbojpeg.lib" "${_IMPORT_PREFIX}/bin/turbojpeg.dll" )

# Import target "libjpeg-turbo::turbojpeg-static" for configuration "Release"
set_property(TARGET libjpeg-turbo::turbojpeg-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::turbojpeg-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "ASM_NASM;C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/turbojpeg-static.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS libjpeg-turbo::turbojpeg-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_libjpeg-turbo::turbojpeg-static "${_IMPORT_PREFIX}/lib/turbojpeg-static.lib" )

# Import target "libjpeg-turbo::jpeg-static" for configuration "Release"
set_property(TARGET libjpeg-turbo::jpeg-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::jpeg-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "ASM_NASM;C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/jpeg-static.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS libjpeg-turbo::jpeg-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_libjpeg-turbo::jpeg-static "${_IMPORT_PREFIX}/lib/jpeg-static.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
