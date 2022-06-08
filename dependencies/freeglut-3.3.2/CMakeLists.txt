CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0 FATAL_ERROR)
PROJECT(freeglut)

# for multiarch LIBDIR support (requires cmake>=2.8.8)
INCLUDE(GNUInstallDirs)

# NOTE: On Windows and Cygwin, the dll's are placed in the
# CMAKE_RUNTIME_OUTPUT_DIRECTORY, while their corresponding import
# libraries end up in CMAKE_ARCHIVE_OUTPUT_DIRECTORY. On other
# platforms, such as Linux, the shared libraries are put in
# CMAKE_ARCHIVE_OUTPUT_DIRECTORY instead.
# Static libraries end up in CMAKE_ARCHIVE_OUTPUT_DIRECTORY on all
# platforms.
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# setup version numbers
# XXX: Update these for each release!
set(VERSION_MAJOR 3)
set(VERSION_MINOR 2)
set(VERSION_PATCH 2)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# Update fg_version.h to match the versions number here in cmake
CONFIGURE_FILE(src/fg_version.h.in src/fg_version.h)

# shared lib version numbers (XXX: change before release)
set(SO_MAJOR 3)		# increment on backwards incompatible API/ABI changes
set(SO_MINOR 11)	# increment on backwards compatible API additions
set(SO_REV 1)		# increment if only internal changes happened between releases

# FREEGLUT_BUILD_SHARED_LIBS is already a standard CMake variable, but we need to
# re-declare it here so it will show up in the GUI.
# by default, we want to build both
OPTION(FREEGLUT_BUILD_SHARED_LIBS "Build FreeGLUT shared library." ON)
OPTION(FREEGLUT_BUILD_STATIC_LIBS "Build FreeGLUT static library." ON)

# option for whether warnings and errors should be printed
OPTION(FREEGLUT_PRINT_ERRORS "Lib prints errors to stderr" ON)
#MARK_AS_ADVANCED(FREEGLUT_PRINT_ERRORS)
OPTION(FREEGLUT_PRINT_WARNINGS "Lib prints warnings to stderr" ON)
#MARK_AS_ADVANCED(FREEGLUT_PRINT_WARNINGS)

# option to also copy .pdb files to install directory when executing
# INSTALL target
IF(MSVC)
    OPTION(INSTALL_PDB "Also install .pdb files" ON)
ELSE()
    SET(INSTALL_PDB OFF)
ENDIF()

# OpenGL ES support
IF(ANDROID)
    OPTION(FREEGLUT_GLES "Use OpenGL ES (requires EGL)" ON)
ELSE()
    OPTION(FREEGLUT_GLES "Use OpenGL ES (requires EGL)" OFF)
ENDIF()

# option to build either as "glut" (ON) or "freeglut" (OFF)
IF(WIN32)
    OPTION(FREEGLUT_REPLACE_GLUT "Be a replacement for GLUT" OFF)
ELSE()
    OPTION(FREEGLUT_REPLACE_GLUT "Be a replacement for GLUT" ON)
ENDIF()

IF(NOT WIN32)
    # Wayland support
    OPTION(FREEGLUT_WAYLAND "Use Wayland (no X11)" OFF)
ENDIF()


SET(FREEGLUT_HEADERS
    include/GL/freeglut.h
    include/GL/freeglut_ucall.h
    include/GL/freeglut_ext.h
    include/GL/freeglut_std.h
)
IF(FREEGLUT_REPLACE_GLUT)
    LIST(APPEND FREEGLUT_HEADERS
        include/GL/glut.h
    )
ENDIF()
SET(FREEGLUT_SRCS
    ${FREEGLUT_HEADERS}
    src/fg_callbacks.c
    src/fg_cursor.c
    src/fg_display.c
    src/fg_ext.c
    src/fg_font_data.c
    src/fg_gamemode.c
    src/fg_geometry.c
    src/fg_gl2.c
    src/fg_gl2.h
    src/fg_init.c
    src/fg_init.h
    src/fg_internal.h
    src/fg_callback_macros.h
    src/fg_input_devices.c
    src/fg_joystick.c
    src/fg_main.c
    src/fg_misc.c
    src/fg_overlay.c
    src/fg_spaceball.c
    src/fg_state.c
    src/fg_stroke_mono_roman.c
    src/fg_stroke_roman.c
    src/fg_structure.c
    src/fg_teapot.c
    src/fg_teapot_data.h
    src/fg_videoresize.c
    src/fg_window.c
)
# TODO: OpenGL ES requires a compatible version of these files:
IF(NOT FREEGLUT_GLES)
    LIST(APPEND FREEGLUT_SRCS
        src/fg_font.c
        src/fg_menu.c
    )
ELSE()
    LIST(APPEND FREEGLUT_SRCS
        src/gles_stubs.c
    )
ENDIF()

IF(WIN32)
    LIST(APPEND FREEGLUT_SRCS
        src/mswin/fg_cursor_mswin.c
        src/mswin/fg_display_mswin.c
        src/mswin/fg_ext_mswin.c
        src/mswin/fg_gamemode_mswin.c
        src/mswin/fg_init_mswin.c
        src/mswin/fg_internal_mswin.h
        src/mswin/fg_input_devices_mswin.c
        src/mswin/fg_joystick_mswin.c
        src/mswin/fg_main_mswin.c
        src/mswin/fg_menu_mswin.c
        src/mswin/fg_spaceball_mswin.c
        src/mswin/fg_state_mswin.c
        src/mswin/fg_structure_mswin.c
        src/mswin/fg_window_mswin.c
        ${CMAKE_BINARY_DIR}/freeglut.rc # generated below from freeglut.rc.in
    )
    IF (MSVC AND NOT CMAKE_CL_64)
        # .def file only for 32bit Windows builds (TODO: MSVC only right
        # now, needed for any other Windows platform?)
        LIST(APPEND FREEGLUT_SRCS
            ${CMAKE_BINARY_DIR}/freeglutdll.def # generated below from src/freeglutdll.def.in
        )
    ENDIF()

ELSEIF(ANDROID OR BLACKBERRY)
    # BlackBerry and Android share some similar design concepts and ideas, as with many mobile devices.
    # As such, some classes can be shared between the two. XXX: Possibly rename shareable classes to
    # a more generic name. *_stub? *_mobile?
    LIST(APPEND FREEGLUT_SRCS
        src/android/fg_cursor_android.c
        src/android/fg_ext_android.c
        src/android/fg_gamemode_android.c
        src/android/fg_joystick_android.c
        src/android/fg_spaceball_android.c
    )
    IF(ANDROID)
        LIST(APPEND FREEGLUT_SRCS
            src/android/native_app_glue/android_native_app_glue.c
            src/android/native_app_glue/android_native_app_glue.h
            src/android/fg_internal_android.h
            src/android/fg_init_android.c
            src/android/fg_input_devices_android.c
            src/android/fg_main_android.c
            src/android/fg_main_android.h
            src/android/fg_runtime_android.c
            src/android/fg_state_android.c
            src/android/fg_structure_android.c
            src/android/fg_window_android.c
        )
    ELSE()
        LIST(APPEND FREEGLUT_SRCS
            src/blackberry/fg_internal_blackberry.h
            src/blackberry/fg_init_blackberry.c
            src/x11/fg_input_devices_x11.c
            src/blackberry/fg_main_blackberry.c
            src/blackberry/fg_state_blackberry.c
            src/blackberry/fg_structure_blackberry.c
            src/blackberry/fg_window_blackberry.c
        )
    ENDIF()

ELSE()
    # UNIX (Wayland)
    IF(FREEGLUT_WAYLAND)
        LIST(APPEND FREEGLUT_SRCS
            src/wayland/fg_cursor_wl.c
            src/wayland/fg_ext_wl.c
            src/wayland/fg_gamemode_wl.c
            src/wayland/fg_init_wl.c
            src/wayland/fg_internal_wl.h
            src/wayland/fg_input_devices_wl.c
            src/wayland/fg_main_wl.c
            src/wayland/fg_state_wl.c
            src/wayland/fg_structure_wl.c
            src/wayland/fg_window_wl.c
            # font, serial port & joystick code are agnostic
            src/x11/fg_glutfont_definitions_x11.c
            src/x11/fg_input_devices_x11.c
            src/x11/fg_joystick_x11.c
        )
    # UNIX (X11)
    ELSE()
        LIST(APPEND FREEGLUT_SRCS
            src/x11/fg_cursor_x11.c
            src/x11/fg_ext_x11.c
            src/x11/fg_gamemode_x11.c
            src/x11/fg_glutfont_definitions_x11.c
            src/x11/fg_init_x11.c
            src/x11/fg_internal_x11.h
            src/x11/fg_input_devices_x11.c
            src/x11/fg_joystick_x11.c
            src/x11/fg_main_x11.c
            src/x11/fg_menu_x11.c
            src/x11/fg_spaceball_x11.c
            src/x11/fg_state_x11.c
            src/x11/fg_structure_x11.c
            src/x11/fg_window_x11.c
            src/x11/fg_xinput_x11.c
        )
        IF(NOT(FREEGLUT_GLES))
            LIST(APPEND FREEGLUT_SRCS
                src/x11/fg_internal_x11_glx.h
                src/x11/fg_display_x11_glx.c
                src/x11/fg_state_x11_glx.c
                src/x11/fg_state_x11_glx.h
                src/x11/fg_window_x11_glx.c
                src/x11/fg_window_x11_glx.h
            )
        ENDIF()
    ENDIF()
ENDIF()

# OpenGL ES requires EGL, and so does Wayland
IF(FREEGLUT_GLES OR FREEGLUT_WAYLAND)
    LIST(APPEND FREEGLUT_SRCS
        src/egl/fg_internal_egl.h
        src/egl/fg_display_egl.c
        src/egl/fg_ext_egl.c
        src/egl/fg_init_egl.c
        src/egl/fg_init_egl.h
        src/egl/fg_state_egl.c
        src/egl/fg_state_egl.h
        src/egl/fg_structure_egl.c
        src/egl/fg_structure_egl.h
        src/egl/fg_window_egl.c
        src/egl/fg_window_egl.h
    )
ENDIF()

# For OpenGL ES (GLES): compile with -DFREEGLUT_GLES to cleanly
# bootstrap headers inclusion in freeglut_std.h; this constant also
# need to be defined in client applications (e.g. through pkg-config),
# but do use GLES constants directly for all other needs
# GLES1 and GLES2 libraries are compatible and can be co-linked.
IF(FREEGLUT_GLES)
    LIST(APPEND PUBLIC_DEFINITIONS -DFREEGLUT_GLES)
    LIST(APPEND LIBS GLESv2 GLESv1_CM EGL)
ELSE()
  FIND_PACKAGE(OpenGL REQUIRED)
  LIST(APPEND LIBS ${OPENGL_gl_LIBRARY})
  INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
ENDIF()

# For Wayland: compile with -DFREEGLUT_WAYLAND and pull EGL
IF(FREEGLUT_WAYLAND)
  ADD_DEFINITIONS(-DFREEGLUT_WAYLAND)
  LIST(APPEND LIBS wayland-client wayland-cursor wayland-egl EGL xkbcommon)
ENDIF()

# lib m for math, not needed on windows
IF (NOT WIN32)
    # For compilation:
    LIST(APPEND LIBS m)
    # For CHECK_FUNCTION_EXISTS:
    LIST(APPEND CMAKE_REQUIRED_LIBRARIES m)
ENDIF()

IF(WIN32)
    # hide insecure CRT warnings, common practice
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    IF(MSVC)
        SET( CMAKE_DEBUG_POSTFIX "d" )
    ENDIF(MSVC)
    
    IF(NOT(MSVC_VERSION LESS "1600"))
    # minimum requirement for WM_TOUCH device
        ADD_DEFINITIONS(-D_WIN32_WINNT=0x0601)
        ADD_DEFINITIONS(-DWINVER=0x0601)
    ELSEIF(NOT(MSVC_VERSION LESS "1300"))
    # minimum requirement for spaceball device
        ADD_DEFINITIONS(-D_WIN32_WINNT=0x0501)
        ADD_DEFINITIONS(-DWINVER=0x0501)
    ELSE()
    # enable the use of Win2000 APIs (needed for really old compilers like MSVC6)
        ADD_DEFINITIONS(-D_WIN32_WINNT=0x0500)
        ADD_DEFINITIONS(-DWINVER=0x0500)
    ENDIF()

ENDIF()

IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  IF(NOT(ANDROID OR BLACKBERRY OR FREEGLUT_WAYLAND))
    # not setting -ansi as EGL/KHR headers doesn't support it
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic")
  ENDIF()
ENDIF(CMAKE_COMPILER_IS_GNUCC)

INCLUDE(CheckIncludeFiles)
IF(UNIX AND NOT(ANDROID OR BLACKBERRY OR FREEGLUT_WAYLAND))
    FIND_PACKAGE(X11 REQUIRED)
    INCLUDE_DIRECTORIES(${X11_X11_INCLUDE_PATH})
    LIST(APPEND LIBS ${X11_X11_LIB})
    IF(X11_Xrandr_FOUND)
        SET(HAVE_X11_EXTENSIONS_XRANDR_H TRUE)
        LIST(APPEND LIBS ${X11_Xrandr_LIB})
    ENDIF()
    IF(X11_xf86vmode_FOUND)
        SET(HAVE_X11_EXTENSIONS_XF86VMODE_H TRUE)
        LIST(APPEND LIBS ${X11_Xxf86vm_LIB})
    ENDIF()
    IF(X11_Xinput_FOUND)
        # Needed for multi-touch:
        CHECK_INCLUDE_FILES("${X11_Xinput_INCLUDE_PATH}/X11/extensions/XInput2.h" HAVE_X11_EXTENSIONS_XINPUT2_H)
        LIST(APPEND LIBS ${X11_Xinput_LIB})
    ELSE()
        MESSAGE(FATAL_ERROR "Missing X11's XInput.h (X11/extensions/XInput.h)")
    ENDIF()
ENDIF()
IF(ANDROID)
    # -landroid for ANativeWindow
    # -llog for native Android logging
    LIST(APPEND LIBS android log)
ELSEIF(BLACKBERRY)
    # -lbps for event loop
    # -screen for native screen
    LIST(APPEND LIBS bps screen)

    if(NOT PLAYBOOK)
        # -lslog2 for logging
        # -pps for low-level screen manipulation
        LIST(APPEND LIBS slog2 pps)
    ENDIF()
ENDIF()

INCLUDE(CheckFunctionExists)
INCLUDE(CheckTypeSize)
INCLUDE(CheckCCompilerFlag)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(unistd.h 	HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(sys/time.h 	HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES(stdbool.h 	HAVE_STDBOOL_H)
CHECK_INCLUDE_FILES(sys/param.h HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILES(sys/ioctl.h HAVE_SYS_IOCTL_H)
CHECK_INCLUDE_FILES(fcntl.h 	HAVE_FCNTL_H)
CHECK_INCLUDE_FILES(usbhid.h 	HAVE_USBHID_H)
CHECK_FUNCTION_EXISTS(gettimeofday HAVE_GETTIMEOFDAY)
CHECK_FUNCTION_EXISTS(XParseGeometry   HAVE_XPARSEGEOMETRY)
IF (NOT HAVE_XPARSEGEOMETRY)
   LIST(APPEND FREEGLUT_SRCS
        src/util/xparsegeometry_repl.c
        src/util/xparsegeometry_repl.h)
   SET(NEED_XPARSEGEOMETRY_IMPL TRUE)
ENDIF()
# decide on suitable type for internal time keeping, 64-bit if possible
CHECK_INCLUDE_FILES(stdint.h    HAVE_STDINT_H)
CHECK_INCLUDE_FILES(inttypes.h  HAVE_INTTYPES_H)
IF (NOT (HAVE_STDINT_H OR HAVE_INTTYPES_H))
    IF (MSVC)
        # Some old Microsoft VC releases don't support unsigned long
        # long, but all we care about is support for unsigned __int64 on
        # MSVC, so test for presence of that type
        CHECK_TYPE_SIZE("unsigned __int64" U__INT64 BUILTIN_TYPES_ONLY)
    ELSEIF()
        CHECK_TYPE_SIZE("unsigned long long" ULONG_LONG BUILTIN_TYPES_ONLY)
    ENDIF()
ENDIF()

# The generated config.h is placed in the project's build directory, just to
# ensure that all CMake-generated files are kept away from the main source tree.
# As a result, the build directory must to be added to the include path list.
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_BINARY_DIR}/config.h)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
ADD_DEFINITIONS(-DHAVE_CONFIG_H)
IF(WIN32)
    # we also have to generate freeglut.rc, which contains the version
    # number
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/freeglut.rc.in ${CMAKE_BINARY_DIR}/freeglut.rc)
    IF (MSVC AND NOT CMAKE_CL_64)
        # .def file only for 32bit Windows builds with Visual Studio
        CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/freeglutdll.def.in ${CMAKE_BINARY_DIR}/freeglutdll.def)
    ENDIF()
ENDIF()

IF(FREEGLUT_BUILD_SHARED_LIBS)
    ADD_LIBRARY(freeglut SHARED ${FREEGLUT_SRCS})
ENDIF()
IF(FREEGLUT_BUILD_STATIC_LIBS)
    ADD_LIBRARY(freeglut_static STATIC ${FREEGLUT_SRCS})
ENDIF()


SET(LIBNAME freeglut)
IF(WIN32)
    IF(FREEGLUT_REPLACE_GLUT)
        SET(LIBNAME glut)
    ENDIF()

    LIST(APPEND LIBS winmm)
    IF(FREEGLUT_BUILD_SHARED_LIBS)
        TARGET_COMPILE_DEFINITIONS(freeglut PRIVATE FREEGLUT_EXPORTS)
        SET_TARGET_PROPERTIES(freeglut PROPERTIES OUTPUT_NAME ${LIBNAME})
    ENDIF()
    IF(FREEGLUT_BUILD_STATIC_LIBS)
        TARGET_COMPILE_DEFINITIONS(freeglut_static PUBLIC FREEGLUT_STATIC)
        IF(FREEGLUT_REPLACE_GLUT)
            SET_TARGET_PROPERTIES(freeglut_static PROPERTIES OUTPUT_NAME ${LIBNAME})
        ENDIF()
        # need to set machine:x64 for linker, at least for VC10, and
        # doesn't hurt for older compilers:
        # http://public.kitware.com/Bug/view.php?id=11240#c22768
        IF (CMAKE_CL_64)
            SET_TARGET_PROPERTIES(freeglut_static PROPERTIES STATIC_LIBRARY_FLAGS "/machine:x64")
        ENDIF()
    ENDIF()
ELSE()
    # on UNIX we need to make sure:
    # - all shared libraries must have a soname/version, see :
    #   http://sourceware.org/autobook/autobook/autobook_91.html#SEC91
    #   http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
    #   Current: -version-info 12:0:9 -> 3.9.0
    #   Note: most platforms now prefer the latter major.minor.revision form
    #   (e.g. FreeBSD, cf. http://debbugs.gnu.org/cgi/bugreport.cgi?bug=8765),
    #   or special-cased FreeGLUT long ago (e.g. .so.4 on OpenBSD), so
    #   the lack of support for current:revision:age in CMake should
    #   not be a problem.
    # - the output library should be named glut so it'll be linkable with -lglut
    #   (unless FREEGLUT_REPLACE_GLUT is false).
    # - the shared library should link to the dependency libraries so that the user
    #   won't have to link them explicitly (they shouldn't have to know that we depend
    #   on Xrandr or Xxf86vm)
    IF(FREEGLUT_GLES)
        SET(LIBNAME freeglut-gles)
    ELSE()
        IF(FREEGLUT_REPLACE_GLUT)
            SET(LIBNAME glut)
        ENDIF()
    ENDIF()

    IF(FREEGLUT_BUILD_SHARED_LIBS)
      SET_TARGET_PROPERTIES(freeglut PROPERTIES VERSION ${SO_MAJOR}.${SO_MINOR}.${SO_REV} SOVERSION ${SO_MAJOR} OUTPUT_NAME ${LIBNAME})
    ENDIF()
    IF(FREEGLUT_BUILD_STATIC_LIBS)
      SET_TARGET_PROPERTIES(freeglut_static PROPERTIES OUTPUT_NAME ${LIBNAME})
    ENDIF()
    IF(ANDROID)
        # Not in CMake toolchain file, because the toolchain
        # file is called several times and generally doesn't
        # seem to be meant for modifying CFLAGS:
        # '-mandroid' is not mandatory but doesn't hurt
        # '-O0 -gstabs+' helps the currently buggy GDB port
        # Too late to manipulate ENV: SET(ENV{CFLAGS} "$ENV{CFLAGS} -mandroid")
        # Not using _INIT variables, they seem to be used internally only
        IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mandroid")
        ENDIF()
        CHECK_C_COMPILER_FLAG(-O0 HAVE_O0_FLAG)
        IF(HAVE_O0_FLAG)
            SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
        ENDIF()
        CHECK_C_COMPILER_FLAG(-gstabs+ HAVE_GSTABSP_FLAG)
        IF(HAVE_gstabsp_FLAG)
            SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -gstabs+")
        ENDIF()
    ENDIF()
ENDIF()
IF(FREEGLUT_BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(freeglut ${LIBS})
    TARGET_COMPILE_DEFINITIONS(freeglut PUBLIC ${PUBLIC_DEFINITIONS})
ENDIF()
IF(FREEGLUT_BUILD_STATIC_LIBS)
    TARGET_LINK_LIBRARIES(freeglut_static ${LIBS})
    TARGET_COMPILE_DEFINITIONS(freeglut_static PUBLIC ${PUBLIC_DEFINITIONS})
ENDIF()

IF(FREEGLUT_BUILD_SHARED_LIBS)
    INSTALL(TARGETS freeglut EXPORT FreeGLUTTargets
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION include
    )
    IF(INSTALL_PDB)
        IF(CMAKE_GENERATOR MATCHES "^Visual Studio")
            INSTALL(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/freeglut${CMAKE_DEBUG_POSTFIX}.pdb
                DESTINATION bin
                CONFIGURATIONS Debug
                COMPONENT Devel)
        ELSE()
            INSTALL(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/freeglut${CMAKE_DEBUG_POSTFIX}.pdb
                DESTINATION bin
                CONFIGURATIONS Debug
                COMPONENT Devel)
        ENDIF()
    ENDIF()
ENDIF()
IF(FREEGLUT_BUILD_STATIC_LIBS)
    INSTALL(TARGETS freeglut_static EXPORT FreeGLUTTargets
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION include
    )
ENDIF()
INSTALL(FILES ${FREEGLUT_HEADERS} DESTINATION include/GL COMPONENT Devel)



# Optionally build demos, on by default.
option( FREEGLUT_BUILD_DEMOS "Build FreeGLUT demos." ON )

SET(DEMO_LIBS ${OPENGL_glu_LIBRARY} ${LIBS})
# lib m for math, not needed on windows
IF (NOT WIN32)
    LIST(APPEND DEMO_LIBS m)
ENDIF()

MACRO(ADD_DEMO name)
    IF( FREEGLUT_BUILD_DEMOS )
        IF(FREEGLUT_BUILD_SHARED_LIBS)
            ADD_EXECUTABLE(${name} ${ARGN})
            TARGET_LINK_LIBRARIES(${name} ${DEMO_LIBS} freeglut)
            IF(WIN32 AND MSVC)
                SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX d)
            ENDIF()
        ENDIF()
        IF(FREEGLUT_BUILD_STATIC_LIBS)
            ADD_EXECUTABLE(${name}_static ${ARGN})
            TARGET_LINK_LIBRARIES(${name}_static ${DEMO_LIBS} freeglut_static)
            IF(WIN32 AND MSVC)
                SET_TARGET_PROPERTIES(${name}_static PROPERTIES DEBUG_POSTFIX d)
            ENDIF()
        ENDIF()
    ENDIF()
ENDMACRO()

ADD_DEMO(CallbackMaker   progs/demos/CallbackMaker/CallbackMaker.c)
ADD_DEMO(Fractals        progs/demos/Fractals/fractals.c)
ADD_DEMO(Fractals_random progs/demos/Fractals_random/fractals_random.c)
ADD_DEMO(Lorenz          progs/demos/Lorenz/lorenz.c)
IF (NOT WIN32)
    ADD_DEMO(One             progs/demos/One/one.c)
ELSE()
    ADD_DEMO(One             progs/demos/One/one.c
                             progs/demos/One/one.rc)
ENDIF()
ADD_DEMO(Resizer         progs/demos/Resizer/Resizer.cpp)
ADD_DEMO(multi-touch     progs/demos/multi-touch/multi-touch.c)
ADD_DEMO(shapes          progs/demos/shapes/shapes.c
                         progs/demos/shapes/glmatrix.h
                         progs/demos/shapes/glmatrix.c)
ADD_DEMO(smooth_opengl3  progs/demos/smooth_opengl3/smooth_opengl3.c)
IF(UNIX)
    ADD_DEMO(spaceball       progs/demos/spaceball/spaceball.c
                             progs/demos/spaceball/vmath.c
                             progs/demos/spaceball/vmath.h)
ENDIF()
ADD_DEMO(subwin          progs/demos/subwin/subwin.c)
ADD_DEMO(timer           progs/demos/timer/timer.c)
ADD_DEMO(timer_callback  progs/demos/timer_callback/timer.c)



# pkg-config support, to install at $(libdir)/pkgconfig
# Define static build dependencies
IF(WIN32)
  SET(PC_LIBS_PRIVATE "-lopengl32 -lwinmm -lgdi32 -lm")
ELSEIF(FREEGLUT_GLES)
  IF(ANDROID)
    SET(PC_LIBS_PRIVATE "-llog -landroid -lGLESv2 -lGLESv1_CM -lEGL -lm")
  ELSEIF(BLACKBERRY)
    IF(PLAYBOOK)
      SET(PC_LIBS_PRIVATE "-lbps -lscreen -lGLESv2 -lGLESv1_CM -lEGL -lm")
    ELSE()
      SET(PC_LIBS_PRIVATE "-lbps -lslog2 -lscreen -lGLESv2 -lGLESv1_CM -lEGL -lm")
    ENDIF()
  ELSEIF(FREEGLUT_WAYLAND)
    SET(PC_LIBS_PRIVATE "-lwayland-client -lwayland-cursor -lwayland-egl -lGLESv2 -lGLESv1_CM -lEGL -lxkbcommon -lm")
  ELSE()
    SET(PC_LIBS_PRIVATE "-lX11 -lXxf86vm -lXrandr -lGLESv2 -lGLESv1_CM -lEGL -lm")
  ENDIF()
ELSE()
  IF(FREEGLUT_WAYLAND)
    SET(PC_LIBS_PRIVATE "-lwayland-client -lwayland-cursor -lwayland-egl -lGL -lxkbcommon -lm")
  ELSE()
    SET(PC_LIBS_PRIVATE "-lX11 -lXxf86vm -lXrandr -lGL -lm")
  ENDIF()
ENDIF()
# Client applications need to define FreeGLUT GLES version to
# bootstrap headers inclusion in freeglut_std.h:
SET(PC_LIBNAME ${LIBNAME})
SET(PC_FILENAME "${LIBNAME}.pc")
IF(FREEGLUT_GLES)
  SET(PC_CFLAGS "-DFREEGLUT_GLES")
ENDIF()
IF(FREEGLUT_BUILD_STATIC_LIBS)
  LIST(APPEND PC_CFLAGS -DFREEGLUT_STATIC)
ENDIF()
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/freeglut.pc.in ${CMAKE_BINARY_DIR}/freeglut.pc @ONLY)
INSTALL(FILES ${CMAKE_BINARY_DIR}/freeglut.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/ RENAME ${PC_FILENAME} COMPONENT Devel)
# TODO: change the library and .pc name when building for GLES,
# e.g. -lglut-GLES

INCLUDE(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(
    "${CMAKE_CURRENT_BINARY_DIR}/FreeGLUT/FreeGLUTConfigVersion.cmake"
    VERSION ${VERSION}
    COMPATIBILITY AnyNewerVersion
)

# needs cmake 3.0 (as does the "INCLUDES DIRECTORY" argument to install(TARGETS)):
EXPORT(EXPORT FreeGLUTTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/FreeGLUT/FreeGLUTTargets.cmake"
    NAMESPACE FreeGLUT::
)
CONFIGURE_FILE(FreeGLUTConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/FreeGLUT/FreeGLUTConfig.cmake"
    @ONLY
)
SET(ConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/FreeGLUT)
INSTALL(EXPORT FreeGLUTTargets
    FILE FreeGLUTTargets.cmake
    NAMESPACE FreeGLUT::
    DESTINATION ${ConfigPackageLocation}
)
INSTALL(
    FILES
    "${CMAKE_CURRENT_BINARY_DIR}/FreeGLUT/FreeGLUTConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/FreeGLUT/FreeGLUTConfigVersion.cmake"
    DESTINATION ${ConfigPackageLocation}
    COMPONENT Devel
)
