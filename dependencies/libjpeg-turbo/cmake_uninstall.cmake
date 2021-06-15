# This code is from the CMake FAQ

if (NOT EXISTS "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo/install_manifest.txt")
  message(FATAL_ERROR "Cannot find install manifest: \"E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo/install_manifest.txt\"")
endif(NOT EXISTS "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo/install_manifest.txt")

file(READ "E:/Frank/Desktop/Frank/CS184/3DWorld_git/dependencies/libjpeg-turbo/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")
list(REVERSE files)
foreach (file ${files})
  message(STATUS "Uninstalling \"$ENV{DESTDIR}${file}\"")
    if (EXISTS "$ENV{DESTDIR}${file}")
      execute_process(
        COMMAND "E:/Program Files (x86)/CMake/bin/cmake.exe" -E remove "$ENV{DESTDIR}${file}"
        OUTPUT_VARIABLE rm_out
        RESULT_VARIABLE rm_retval
      )
    if(NOT ${rm_retval} EQUAL 0)
      message(FATAL_ERROR "Problem when removing \"$ENV{DESTDIR}${file}\"")
    endif (NOT ${rm_retval} EQUAL 0)
  else (EXISTS "$ENV{DESTDIR}${file}")
    message(STATUS "File \"$ENV{DESTDIR}${file}\" does not exist.")
  endif (EXISTS "$ENV{DESTDIR}${file}")
endforeach(file)
