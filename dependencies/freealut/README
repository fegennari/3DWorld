The freealut library
=====================

freealut is a free implementation of OpenAL's ALUT standard. See the file
AUTHORS for the people involved.

This is a slightly-updated version maintained at
https://github.com/vance-group/freealut that contains an updated
CMake build system, and potentially other fixes over time.


Building freealut
=================

While the "official" build system of upstream may be autotools,
that build is sufficiently outdated to require patching on modern
Linux operating systems. The CMake build system has been updated
to work smoothly on all platforms (tested on Linux and Windows,
including cross-building), so I'd recommend using it.

(Instructions for older build systems are still included for completeness.)


This document explains briefly how to build with CMake via out-
of-tree build:

- Change to the alut dir.
- Create a dir, eg "build", and change into it.
- Now (eg.) run:

cmake .. -DCMAKE_INSTALL_PREFIX:STRING="/usr" -DCMAKE_C_FLAGS:STRING="-march=athlon-xp -O2"
make
make install

  If you exported your CFLAGS, there is no need to specify them
  explicitly, as CMake will pick them up.

- OpenAL should get installed as you got used to it.

You could also use ccmake for a console gui, or the Qt CMake gui.
An out-of-source build is recommended. This works on Windows, too,
generating a Visual Studio solution for whatever version you
like.

Some Tips:

- If you want more output at compile time, use

make VERBOSE=1

- If you want to install to a different directory (using same prefix),
  use

make install DESTDIR=/foo/bar

- CMake doesn't has a distclean target by default, so you better
  really do an out-of-tree build, then you can simply delete its
  content when you want a distclean... Furthermore it is easier to
  have different builds using different parameters via out-of-tree
  builds.


Building and installing on *nix-alikes using autotools
======================================================

To build and run ALUT, you need an OpenAL somewhere (and a sound card, of course
:-). If OpenAL is not installed in a standard place (let's say ${OPENAL_PREFIX})
you'll have to tell the build system about it:

   CPPFLAGS="-I${OPENAL_PREFIX}/include"
   LDFLAGS="-L${OPENAL_PREFIX}/lib"
   export CPPFLAGS LDFLAGS

If you install ALUT from the CVS repository, you'll have to call

   ./autogen.sh

first, this will generate some necessary files for building ALUT.

If these initial preparations are done, you can simply run:

   ./configure
   make
   make install

The last step is optional and must be done as 'root', you can build and use the
test programs in 'test_suite' and the example programs in 'examples' without
installing ALUT. To build these, just go to the respective subdirectories and
run 'make'.

Alternatively, there is a spec file in the admin/RPM subdirectory, which can be
used to build an RPM.


Building with VisualStudio
==========================

The admin subdirectory contains VisualStudio6 and VisualStudioDotNET
subdirectories, which contain a Visual Studio 6 workspace and a Visual Studio
.NET solution, respectively. These expect the OpenAL SDK (header and import
library) at the standard path "C:\Program Files\OpenAL 1.1 with EFX SDK". If
this is not the case for your setup, you can change this either directly with a
text editor in the project files (the only way for *real* men ;-) or via the
project settings within Visual Studio. When everything is OK, just build
e.g. the hello_world project and you should get the ALUT library itself and a
demo with Steve Baker saying an immortal phrase.

Support
=======

Please don't email me the authors directly, sign up to the OpenAL user's mailing
list. Instructions can be found here:

   http://www.openal.org/


Licensing
=========

These libraries are offered as OpenSource freeware under the provisions of the
GNU Library Public License (LGPL). Those terms and conditions are detailed in
the file 'COPYING' in this directory.


Portability and Dependencies
============================

This library is designed to be easily portable (and more importantly so that
programs using it can be UTTERLY portable).


Compiling/Linking with ALUT
===========================

In your source code, simply add:

   #include <AL/alut.h>

On *nix-alikes, use

   openal-config --cflags
   openal-config --libs

to find out the compiler and linker flags needed for OpenAL, the flags needed
for ALUT itself depend on the --prefix option given to configure.

In a VisualStudio project, add the path to the ALUT header, import library and
DLL to your project settings as usual.


Directories
===========

Within this distribution, you will find the following directories:

admin      -- Administrative files needed for building ALUT
doc        -- The ALUT documentation
examples   -- For application authors, here is a small collection of example
              programs that you can use to get you started
include    -- The ALUT header
src        -- The source code for the ALUT library itself
test_suite -- A small suite of test programs to check that ALUT is installed and
              working correctly
