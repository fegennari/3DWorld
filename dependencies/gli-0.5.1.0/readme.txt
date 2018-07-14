================================================================================
OpenGL Mathematics (GLI)
--------------------------------------------------------------------------------
www.g-truc.net
mail@g-truc.net

================================================================================
The MIT License
--------------------------------------------------------------------------------
Copyright (c) 2010 - 2013 G-Truc Creation (www.g-truc.net)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

================================================================================
GLI Usage
--------------------------------------------------------------------------------
GLI is a header only library, there is nothing to build, just include it.
#include <gli/gli.hpp>

================================================================================
GLI 0.5.1.0: 2014-01-18
--------------------------------------------------------------------------------
- Added flip function
- Added level_count function
- Fixed interaction with std::map (#33)
- Added texelFetch and texelWrite functions

================================================================================
GLI 0.5.0.0: 2013-11-24
--------------------------------------------------------------------------------
- Essencially a rewrite of the library
- Added explicit copies
- Added single memory allocation per texture storage
- Added texture views
- Added texture copies
- Added comparison operators
- Added clear

================================================================================
GLI 0.4.1.0: 2013-03-10
--------------------------------------------------------------------------------
- Added DDS saving
- Fixed GCC build
- Fixed XCode build

================================================================================
GLI 0.4.0.0: 2013-01-28
--------------------------------------------------------------------------------
- Large API refactoring
- Performance improvements at loading: 50x in debug and 50% in release build
- Added texture2DArray
- Added textureCube and textureCubeArray
- Added texture3D
- Added texture1D and texture1DArray
- Improved DDS loading support

================================================================================
GLI 0.3.0.3: 2011-04-05
--------------------------------------------------------------------------------
- Fixed bugs

================================================================================
GLI 0.3.0.2: 2011-02-08
--------------------------------------------------------------------------------
- Fixed bugs

================================================================================
GLI 0.3.0.1: 2010-10-15
--------------------------------------------------------------------------------
- Fixed bugs

================================================================================
GLI 0.3.0.0: 2010-10-01
--------------------------------------------------------------------------------
- Added DDS10 load and store (BC1 - BC7)
- Added extension system
- Added automatic OpenGL texture object creation from file

================================================================================
GLI 0.2.2.0: 2010-09-07
--------------------------------------------------------------------------------
- Added DDS exporter

================================================================================
GLI 0.2.1.1: 2010-05-12
--------------------------------------------------------------------------------
- Fixed GCC build

================================================================================
GLI 0.2.1.0: 2010-02-15
--------------------------------------------------------------------------------
- Added texelWrite function
- Fixed Visual Studio 2010 warnings
- Added readme.txt and copying.txt

================================================================================
GLI 0.2.0.0: 2010-01-10
--------------------------------------------------------------------------------
- Updated API
- Removed Boost dependency

================================================================================
GLI 0.1.1.0: 2009-09-18
--------------------------------------------------------------------------------
- Fixed DDS loader
- Added RGB8 to DDS loader
- Added component swizzle operation
- Added 32 bits integer components support
- Fixed texel fetch

================================================================================
GLI 0.1.0.0: 2009-03-28
--------------------------------------------------------------------------------
- Added TGA loader
- Added DDS loader
- Added duplicate, crop, partial copy
- Added mipmaps generation.

================================================================================


