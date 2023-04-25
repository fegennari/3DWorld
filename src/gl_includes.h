// 3D World - Cross-Platform OpenGL Includes
// by Frank Gennari
// 8/26/11
#pragma once

#if ((defined(__MACH__))&&(defined(__APPLE__)))
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <OpenGL/glext.h>
#include <GLU/glu.h>
#else
#include <GL/glew.h>
#include <GL/gl.h>

#ifdef _WIN32
#ifdef __MINGW32__
#include <GL/freeglut.h>
#else
#include <freeglut.h>
#endif
#undef FAR // undefine conflicting defines picked up from windows headers
#else // linux
#include <GL/freeglut.h> // standard glut
#endif

#include <GL/glu.h>

//#include <GL/glext.h>
#endif

