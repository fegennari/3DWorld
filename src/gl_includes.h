// 3D World - Cross-Platform OpenGL Includes
// by Frank Gennari
// 8/26/11

#ifndef _GL_INCLUDES_H_
#define _GL_INCLUDES_H_

#if ((defined(__MACH__))&&(defined(__APPLE__)))
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <OpenGL/glext.h>
#include <GLU/glu.h>
#else
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/glu.h>
//#include <GL/glext.h>
#endif

#endif // _GL_INCLUDES_H_

