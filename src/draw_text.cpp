// 3D World - Drawing Code
// by Frank Gennari
// 3/10/02
#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"


extern int window_height;


void draw_text(colorRGBA const &color, float x, float y, float z, char const *text, float tsize, bool bitmap_font) {

	//bitmap_font |= ((display_mode & 0x80) != 0);
	shader_t s;
	glDisable(GL_DEPTH_TEST);

	if (bitmap_font) { // FIXME: uses fixed function pipeline, need to replace glRasterPos and glutBitmapCharacter with something else
		color.do_glColor();
		glRasterPos3f(x, y, z);
	}
	else {
		s.begin_color_only_shader(color);
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		glPushMatrix();
		glTranslatef(x, y, z);
		uniform_scale(0.000005*tsize);
	}
	unsigned line_num(0);

	while (*text) {
		if (*text == '\n') { // newline (CR/LF)
			++line_num;

			if (bitmap_font) {
				glRasterPos3f(x, y-(0.5*line_num)/window_height, z);
			}
			else {
				glPopMatrix();
				glPushMatrix();
				glTranslatef(x, y-0.001*line_num*tsize, z);
				uniform_scale(0.000005*tsize);
			}
		}
		else {
			if (bitmap_font) {
				glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *text); // other fonts available
			}
			else {
				glutStrokeCharacter(GLUT_STROKE_ROMAN, *text); // GLUT_STROKE_MONO_ROMAN
			}
		}
		text++;
	}
	if (!bitmap_font) {
		glPopMatrix();
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
		s.end_shader();
	}
	glEnable(GL_DEPTH_TEST);
}

