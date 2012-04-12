// 3D World - FILE Utility Functions
// by Frank Gennari
// 3/6/12

#ifndef _FILE_UTILS_H_
#define _FILE_UTILS_H_

#include <fstream>
#include "3DWorld.h"

inline bool read_int  (FILE *fp, int      &val) {return (fscanf(fp, "%i", &val) == 1);}
inline bool read_uint (FILE *fp, unsigned &val) {return (fscanf(fp, "%u", &val) == 1);}
inline bool read_nonzero_uint(FILE *fp, unsigned &val) {return (fscanf(fp, "%u", &val) == 1 && val > 0);}
inline bool read_float(FILE *fp, float    &val) {return (fscanf(fp, "%f", &val) == 1);}
inline bool read_str  (FILE *fp, char     *val) {return (fscanf(fp, "%s",  val) == 1);}

inline bool read_vector(FILE *fp, vector3d &v) { // or point
	return (fscanf(fp, "%f%f%f", &v.x, &v.y, &v.z) == 3);
}

inline bool read_color(FILE *fp, colorRGBA &c) {
	c.A = 1.0; // default
	return (fscanf(fp, "%f%f%f%f", &c.R, &c.G, &c.B, &c.A) >= 3); // alpha is optional
}


inline bool read_bool (FILE *fp, bool     &val) {
	int tmp;
	if (fscanf(fp, "%i", &tmp) != 1) return 0;
	val = (tmp != 0);
	return 1;
}

#endif // _FILE_UTILS_H_

