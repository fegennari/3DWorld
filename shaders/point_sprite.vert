uniform float point_scale = 1.0;
attribute float point_size;

void main()
{
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
	vec4 eye      = gl_ModelViewMatrixInverse[3]; // world space
#ifdef CONSTANT_PT_SIZE
	gl_PointSize  = point_scale;
#else
	gl_PointSize  = clamp(point_scale*point_size/distance(eye.xyz, gl_Vertex.xyz), 1.0, 64.0);
#endif
} 
