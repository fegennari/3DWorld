uniform float point_scale = 1.0;
attribute float point_size;

void main()
{
	gl_Position  = ftransform();
	vec4 eye     = gl_ModelViewMatrixInverse[3]; // world space
#ifdef CONSTANT_PT_SIZE
	gl_PointSize = point_scale;
#else
	gl_PointSize = clamp(point_scale*point_size/distance(eye.xyz, gl_Vertex.xyz), 1.0, 64.0);
#endif

#ifdef ENABLE_LIGHTING
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	vec3 color  = vec3(0.0);
	if (enable_light0) color += add_light_comp0(normal).rgb;
	if (enable_light1) color += add_light_comp1(normal).rgb;
	gl_FrontColor = vec4(color, gl_Color.a);
#else
	gl_FrontColor = gl_Color;
#endif
} 
