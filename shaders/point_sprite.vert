uniform float point_scale = 1.0; // or point size in pixels, for constant size
in float point_size; // world space

void main()
{
#ifdef CONSTANT_PT_SIZE
	gl_Position  = fg_ftransform();
	gl_PointSize = point_scale;
#else
	vec4 epos    = fg_ModelViewMatrix * fg_Vertex;
	gl_Position  = fg_ProjectionMatrix * epos;
	gl_PointSize = clamp(point_scale*point_size/length(epos.xyz), 1.0, 64.0);
#endif

#ifdef ENABLE_LIGHTING
	vec3 normal = normalize(fg_NormalMatrix * fg_Normal); // eye space
	vec3 color  = vec3(0.0);
	if (enable_light0) color += add_light_comp0(normal).rgb;
	if (enable_light1) color += add_light_comp1(normal).rgb;
	gl_FrontColor = vec4(color, fg_Color.a);
#else
	gl_FrontColor = fg_Color;
#endif
} 
