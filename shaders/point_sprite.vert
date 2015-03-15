uniform float point_scale = 1.0; // or point size in pixels, for constant size
in float point_size; // world space

void main()
{
	vec4 epos    = fg_ModelViewMatrix * fg_Vertex;
	gl_Position  = fg_ProjectionMatrix * epos;
#ifdef CONSTANT_PT_SIZE
	gl_PointSize = point_scale;
#else
	gl_PointSize = clamp(point_scale*point_size/length(epos.xyz), 1.0, 64.0);
#endif

#ifdef ENABLE_LIGHTING
	vec3 normal = normalize(fg_NormalMatrix * fg_Normal); // eye space
	vec3 color  = vec3(0.0);
	if (enable_light0) color += add_light_comp_pos0(normal, epos).rgb;
	if (enable_light1) color += add_light_comp_pos1(normal, epos).rgb;
	fg_Color_vf = vec4(color, fg_Color.a);
#else
	fg_Color_vf = fg_Color;
#endif
} 
