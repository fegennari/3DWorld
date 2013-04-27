uniform int num_dlights = 0;
uniform float normal_scale = 1.0;
uniform vec4 color_scale = vec4(1.0);

void calc_leaf_lighting()
{
	// transform the normal into eye space, but don't normalize because it may be scaled for shadows
	vec3 normal = gl_NormalMatrix * gl_Normal * normal_scale;
	
	vec4 eye_space_pos = gl_ModelViewMatrix * gl_Vertex;
	if (dot(normal, eye_space_pos.xyz) > 0.0) normal = -normal; // facing away from the eye, so reverse (could use faceforward())
	
	// Compute the globalAmbient term
	bool shadowed = (sqrt(dot(gl_Normal, gl_Normal)) < 0.4);
	vec4 color    = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_leaf_light_comp(shadowed, normal, eye_space_pos, 0);
	if (enable_light1) color += add_leaf_light_comp(shadowed, normal, eye_space_pos, 1);
#ifndef NO_LEAF_DLIGHTS
	vec3 n = normalize(normal);

	for (int i = 2; i < num_dlights+2; ++i) { // add N dynamic point light sources
		color += add_pt_light_comp(n, eye_space_pos, i);
	}
#endif
	gl_FrontColor   = min(2*gl_Color, clamp(color*color_scale, 0.0, 1.0)); // limit lightning color
	gl_FogFragCoord = length(eye_space_pos.xyz);
}