uniform int num_dlights = 0;
uniform float normal_scale = 1.0;
uniform vec4 color_scale = vec4(1.0);
uniform vec3 world_space_offset = vec3(0.0);

void calc_leaf_lighting()
{
	// transform the normal into eye space, but don't normalize because it may be scaled for shadows
	vec3 normal = gl_NormalMatrix * fg_Normal * normal_scale;
	
	vec4 eye_space_pos = gl_ModelViewMatrix * fg_Vertex;
	//if (dot(normal, eye_space_pos.xyz) > 0.0) normal = -normal; // facing away from the eye, so reverse (could use faceforward())
	float nscale = ((dot(normal, eye_space_pos.xyz) > 0.0) ? -1.0 : 1.0);
	normal *= nscale;
	
	// Compute the globalAmbient term
	bool shadowed = (sqrt(dot(fg_Normal, fg_Normal)) < 0.4);
	vec4 color    = vec4(0,0,0,1);
	if (enable_light0) color.rgb += add_leaf_light_comp(shadowed, normal,  eye_space_pos, 0).rgb;
	if (enable_light1) color.rgb += add_leaf_light_comp(shadowed, normal,  eye_space_pos, 1).rgb;
	if (enable_light2) color.rgb += add_pt_light_comp  (normalize(normal), eye_space_pos, 2).rgb; // lightning

	if (enable_dlights) {
		vec3 vpos  = fg_Vertex.xyz + world_space_offset;
		color.rgb += add_dlights(vpos, nscale*normalize(fg_Normal), vec3(1.0)).rgb;
	}
	gl_FrontColor   = min(2*fg_Color, clamp(color*color_scale, 0.0, 1.0)); // limit lightning color
	gl_FogFragCoord = length(eye_space_pos.xyz);
}