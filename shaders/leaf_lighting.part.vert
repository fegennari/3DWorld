uniform float normal_scale = 1.0;
uniform float water_depth = 0.0;
uniform vec4 color_scale = vec4(1.0);
uniform vec3 world_space_offset = vec3(0.0);

// for indir lighting
uniform sampler3D smoke_and_indir_tex;
uniform vec3 const_indir_color = vec3(0.0);

void add_indir_lighting(inout vec3 lit_color, in vec3 vpos) {
	vec3 indir_color = const_indir_color; // add constant indir

	if (indir_lighting) {
		vec3 spos    = clamp((vpos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
		indir_color += texture(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
	}
	lit_color += gl_Color.rgb * indir_color;
}

void calc_leaf_lighting() {
	// transform the normal into eye space, but don't normalize because it may be scaled for shadows
	vec3 normal = fg_NormalMatrix * fg_Normal * normal_scale;
	
	vec4 eye_space_pos = fg_ModelViewMatrix * fg_Vertex;
	float nscale = ((dot(normal, eye_space_pos.xyz) > 0.0) ? -1.0 : 1.0); // facing away from the eye, so reverse (could use faceforward())
	normal *= nscale;
	
	vec3 color    = vec3(0.0);
	vec3 vpos     = fg_Vertex.xyz + world_space_offset;
	float ambient_scale = (indir_lighting ? 0.0 : 1.0);
	float diffuse_scale = sqrt(dot(fg_Normal, fg_Normal));
	add_indir_lighting(color, vpos);
	if (enable_light0)  {color += add_leaf_light_comp(normal,  eye_space_pos, 0, ambient_scale, diffuse_scale).rgb;}
	if (enable_light1)  {color += add_leaf_light_comp(normal,  eye_space_pos, 1, ambient_scale, diffuse_scale).rgb;}
	if (enable_light2)  {color += add_pt_light_comp  (normalize(normal), eye_space_pos, 2).rgb;} // lightning
	if (enable_dlights) {add_dlights(color, vpos, nscale*normalize(fg_Normal), vec3(1.0));}
	fg_Color_vf = vec4(min(2.0*fg_Color.rgb, clamp(color*color_scale.rgb, 0.0, 1.0)), 1.0); // limit lightning color

	if (underwater) {
		vec3 eye        = fg_ModelViewMatrixInverse[3].xyz;
		gl_FogFragCoord = length(mix(eye, fg_Vertex.xyz, min(1.0, water_depth/max(0.0001, (fg_Vertex.z - eye.z)))) - eye);
	}
	else {
		gl_FogFragCoord = length(eye_space_pos.xyz);
	}
}