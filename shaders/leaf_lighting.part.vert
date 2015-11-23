uniform float normal_scale = 1.0;
uniform float water_depth = 0.0;
uniform vec4 color_scale = vec4(1.0);
uniform vec3 world_space_offset = vec3(0.0);

void calc_leaf_lighting() {
	// transform the normal into eye space, but don't normalize because it may be scaled for shadows
	vec3 normal = fg_NormalMatrix * fg_Normal * normal_scale;
	
	vec4 epos    = fg_ModelViewMatrix * fg_Vertex;
	float nscale = ((dot(normal, epos.xyz) > 0.0) ? -1.0 : 1.0); // facing away from the eye, so reverse (could use faceforward())
	normal *= nscale;

	vec3 vpos      = fg_Vertex.xyz + world_space_offset;
	vec3 ws_normal = nscale*normalize(fg_Normal);
	vec3 color     = get_tree_leaf_lighting(epos, normal, vpos, ws_normal);
	fg_Color_vf    = vec4(min(2.0*fg_Color.rgb, clamp(color*color_scale.rgb, 0.0, 1.0)), 1.0); // limit lightning color

	if (underwater) {
		vec3 eye        = fg_ModelViewMatrixInverse[3].xyz;
		gl_FogFragCoord = length(mix(eye, fg_Vertex.xyz, min(1.0, water_depth/max(0.0001, (fg_Vertex.z - eye.z)))) - eye);
	}
	else {
		gl_FogFragCoord = length(epos.xyz);
	}
}