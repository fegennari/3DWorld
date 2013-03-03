uniform sampler2D color_map, normal_map;
uniform vec4 color_scale = vec4(1,1,1,1);

varying vec4 eye_space_pos;

void main()
{
	vec4 texel = texture2D(color_map, gl_TexCoord[0].st);
	if (texel.a < 0.75) discard; // transparent
	//if (normal.w == 0.0) discard; // normal not written to (uses nearest filter)
	check_noise_and_maybe_discard(0.0, gl_Color.a);

	// transform the normal into eye space, but don't normalize because it may be scaled for shadows
	vec3 normal = gl_NormalMatrix * (2.0*texture2D(normal_map, gl_TexCoord[0].st).xyz - vec3(1,1,1));
	if (dot(normal, eye_space_pos.xyz) > 0.0) normal = -normal; // facing away from the eye, so reverse (could use faceforward())
	
	vec4 color = gl_Color * gl_LightModel.ambient;
	const bool shadowed = false;
	if (enable_light0) color += add_leaf_light_comp(shadowed, normal, eye_space_pos, 0);
	if (enable_light1) color += add_leaf_light_comp(shadowed, normal, eye_space_pos, 1);
	gl_FragColor = apply_fog(clamp(color*color_scale, 0.0, 1.0)*vec4(texel.rgb, 1.0));
}
