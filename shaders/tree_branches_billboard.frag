uniform sampler2D color_map, normal_map;
uniform vec3 ref_dir = vec3(0,1,0);

varying vec4 world_space_pos;

void main()
{
	vec4 texel = texture2D(color_map, gl_TexCoord[0].st);
	if (texel.a < 0.9) discard; // transparent
	check_noise_and_maybe_discard(0.0, gl_Color.a);

	// transform normal into billboard orientation 
	vec3 normal = 2.0*texture2D(normal_map, gl_TexCoord[0].st).xyz - vec3(1,1,1);
	vec4 eye    = gl_ModelViewMatrixInverse[3]; // world space
	vec3 vdir   = normalize(eye.xyz - world_space_pos.xyz);
	float angle = acos(dot(normalize(ref_dir.xy), normalize(vdir.xy)));
	if (cross(vdir, ref_dir).z < 0.0) {angle = -angle;} // rotate the other direction
	mat3 mrot   = mat3(cos(angle), -sin(angle), 0.0,
			           sin(angle),  cos(angle), 0.0,
			           0.0,         0.0,        1.0);
	normal      = mrot * normal;
	normal      = normalize(gl_NormalMatrix * normal); // convert to eye space
	
	vec4 color  = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	gl_FragColor = apply_fog(clamp(color, 0.0, 1.0)*vec4(texel.rgb, 1.0));
}
