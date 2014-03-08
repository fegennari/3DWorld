uniform sampler2D color_map, normal_map;
uniform vec3 ref_dir         = vec3(0,1,0);
uniform vec2 normal_tc_off   = vec2(0,0);
uniform vec2 normal_tc_scale = vec2(1,1);

varying vec4 world_space_pos, eye_space_pos;
varying vec2 tc;

void main()
{
	vec2 tc_scaled = normal_tc_scale*tc;
	vec4 texel = texture2D(color_map, tc_scaled);
	if (texel.a < 0.5) discard; // transparent
	check_noise_and_maybe_discard(0.0, gl_Color.a);

	// transform normal into billboard orientation 
	vec3 normal = 2.0*texture2D(normal_map, (tc_scaled + normal_tc_off)).xyz - vec3(1.0);
	normal.y *= -1.0; // texture is rendered with ybot < ytop
	vec4 eye  = gl_ModelViewMatrixInverse[3]; // world space
	vec3 vdir = eye.xyz - world_space_pos.xyz;
	vec2 rd_n = normalize(ref_dir.xy);
	vec2 vd_n = normalize(vdir.xy);
	float dp  = dot(rd_n, vd_n);
	float s   = length(cross(vec3(rd_n, 0), vec3(vd_n, 0))) * ((cross(vdir, ref_dir).z < 0.0) ? -1.0 : 1.0);
	mat3 mrot = mat3(dp, -s, 0.0,  s, dp, 0.0,  0.0, 0.0, 1.0);
	normal    = normalize(gl_NormalMatrix * (mrot * normal)); // convert to eye space
	
	vec4 color  = vec4(0,0,0,1);
	if (enable_light0) color.rgb += add_light_comp0(normal).rgb;
	if (enable_light1) color.rgb += add_light_comp1(normal).rgb;
	if (enable_light2) color.rgb += add_light_comp (normal, 2).rgb * calc_light_atten(eye_space_pos, 2);
	gl_FragColor = apply_fog_epos(clamp(color, 0.0, 1.0)*vec4(texel.rgb, 1.0), eye_space_pos);
}
