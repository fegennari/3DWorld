varying vec3 normal;
varying vec4 epos, proj_pos;
uniform sampler2D reflection_tex, water_normal_tex, height_tex;
uniform vec4 water_color, reflect_color;
uniform float noise_time, wave_time, wave_amplitude, water_plane_z, water_green_comp, reflect_scale;
uniform float x1, y1, x2, y2, zmin, zmax;

vec3 water_normal_lookup(in vec2 tc) {
	return 2.0*(texture2D(water_normal_tex, 0.5*tc).rgb - 0.5);
}

vec3 get_wave_normal(in vec2 tc) {
	float ntime = 2.0*abs(fract(0.005*wave_time) - 0.5);
	return mix(water_normal_lookup(tc), water_normal_lookup(tc + vec2(0.3, 0.6)), ntime);
}

void main()
{
	float mesh_z = zmin + (zmax - zmin)*texture2D(height_tex, gl_TexCoord[1].st).r;
	if (water_plane_z < mesh_z) discard;

	vec3 norm   = normalize(normal); // renormalize
	vec2 ripple = vec2(0,0);

	if (add_noise) { // for rain
		vec3 wave_n = get_wave_normal(fract(4.61*noise_time)*3.0*proj_pos.xy/proj_pos.w);
		norm        = normalize(norm + 0.06*gl_NormalMatrix*wave_n);
		ripple     += 0.025*wave_n.xy;
	}
	vec3 light_norm = norm;

	if (add_waves) {
		// calculate ripple adjustment of normal and reflection based on scaled water normal map texture
		vec3 wave_n     = wave_amplitude*get_wave_normal(gl_TexCoord[0].st);
		vec3 wave_n_eye = gl_NormalMatrix*wave_n;
		light_norm  = normalize(norm + wave_n_eye);
		norm        = normalize(norm + 0.1*wave_n_eye); // lower scale for fresnel term
		ripple     += 0.05*wave_n.xy;
	}

	// calculate lighting
	vec3 epos_n   = normalize(epos.xyz);
	float cos_view_angle = abs(dot(epos_n, norm));
	vec4 color    = water_color;
	color.a      *= mix(1.0, min(1.0, 20.0*(water_plane_z - mesh_z)), min(1.0, 2.5*cos_view_angle)); // blend to alpha=0 near the shore
	vec4 lighting = gl_FrontMaterial.emission + gl_FrontMaterial.ambient * gl_LightModel.ambient;
	if (enable_light0) {lighting += add_light_comp_pos(light_norm, epos, 0);}
	if (enable_light1) {lighting += add_light_comp_pos(light_norm, epos, 1);}

	if (reflections) {
		// add some green at shallow view angles
		color = mix(color, vec4(0.0, 1.0, 0.5, color.a), water_green_comp*(1.0 - cos_view_angle));

		// calculate reflections
		float reflect_w  = reflect_scale*get_fresnel_reflection(-epos_n, norm, 1.0, 1.333);
		vec2 ref_tex_st  = clamp(0.5*proj_pos.xy/proj_pos.w + 0.3*ripple + vec2(0.5, 0.5), 0.0, 1.0);
		vec4 reflect_tex = vec4(texture2D(reflection_tex, ref_tex_st).rgb, 1.0);
		color = mix(color, reflect_color * reflect_tex, reflect_w);
	}

	// determine final color with fog
	vec4 frag_color = vec4(color.rgb * lighting.rgb, color.a * gl_FrontMaterial.diffuse.a); // use diffuse alpha directly
	gl_FragColor = apply_fog_epos(frag_color, epos);
}
