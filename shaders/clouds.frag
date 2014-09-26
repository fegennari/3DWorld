uniform float water_plane_z;
uniform vec3 cloud_offset = vec3(0.0);
uniform vec4 sun_color;
uniform vec3 sun_pos, camera_pos;

in vec3 vertex;
in vec4 color;

// Note: these two functions assume the sun is a directional light
vec4 apply_fog_colored(in vec4 color, in vec3 vertex, in float cscale) {
	float fog_scale = gl_FogFragCoord*get_custom_fog_scale(vertex.z);
	return apply_fog_colored(color, fog_scale, normalize(vertex - camera_pos), normalize(sun_pos - camera_pos), cscale);
}

void main()
{
	vec3 light_dir = normalize(sun_pos - vertex);
	float dp    = max(0.0, dot(normalize(vertex - camera_pos), light_dir));
	vec4 color2 = 0.75*color + vec4(0.5*sun_color.rgb*pow(dp, 8.0), 0.0); // add sun glow
	vec3 pos    = vertex + cloud_offset;
	float eye_z = camera_pos.z; // world space
	float t     = min(1.0, (eye_z - water_plane_z)/max(0.0, (eye_z - pos.z)));
	float black_mix = ((underwater_atten && t > 0.0 && t < 1.0) ? 1.0 : 0.0);
	float alpha     = gen_cloud_alpha(pos.xy);
	float lcolor    = 1.0 - 0.1*alpha;

#ifdef CLOUD_LIGHTING // looks nice, but very slow (but could use fewer octaves if needed)
	float shadow = 0.0;
	for (int i = 0; i < 4; ++i) {
		shadow += gen_cloud_alpha(pos.xy + 0.5*i*light_dir.xy);
	}
	lcolor *= 1.0 - 0.05*shadow;
#endif
	vec4 cscale     = mix(vec4(lcolor, lcolor, lcolor, clamp(alpha, 0.0, 1.0)), vec4(0,0,0,1), black_mix);
	//fg_FragColor    = apply_fog_scaled(color2*cscale, vertex.z);
	fg_FragColor    = apply_fog_colored(color2*cscale, vertex, clamp((1.0 + 0.05*(vertex.z - eye_z)), 0.0, 1.0));
}
