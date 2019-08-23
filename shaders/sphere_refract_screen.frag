uniform sampler2D frame_buffer_tex;
uniform float ref_ix = 1.1; // index of refraction
uniform float radius = 0.5; // sphere radius in screen space
uniform vec3  center = vec3(0.5, 0.5, 0.0); // sphere center in screen space
uniform float intensity    = 1.0;
uniform float aspect_ratio = 1.0;

in vec2 tc;

void main() { // ignore z-value (depth) for now

	vec2  scaled_sr = radius * vec2(1.0, aspect_ratio);
	vec2  dir = (tc - center.xy)/scaled_sr; // distance to sphere center relative to sphere radius
	float dsq = dot(dir, dir); // squared distance from sphere center

	if (dsq < 1.0) { // inside sphere - slow case (uncommon)
		float dist = sqrt(dsq);
		vec2 scaled_dist = dist*scaled_sr;
		vec3 normal = vec3(dir, sqrt(1.0 - dsq));
		vec3 refract_dir = -refract(vec3(0,0,1), normal, 1.0/ref_ix); // refraction into the sphere, inverted
		vec2 delta = (tc - center.xy - refract_dir.xy*scaled_dist)*intensity;
#ifdef CHROMATIC_REFRACT
		float dist_sqrt = sqrt(dist);
		vec2 R_pos = tc + (1.0 - 0.5*dist_sqrt)*delta;
		vec2 G_pos = tc + (1.0 + 0.0*dist_sqrt)*delta;
		vec2 B_pos = tc + (1.0 + 0.5*dist_sqrt)*delta;
		fg_FragColor = vec4(texture(frame_buffer_tex, R_pos).r, texture(frame_buffer_tex, G_pos).g, texture(frame_buffer_tex, B_pos).b, 1.0);
#else
		fg_FragColor = vec4(texture(frame_buffer_tex, (tc + delta)).rgb, 1.0);
#endif // CHROMATIC_REFRACT
	}
	else { // outside sphere - fast case (common)
		fg_FragColor = vec4(texture(frame_buffer_tex, tc).rgb, 1.0);
	}
}
