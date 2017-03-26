uniform sampler2D frame_buffer_tex;
uniform float ref_ix = 1.1; // index of refraction
uniform float radius = 0.5; // sphere radius in screen space
uniform vec3  center = vec3(0.5, 0.5, 0.0); // sphere center in screen space
uniform float intensity    = 1.0;
uniform float aspect_ratio = 1.0;

in vec2 tc;

void apply_sphere_refract(inout vec2 pos, in vec2 sc, in float sr, in float intensity) {
	
	vec2  scaled_sr = sr * vec2(1.0, aspect_ratio);
	vec2  dir = (tc - sc)/scaled_sr; // distance to sphere center relative to sphere radius
	float dsq = dot(dir, dir); // squared distance from sphere center

	if (dsq < 1.0) {
		vec3 normal = vec3(dir, sqrt(1.0 - dsq));
		vec3 refract_dir = -refract(vec3(0,0,1), normal, 1.0/ref_ix); // refraction into the sphere, inverted
		pos += tc - mix(pos, (sc + refract_dir.xy*sqrt(dsq)*scaled_sr), intensity);
	}
}

void main() {

	vec2  pos = tc;
	apply_sphere_refract(pos, center.xy, radius, intensity); // ignore z-value (depth) for now
	fg_FragColor = vec4(texture(frame_buffer_tex, pos).rgb, 1.0);
}
