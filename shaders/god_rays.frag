uniform vec3 sun_pos;
uniform vec4 sun_color;
uniform float aspect_ratio = 1.0;
uniform sampler2D depth_tex;

uniform float exposure = 1.0;
uniform float decay    = 0.988;
const int NUM_SAMPLES  = 64;

in vec2 tc;

float get_sample(in vec2 pos) {
	float depth = texture(depth_tex, pos).r;
	vec2 v      = (sun_pos.xy - pos)*vec2(aspect_ratio, 1.0);
	float dstsq = dot(v, v); // distance squared from light source to current pixel in screen space
	float val   = 0.0; // declare sample color to that of an occluder
	if (depth > 0.999) {val = mix(0.3, 1.0, exp(-40.0*dstsq));} // if not occluded, mix the color intensity between atmosphere and sun with Gaussian sun profile
	return val;
}

void main()
{
	vec2 delta = (tc - sun_pos.xy) / float(NUM_SAMPLES); // tc is screen space pos
	vec2 pos   = tc;
	float illumination_decay = 1.0;
	float weight = 0.0;
	
	for (int i = 0; i < NUM_SAMPLES; i++) {
		pos -= delta;
		weight += get_sample(pos) * illumination_decay;
		illumination_decay *= decay;
	}
	gl_FragColor = sun_color * weight * exposure / float(NUM_SAMPLES);
}
