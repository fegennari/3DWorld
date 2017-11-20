uniform float radius    = 1.0;
uniform float intensity = 1.0;
uniform sampler2D tex0;
uniform vec4 colorA, colorB;

in vec3 world_space_pos;
in vec2 tc;

void main()
{
	float dist   = length(world_space_pos);
	vec3 spos    = world_space_pos*radius/sqrt(dist);
	vec3 tv      = world_space_pos/dist + vec3(100.0); // direction + constant offset
	float tsaw   = fract(0.25*time) - 0.5; // sawtooth wave
	float ttri   = abs(tsaw); // triangle wave
	float time2  = (tsaw - 0.005*sign(ttri)*dist + 6.0); // add constant offset + radius offset
	vec3 ftime   = 0.02*vec3(fract_smooth(tv.x*time2), fract_smooth(tv.y*time2), fract_smooth(tv.z*time2));
	fg_FragColor = texture(tex0, tc) * mix(colorA, colorB, gen_cloud_alpha_time(0.18*spos, ftime));
	fg_FragColor.a *= intensity;
}
