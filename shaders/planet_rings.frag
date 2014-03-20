uniform vec3 planet_pos, sun_pos, camera_pos;
uniform float planet_radius, ring_ri, ring_ro, sun_radius, bf_draw_sign;
uniform sampler2D noise_tex, particles_tex;
uniform sampler1D ring_tex;

varying vec3 normal, world_space_pos, vertex;
varying vec2 tc;


vec3 add_light_rings(in vec3 n, in vec4 epos)
{
	vec3 light_dir = normalize(gl_LightSource[0].position.xyz - epos.xyz); // normalize the light's direction in eye space
	vec3 diffuse   = (gl_Color * gl_LightSource[0].diffuse).rgb;
	vec3 ambient   = (gl_Color * gl_LightSource[0].ambient).rgb;
	vec3 specular  = get_light_specular(n, light_dir, epos.xyz, gl_LightSource[0].specular.rgb);
	float atten    = calc_light_atten0(epos);
	if (sun_radius > 0.0) {atten *= calc_sphere_shadow_atten(world_space_pos, sun_pos, sun_radius, planet_pos, planet_radius);} // sun exists
	return (ambient + (abs(dot(n, light_dir))*diffuse + specular)) * atten;
}

void main()
{
	if (bf_draw_sign*(length(world_space_pos - camera_pos) - length(planet_pos - camera_pos)) < 0.0) discard; // on the wrong side of the planet

	float rval = clamp((length(vertex) - ring_ri)/(ring_ro - ring_ri), 0.0, 1.0);
	vec4 texel = texture1D(ring_tex, rval);
	if (texel.a == 0.0) discard;

	// alpha lower when viewing edge
	vec4 epos = gl_ModelViewMatrix * vec4(vertex, 1.0);
	texel.a *= pow(abs(dot(normal, normalize(epos.xyz))), 0.2); // 5th root

	vec2 tcs   = 16*tc;
	vec3 norm2 = normalize(normal + vec3(texture2D(noise_tex, tcs).r-0.5, texture2D(noise_tex, tcs+vec2(0.4,0.7)).r-0.5, texture2D(noise_tex, tcs+vec2(0.3,0.8)).r-0.5));
	vec4 color = vec4(0,0,0,1);
	color.rgb += add_light_rings(norm2, epos); // ambient, diffuse, and specular
	color.rgb += (gl_Color * gl_LightSource[1].ambient).rgb; // ambient only

	float alpha = texture2D(particles_tex, 23 *tc).r;
	alpha      += texture2D(particles_tex, 42 *tc).r;
	alpha      += texture2D(particles_tex, 75 *tc).r;
	alpha      += texture2D(particles_tex, 133*tc).r;
	if (alpha == 0.0) discard;
	color.a    *= min(1.0, alpha); // add 4 octaves of random particles
	color.a     = min(1.0, 2.0*color.a); // increase alpha to make alpha_to_coverage mode look better
	gl_FragColor = color * texel;
}
