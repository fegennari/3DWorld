uniform vec3 planet_pos, sun_pos, camera_pos;
uniform float planet_radius, ring_ri, ring_ro, sun_radius, bf_draw_sign;
uniform sampler2D noise_tex;
uniform sampler1D ring_tex;
varying vec4 epos;
varying vec3 normal, world_space_pos, vertex;


vec4 add_light_rings(in vec3 n)
{
	vec3 light_dir = normalize(gl_LightSource[0].position.xyz - epos.xyz); // normalize the light's direction in eye space
	vec4 diffuse   = gl_Color * gl_LightSource[0].diffuse;
	vec4 ambient   = gl_Color * gl_LightSource[0].ambient;
	vec4 specular  = get_light_specular(n, light_dir, epos.xyz, 0);
	float atten    = calc_light_atten(epos, 0);
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
	texel.a *= pow(abs(dot(normal, normalize(epos.xyz))), 0.2); // 5th root

	vec2 tc = 16*gl_TexCoord[0].st;
	vec3 norm2 = normalize(normal + vec3(texture2D(noise_tex, tc).r-0.5, texture2D(noise_tex, tc+vec2(0.4,0.7)).r-0.5, texture2D(noise_tex, tc+vec2(0.3,0.8)).r-0.5));
	vec4 color = gl_FrontMaterial.emission;
	color.rgb += add_light_rings(norm2).rgb; // ambient, diffuse, and specular
	color.rgb += (gl_Color * gl_LightSource[1].ambient).rgb; // ambient only
	color.a   *= texture2D(noise_tex, 25*gl_TexCoord[0].st).r;
	gl_FragColor = apply_fog(color*texel);
}
