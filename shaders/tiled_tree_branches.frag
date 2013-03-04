uniform sampler2D tex0;
uniform float opacity = 1.0;

varying vec4 epos;
varying vec3 normal;

void main()
{
	check_noise_and_maybe_discard((1.0 - opacity), 1.0); // inverted value
	vec4 texel   = texture2D(tex0, gl_TexCoord[0].st);
	vec3 normal2 = normalize(normal); // renormalize
	vec4 color   = gl_FrontMaterial.emission + gl_FrontMaterial.ambient * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp_pos(normal2, epos, 0); // sun
	if (enable_light1) color += add_light_comp_pos(normal2, epos, 1); // moon
	if (enable_light2) color += add_light_comp_pos(normal2, epos, 2) * calc_light_atten(epos, 2); // lightning
	gl_FragColor = apply_fog(texel * color);
}
