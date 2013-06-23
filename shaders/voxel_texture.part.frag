uniform vec4 color0 = vec4(1.0);
uniform vec4 color1 = vec4(1.0);
uniform sampler2D tex0, tex1, tex_top;

vec4 get_texture_val(in vec3 normal, in vec3 pos)
{
	vec4 texel0 = lookup_triplanar_texture(pos, normal, tex0, tex0, tex0) * color0;
	vec4 texel1 = lookup_triplanar_texture(pos, normal, tex1, tex1, tex1) * color1;
	// interpolate between the two texture/color pairs using a random noise function
	return mix(texel0, texel1, procedural_eval(pos));
}

vec4 get_texture_val_z_test(in vec3 normal, in vec3 pos)
{
	vec4 texel0 = lookup_triplanar_texture_2sz(pos, normal, tex0, tex0, tex0, tex_top) * color0;
	vec4 texel1 = lookup_triplanar_texture_2sz(pos, normal, tex1, tex1, tex1, tex_top) * color1;
	// interpolate between the two texture/color pairs using a random noise function
	return mix(texel0, texel1, procedural_eval(pos));
}
