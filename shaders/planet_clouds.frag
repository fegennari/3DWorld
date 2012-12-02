uniform vec2 light_scale = vec2(1,1);
uniform sampler2D tex0;
varying vec4 epos;
varying vec3 normal; // world space

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st); // FIXME: use procedural volumetric texture
	if (texel.a == 0.0) discard;
	vec3 norm  = normalize(normal); // renormalize
	vec4 color = gl_FrontMaterial.emission;
	color.rgb += light_scale[0]*add_pt_light_comp(norm, epos, 0).rgb;
	color.rgb += light_scale[1]*add_pt_light_comp(norm, epos, 1).rgb;
	gl_FragColor = apply_fog(color * texel);
}
