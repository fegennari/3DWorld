uniform vec2 light_scale = vec2(1,1);
varying vec4 epos;
varying vec3 normal, world_space_pos;

void main()
{
	vec3 norm  = normalize(normal); // renormalize
	vec4 color = gl_FrontMaterial.emission;
	color.rgb += light_scale[0]*add_pt_light_comp(norm, epos, 0).rgb; // ambient, diffuse, and specular
	color.rgb += light_scale[1]*(gl_Color * gl_LightSource[1].ambient).rgb; // ambient only

	// FIXME: color *= some atmospheric stuff
	color *= vec4(0.5, 0.5, 1.0, 1.0);

	gl_FragColor = apply_fog(color);
}
