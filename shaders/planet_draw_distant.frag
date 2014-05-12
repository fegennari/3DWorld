uniform vec3 light_scale = vec3(1.0);

varying vec4 epos;
varying vec3 normal;

// lighting, but no textures and no specular, and light2 indirect contribution is ignored
// Note: most of this is copied from planet_draw.frag
void main()
{
	float atten0 = light_scale[0] * calc_light_atten0(epos);
	vec3 ldir0   = normalize(fg_LightSource[0].position.xyz - epos.xyz);
	vec3 ambient = (fg_LightSource[0].ambient.rgb * atten0) + (fg_LightSource[1].ambient.rgb * light_scale[1]);
	vec3 diffuse = (fg_LightSource[0].diffuse.rgb * max(dot(normal, ldir0), 0.0) * atten0);
	fg_FragColor = gl_Color * vec4((ambient + diffuse), 1.0);
}

