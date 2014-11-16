uniform sampler2D tex0;

in vec2 tc;
in vec3 normal;

void main()
{
	vec4 texel   = texture(tex0, tc);
	vec3 ldir    = normalize(fg_LightSource[4].position.xyz); // only light 4
	vec3 diffuse = (gl_Color.rgb * fg_LightSource[4].diffuse.rgb);
	vec3 ambient = (gl_Color.rgb * (fg_LightSource[4].ambient.rgb + vec3(0.05, 0.0, 0.0))); // add a bit of red ambient
	fg_FragColor = texel * vec4((ambient + clamp(4.0*dot(normal, ldir), 0.0, 1.0)*diffuse), gl_Color.a); // high albedo
}
