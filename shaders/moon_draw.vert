out vec2 tc;

void main()
{
	tc            = fg_TexCoord;
	gl_Position   = fg_ftransform();
	vec3 normal   = normalize(fg_NormalMatrix * fg_Normal); // eye space
	vec3 ldir     = normalize(fg_LightSource[4].position.xyz); // only light 4
	vec3 diffuse  = (gl_Color * fg_LightSource[4].diffuse).rgb;
	vec3 ambient  = (gl_Color * fg_LightSource[4].ambient).rgb;
	gl_FrontColor = vec4((ambient + clamp(4.0*dot(normal, ldir), 0.0, 1.0)*diffuse), gl_Color.a); // high albedo
}
