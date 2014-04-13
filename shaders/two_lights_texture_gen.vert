uniform vec4 emission = vec4(0,0,0,1);

void main()
{
	setup_texgen_st();
	vec4 epos   = fg_ModelViewMatrix * fg_Vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	vec3 normal = normalize(fg_NormalMatrix * fg_Normal); // eye space
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
	vec4 color = emission;
	if (enable_light0) color.rgb += add_light_comp0(normal).rgb;
	if (enable_light1) color.rgb += add_light_comp1(normal).rgb;
	gl_FrontColor = color;
} 
