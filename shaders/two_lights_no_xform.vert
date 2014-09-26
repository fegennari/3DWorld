uniform float color_scale = 1.0; // for snow
uniform vec4 emission = vec4(0,0,0,1);

out vec4 color; // to GS

void main()
{
	gl_Position = fg_Vertex;
	vec3 normal = normalize(fg_NormalMatrix * fg_Normal); // eye space
	color       = emission;
	if (enable_light0) color += add_light_comp0(normal);
	if (enable_light1) color += add_light_comp1(normal);
	color.rgb *= vec3(color_scale);
} 
