uniform float height = 1.0;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space (not normalized)
	vec4 vertex = gl_Vertex;
	// FIXME: add wind?
	vec4 epos   = gl_ModelViewMatrix  * vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	
	// calculate lighting
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp_pos(normal, epos, 0);
	if (enable_light1) color += add_light_comp_pos(normal, epos, 1);
	gl_FrontColor   = color;
	gl_FogFragCoord = length(epos.xyz);
} 
