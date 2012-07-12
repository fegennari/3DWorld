uniform float height = 1.0;
uniform float x1, y1, x2, y2, zmin, zmax, translate_x, translate_y;
uniform sampler2D height_tex;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space (not normalized)

	vec4 vertex = gl_Vertex;
	vertex.xy  += vec2(translate_x, translate_y);
	vertex.z   += zmin + (zmax - zmin)*texture2D(height_tex, vec2((vertex.x - x1)/(x2 - x1), (vertex.y - y1)/(y2 - y1))).r;
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
