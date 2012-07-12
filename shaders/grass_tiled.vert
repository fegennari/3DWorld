uniform float height = 1.0;
uniform float x1, y1, x2, y2, zmin, zmax;
uniform mat4 world_space_mvm;
uniform sampler2D height_tex;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space (not normalized)

	vec4 vertex = gl_Vertex;
	vec4 vpos   = (inverse(world_space_mvm) * (gl_ModelViewMatrix * vertex)); // world space
	vertex.z   += zmin + (zmax - zmin)*texture2D(height_tex, vec2((vpos.x - x1)/(x2 - x1), (vpos.y - y1)/(y2 - y1))).r;
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
