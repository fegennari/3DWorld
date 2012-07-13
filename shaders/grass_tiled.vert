uniform float height = 1.0;
uniform float dist_const = 10.0;
uniform float dist_slope = 0.5;
uniform float x1, y1, x2, y2, wx2, wy2, zmin, zmax, translate_x, translate_y;
uniform sampler2D height_tex, weight_tex, noise_tex;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	vec4 vertex = gl_Vertex;
	vertex.xy  += vec2(translate_x, translate_y);
	float grass_weight = texture2D(weight_tex, vec2(vertex.x/wx2, vertex.y/wy2)).b;
	float noise_weight = texture2D(noise_tex,  vec2(10.0*gl_Color.r, 10.0*gl_Color.g)).r; // "hash" the color

	vertex.z   += zmin + (zmax - zmin)*texture2D(height_tex, vec2((vertex.x - x1)/(x2 - x1), (vertex.y - y1)/(y2 - y1))).r;
	// FIXME: add wind?
	vec4 epos   = gl_ModelViewMatrix  * vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	gl_FogFragCoord = length(epos.xyz);
	grass_weight   *= 1.0 - clamp(dist_slope*(gl_FogFragCoord - dist_const), 0.0, 1.0); // decrease weight far away from camera
	
	// calculate lighting
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space (not normalized)
	vec4 color  = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp_pos(normal, epos, 0);
	if (enable_light1) color += add_light_comp_pos(normal, epos, 1);
	color.a = ((grass_weight < noise_weight) ? 0.0 : color.a); // skip some grass blades by making them transparent
	gl_FrontColor  = color;
} 
