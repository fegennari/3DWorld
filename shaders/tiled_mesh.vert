uniform float x1, y1, dx_inv, dy_inv;
uniform float htex_scale = 1.0;
uniform vec3 tc_xlate    = vec3(0);
uniform vec4 texgen_s, texgen_t, texgen2_s, texgen2_t;
uniform sampler2D height_tex;

varying vec4 vertex;
varying vec2 tc, tc2;

void main()
{
	vec4 vert = fg_Vertex + vec4(tc_xlate, 0.0);
	tc        = vec2(dot(fg_Vertex, texgen_s ), dot(fg_Vertex, texgen_t ));
	tc2       = vec2(dot(vert, texgen2_s), dot(vert, texgen2_t));
	vertex    = fg_Vertex;
	vertex.z += htex_scale*texture2D(height_tex, vec2((vertex.x - x1)*dx_inv, (vertex.y - y1)*dy_inv)).r;
	gl_Position = fg_ModelViewProjectionMatrix * vertex;
} 
