uniform float x1, y1, dx_inv, dy_inv;
uniform float htex_scale = 1.0;
uniform vec4 texgen2_s, texgen2_t;
uniform sampler2D height_tex;

varying vec4 vertex;
varying vec2 tc2;

void main()
{
	setup_texgen_st();
	tc2       = vec2(dot(fg_Vertex, texgen2_s), dot(fg_Vertex, texgen2_t));
	vertex    = fg_Vertex;
	vertex.z += htex_scale*texture2D(height_tex, vec2((vertex.x - x1)*dx_inv, (vertex.y - y1)*dy_inv)).r;
	vec4 epos = gl_ModelViewMatrix * vertex;
	gl_Position = gl_ModelViewProjectionMatrix * vertex;
	set_fog_coord(vertex);
} 
