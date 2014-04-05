varying vec4 vertex;
#ifdef USE_HEIGHT_TEX
uniform float x1, y1, dx_inv, dy_inv;
uniform float htex_scale = 1.0;
uniform sampler2D height_tex;
#endif

void main()
{
	setup_texgen0();
	setup_texgen1();
	setup_texgen2();
	vertex = fg_Vertex;
#ifdef USE_HEIGHT_TEX
	vertex.z += htex_scale*texture2D(height_tex, vec2((vertex.x - x1)*dx_inv, (vertex.y - y1)*dy_inv)).r;
	gl_Position = gl_ModelViewProjectionMatrix * vertex;
#else
	gl_Position = fg_ftransform();
#endif
	set_fog_coord(vertex);
} 
