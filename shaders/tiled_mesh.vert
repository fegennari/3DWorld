uniform float x1, y1, dx_inv, dy_inv;
uniform float htex_scale = 1.0;
uniform sampler2D height_tex;

varying vec4 vertex;
varying vec2 tc3;

void main()
{
	vertex    = fg_Vertex;
	vertex.z += htex_scale*texture2D(height_tex, vec2((vertex.x - x1)*dx_inv, (vertex.y - y1)*dy_inv)).r;
	vec4 epos = gl_ModelViewMatrix * vertex;
	tc  = vec2(dot(epos, gl_EyePlaneS[0]), dot(epos, gl_EyePlaneT[0]));
	tc3 = vec2(dot(epos, gl_EyePlaneS[2]), dot(epos, gl_EyePlaneT[2]));
	gl_Position = gl_ModelViewProjectionMatrix * vertex;
	set_fog_coord(vertex);
} 
