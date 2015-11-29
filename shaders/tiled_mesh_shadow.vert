uniform float x1, y1, dx_inv, dy_inv;
uniform float delta_z = 0.0;
uniform sampler2D height_tex;

void main() {
	vec4 vertex = fg_Vertex;
	vertex.z   += texture(height_tex, vec2((vertex.x - x1)*dx_inv, (vertex.y - y1)*dy_inv)).r + delta_z;
	gl_Position = fg_ModelViewProjectionMatrix * vertex;
} 
