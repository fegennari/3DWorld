uniform vec3 world_space_offset = vec3(0.0);

out vec4 epos;
out vec3 normal; // eye space
out vec2 tc;
out vec3 ws_pos;
out vec3 ws_normal;

void main() {
	tc          = fg_TexCoord;
	ws_pos      = fg_Vertex.xyz + world_space_offset;
	ws_normal   = normalize(fg_Normal);
	normal      = normalize(fg_NormalMatrix * fg_Normal);
	epos        = fg_ModelViewMatrix  * fg_Vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	fg_Color_vf = fg_Color;
}
