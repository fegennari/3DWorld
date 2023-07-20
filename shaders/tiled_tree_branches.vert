out vec4 epos;
out vec3 normal; // eye space
out vec2 tc;
out vec3 ws_pos, ws_normal;

void main() {
	set_ws_pos_and_normal(ws_pos, ws_normal); // from world_space_offset_rot.part.vert
	tc          = fg_TexCoord;
	normal      = normalize(fg_NormalMatrix * fg_Normal);
	epos        = fg_ModelViewMatrix  * fg_Vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	fg_Color_vf = fg_Color;
}
