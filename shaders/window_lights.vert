out vec2 tc;

void main() {
	tc          = fg_TexCoord;
	gl_Position = fg_ProjectionMatrix * (fg_ModelViewMatrix * fg_Vertex); // Note: must multiply this way so that values agree with first pass shader
	fg_Color_vf = fg_Color;
} 
