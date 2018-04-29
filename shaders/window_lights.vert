out vec2 tc;
out vec4 epos;

void main() {
	tc          = fg_TexCoord;
	epos        = fg_ModelViewMatrix * fg_Vertex;
	gl_Position = fg_ProjectionMatrix * epos; // Note: must multiply this way so that values agree with first pass shader
	fg_Color_vf = fg_Color;
} 
