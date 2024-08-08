out vec2 tc;

void main() {
	tc          = fg_TexCoord;
	vec4 epos   = fg_ModelViewMatrix * fg_Vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
	fg_Color_vf = fg_Color;
} 
