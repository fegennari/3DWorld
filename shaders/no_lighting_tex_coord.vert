uniform vec3 xlate = vec3(0);
uniform vec3 scale = vec3(1);

out vec2 tc;

void main()
{
	tc          = fg_TexCoord;
	gl_Position = fg_ModelViewProjectionMatrix * (vec4(xlate, 0.0) + (vec4(scale, 1.0) * fg_Vertex));
	fg_Color_vf = fg_Color;
} 
