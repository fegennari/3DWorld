uniform vec3 center = vec3(0);

out vec3 dir;

void main() {
	gl_Position = fg_ModelViewProjectionMatrix * fg_Vertex;
	fg_Color_vf = fg_Color;
	dir = normalize(fg_Vertex.xyz - center); // world space direction from center
} 
