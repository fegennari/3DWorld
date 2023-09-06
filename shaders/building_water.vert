out vec3 vpos; // world space
out vec4 epos, proj_pos;

void main() {
	epos        = fg_ModelViewMatrix * fg_Vertex;
	proj_pos    = fg_ProjectionMatrix * epos;
	vpos        = fg_Vertex.xyz;
	gl_Position = proj_pos;
} 
