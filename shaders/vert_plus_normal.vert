out vec4 epos;
out vec3 eye_norm;

void main() {
	gl_Position = fg_ftransform();
	fg_Color_vf = fg_Color;
	epos        = fg_ModelViewMatrix * fg_Vertex;
	eye_norm    = normalize(fg_NormalMatrix * fg_Normal); // eye space
}
