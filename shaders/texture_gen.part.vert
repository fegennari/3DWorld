void setup_texgen0() {
	vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
	gl_TexCoord[0] = vec4(dot(ecPosition, gl_EyePlaneS[0]), dot(ecPosition, gl_EyePlaneT[0]), 0.0, 1.0);
}

void setup_texgen1() {
	vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
	gl_TexCoord[1] = vec4(dot(ecPosition, gl_EyePlaneS[1]), dot(ecPosition, gl_EyePlaneT[1]), 0.0, 1.0);
}

void setup_texgen2() {
	vec4 ecPosition = gl_ModelViewMatrix * gl_Vertex;
	gl_TexCoord[2] = vec4(dot(ecPosition, gl_EyePlaneS[2]), dot(ecPosition, gl_EyePlaneT[2]), 0.0, 1.0);
}
