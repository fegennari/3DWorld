uniform vec4 color_modulate = vec4(1.0);
uniform float animation_time = 0.0;

out vec2 tc;

const float PI = 3.14159;

mat3 do_rotation(vec3 a, float angle) { // Note: axis does not need to be normalized
	float s = sin(angle);
	float c = cos(angle);
	float oc = 1.0 - c;
	float sx = s * a.x;
	float sy = s * a.y;
	float sz = s * a.z;
	float ocx = oc * a.x;
	float ocy = oc * a.y;
	float ocz = oc * a.z;
	float ocxx = ocx * a.x;
	float ocxy = ocx * a.y;
	float ocxz = ocx * a.z;
	float ocyy = ocy * a.y;
	float ocyz = ocy * a.z;
	float oczz = ocz * a.z;
	return mat3(vec3(ocxx+c, ocxy-sz, ocxz+sy), vec3(ocxy+sz, ocyy+c, ocyz-sx), vec3(ocxz-sy, ocyz+sx, oczz+c));
}

void rotate_about(inout vec3 vertex, in vec3 about, in vec3 axis, in float angle) {
	vertex -= about; // rotate about this point
	vertex *= do_rotation(axis, angle);
	vertex += about;
}

void main() { // not lit, no normal
	tc          = fg_TexCoord;
	vec3 axis   = vec3(0.7, 0.0, 0.7); // vertical axis of rotation, fish model is oriented 45 degrees upward
	vec3 center = vec3(3.15, -1.36, 0.24); // center of the model, determined experimentally (close to bcube center)
	vec4 vertex = fg_Vertex;
	float val   = 0.8*sin(40.0*animation_time);
	float rot   = clamp((vertex.x - 0.3), 0.0, 1.0)*val;
	rotate_about(vertex.xyz, center, axis, rot);
	gl_Position = fg_ModelViewProjectionMatrix * vertex;
	fg_Color_vf = fg_Color * color_modulate;
}
