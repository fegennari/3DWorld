layout(quads, equal_spacing, cw) in;

uniform float wave_time;

in vec3 vertex_ES[];
in vec2 tc_ES[], tc2_ES[];

out vec4 epos, proj_pos;
out vec2 tc, tc2;
out vec4 fg_Color_vf;
out float water_zval;

vec2 interpolate2D(vec2 v0, vec2 v1, vec2 v2, vec2 v3) {
	return mix(mix(v0, v1, gl_TessCoord.x), mix(v3, v2, gl_TessCoord.x), gl_TessCoord.y);
}
vec3 interpolate3D(vec3 v0, vec3 v1, vec3 v2, vec3 v3) {
	return mix(mix(v0, v1, gl_TessCoord.x), mix(v3, v2, gl_TessCoord.x), gl_TessCoord.y);
}

void main() {
	// Interpolate the attributes of the output vertex using the barycentric coordinates
	fg_Color_vf = vec4(1.0);
	tc          = interpolate2D( tc_ES[0],  tc_ES[1],  tc_ES[2],  tc_ES[3]);
	tc2         = interpolate2D(tc2_ES[0], tc2_ES[1], tc2_ES[2], tc2_ES[3]);
	vec3 vertex = interpolate3D(vertex_ES[0], vertex_ES[1], vertex_ES[2], vertex_ES[3]);
	vec2 val    = 6.0*3.14159*tc.st + vec2(0.060*wave_time, 0.063*wave_time); // 6*PI
	vertex.z   += 0.025*sin(val.x)*sin(val.y);
	water_zval  = vertex.z;
	epos        = fg_ModelViewMatrix * vec4(vertex, 1.0);
	proj_pos    = fg_ProjectionMatrix * epos;
	gl_Position = proj_pos;
}
