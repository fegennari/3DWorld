layout(quads, equal_spacing, cw) in;

uniform float wave_time = 0.0;
uniform float normal_z  = 1.0;
uniform vec3 camera_pos;

in vec3 vertex_ES[];
in vec2 tc_ES[], tc2_ES[];

out vec4 epos, proj_pos;
out vec2 tc, tc2;
out vec4 fg_Color_vf;
out float water_zval;
out vec3 normal;

vec2 interpolate2D(vec2 v0, vec2 v1, vec2 v2, vec2 v3) {
	return mix(mix(v0, v1, gl_TessCoord.x), mix(v3, v2, gl_TessCoord.x), gl_TessCoord.y);
}
vec3 interpolate3D(vec3 v0, vec3 v1, vec3 v2, vec3 v3) {
	return mix(mix(v0, v1, gl_TessCoord.x), mix(v3, v2, gl_TessCoord.x), gl_TessCoord.y);
}

float get_delta_z(in vec2 v) {
	vec2 val = 4.0*3.14159*v.st; // 6*PI
	float dz = 0.0;
	dz += sin( 1.0*val.x + 0.0*val.y + 0.060*wave_time + 0.0)*sin( 0.0*val.x + 1.0*val.y + 0.063*wave_time + 0.0);
	dz += sin( 0.3*val.x - 0.9*val.y - 0.053*wave_time + 0.5)*sin(-0.7*val.x - 0.4*val.y - 0.050*wave_time + 0.7);
	dz += sin(-0.7*val.x + 0.4*val.y + 0.048*wave_time + 0.7)*sin(-0.3*val.x + 0.6*val.y - 0.051*wave_time + 0.3);
	dz += sin(-0.2*val.x - 0.6*val.y - 0.065*wave_time + 0.2)*sin( 0.9*val.x - 0.3*val.y + 0.059*wave_time + 0.8);
	return 0.02*dz;
}

void main() {
	// Interpolate the attributes of the output vertex
	fg_Color_vf = vec4(1.0);
	tc          = interpolate2D(    tc_ES[0],     tc_ES[1],     tc_ES[2],     tc_ES[3]);
	tc2         = interpolate2D(   tc2_ES[0],    tc2_ES[1],    tc2_ES[2],    tc2_ES[3]);
	vec3 vertex = interpolate3D(vertex_ES[0], vertex_ES[1], vertex_ES[2], vertex_ES[3]);
	float dz    = get_delta_z(tc);
	vertex.z   += dz;

	float intensity = 200.0 * clamp((40.0/distance(vertex.xyz, camera_pos) - 1.0), 0.0, 1.0);
	float dzx   = intensity*(get_delta_z(tc + vec2(0.001, 0.0)) - dz);
	float dzy   = intensity*(get_delta_z(tc + vec2(0.0, 0.001)) - dz);
	normal      = fg_NormalMatrix * normalize(vec3(-dzx, -dzy, normal_z));

	water_zval  = vertex.z;
	epos        = fg_ModelViewMatrix * vec4(vertex, 1.0);
	proj_pos    = fg_ProjectionMatrix * epos;
	gl_Position = proj_pos;
}
