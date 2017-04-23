layout(triangles, equal_spacing, cw) in;

in vec3 vertex_ES[];
in vec2 tc_ES[];
in vec4 fg_Color_vf_ES[];
#ifdef INCLUDE_NORMALS
in vec3 normal_ES[], eye_norm_ES[];
#endif

out vec2 tc;
out vec4 epos;
out vec4 fg_Color_vf;
#ifdef INCLUDE_NORMALS
out vec3 vertex_from_vs, normal, eye_norm;
#endif

vec2 interpolate2D(vec2 v0, vec2 v1, vec2 v2) {
	return gl_TessCoord.x * v0 + gl_TessCoord.y * v1 + gl_TessCoord.z * v2;
}
vec3 interpolate3D(vec3 v0, vec3 v1, vec3 v2) {
	return gl_TessCoord.x * v0 + gl_TessCoord.y * v1 + gl_TessCoord.z * v2;
}
vec4 interpolate4D(vec4 v0, vec4 v1, vec4 v2) {
	return gl_TessCoord.x * v0 + gl_TessCoord.y * v1 + gl_TessCoord.z * v2;
}

void main() {
	// Interpolate the attributes of the output vertex; vertex 2 is the tip of the grass blade
	tc          = interpolate2D(    tc_ES[0],     tc_ES[1],     tc_ES[2]);
	vec3 vertex = interpolate3D(vertex_ES[0], vertex_ES[1], vertex_ES[2]);
	epos        = fg_ModelViewMatrix * vec4(vertex, 1.0); // use epos before modifying vertex
	vec3 delta  = vertex - 0.5*(vertex_ES[0] + vertex_ES[1]); // delta from base of grass
	vertex     -= delta;
	float dist  = length(delta);
	vec3 dir    = normalize(cross(delta, (vertex_ES[0] - vertex_ES[1]))); // direction of bend motion
	float mval  = clamp((2.0 - length(epos.xyz)), 0.0, 1.0)*(0.6 - pow(tc[0], 0.4)); // more curvature at the tip and near the camera
	delta       = mix(delta, dist*dir, mval);
	delta      *= dist/length(delta); // renormalize to original length
	vertex     += delta;
	gl_Position = fg_ModelViewProjectionMatrix * vec4(vertex, 1.0);
	fg_Color_vf = interpolate4D(fg_Color_vf_ES[0], fg_Color_vf_ES[1], fg_Color_vf_ES[2]); // could also use tip color

#ifdef INCLUDE_NORMALS
	vertex_from_vs = vertex;
	normal   = interpolate3D(normal_ES  [0], normal_ES  [1], normal_ES  [2]); // FIXME: update based on new vertex?
	eye_norm = interpolate3D(eye_norm_ES[0], eye_norm_ES[1], eye_norm_ES[2]); // FIXME: calculate from normal?
#endif
}
