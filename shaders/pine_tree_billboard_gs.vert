out vec4 vertex_vs;
out vec4 color_vs;
out vec2 xz_sz_vs;

void main() {
	vertex_vs = fg_Vertex;
	xz_sz_vs  = fg_Normal.xy; // Note: could use vec2 xy_sz attribute

	// lighting
	//const float normal_z = 0.816; // high detail tree polygon normal
	const float normal_z = 0.95; // blends better in practice
	vec3 normal = normalize(normal_z*fg_NormalMatrix[2] + vec3(0.0, 0.0, sqrt(1.0 - normal_z*normal_z))) / ambient_scale;
	vec3 color  = vec3(0.0);
	if (enable_light0) {color += add_light_comp0(normal).rgb;}
	if (enable_light1) {color += add_light_comp1(normal).rgb;}
	color_vs = vec4(color, fg_Color.a);
}
