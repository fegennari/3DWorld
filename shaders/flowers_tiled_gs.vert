uniform float x1, y1, dx_inv, dy_inv;
uniform sampler2D height_tex;

in float point_size;

out vec4 color_vs, vertex_vs;
out vec3 normal_vs;
out float size_vs;
out int dmin_vs;

int get_min_dim(vec3 v) {
	return ((abs(v.x) < abs(v.y)) ? ((abs(v.x) < abs(v.z)) ? 0:2) : ((abs(v.y) < abs(v.z)) ? 1:2));
}

void main()
{
	vec3 gwdelta = get_grass_wind_delta(fg_Vertex.xyz, 0.5); // wind_amt=0.5
	normal_vs    = normalize(normalize(fg_Normal) + gwdelta/height); // eye space, height comes from wind.part
	vertex_vs    = fg_Vertex + vec4(gwdelta, 0.0);
	vertex_vs.z += texture(height_tex, vec2((fg_Vertex.x - x1)*dx_inv, (fg_Vertex.y - y1)*dy_inv)).r;
	dmin_vs      = get_min_dim(fg_Normal);
	size_vs      = point_size;
	color_vs     = fg_Color;
}
