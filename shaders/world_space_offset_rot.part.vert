uniform vec4 world_space_offset = vec4(0.0); // {x, y, z, rot_angle}

mat3 get_z_rotation(float angle) { // angle is negated?
	float s = sin(-angle);
	float c = cos(-angle);
	return mat3(vec3(c, -s, 0), vec3(s, c, 0), vec3(0, 0, 1));
}
void set_ws_pos_and_normal(out vec3 pos, out vec3 normal) {
#ifdef ENABLE_ROTATIONS
	mat3 rot_matrix = get_z_rotation(world_space_offset.w);
	pos    = rot_matrix * fg_Vertex.xyz + world_space_offset.xyz;
	normal = normalize(rot_matrix * fg_Normal);
#else // translate only
	pos    = fg_Vertex.xyz + world_space_offset.xyz;
	normal = normalize(fg_Normal);
#endif
}
