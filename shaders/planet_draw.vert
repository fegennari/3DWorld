uniform mat4 fg_ViewMatrix;
uniform float water_val = 0.0;
uniform float lava_val  = 0.0;
uniform float obj_radius;

out vec3 normal, world_space_pos, vertex;
out vec2 tc;

#ifdef PROCEDURAL_DETAIL
float get_delta_h(in vec3 pos) {
	float hval   = eval_terrain_noise(pos, 8);
	float cutoff = max(water_val, lava_val);
	float height = pow((1.0 - cutoff), 0.1)*max(0.0, 1.8*(hval-0.7)); // can go outside the [0,1] range
	float omcinv = 1.0/max(0.01, (1.0 - cutoff)); // avoid div-by-zero
	return omcinv*(max(cutoff, height) - cutoff) - 0.5;
}
#endif

void main()
{
	tc     = fg_TexCoord;
	vertex = obj_radius*fg_Vertex.xyz; // fg_Vertex is typically normalized to a length of 1.0
#ifdef PROCEDURAL_DETAIL
	// since we know the incoming vertices represent a perfect sphere around (0,0,0),
	// we can calculate the normal from the vertex, allowing us to skip sending the normals to the GPU
	vec3 npos  = terrain_scale*fg_Vertex.xyz + vec3(noise_offset);
	float hval = get_delta_h(npos);
	vertex    *= 1.0 + 0.01*hval;
	normal     = normalize(fg_NormalMatrix * vertex);
#if 1 // use derivative at this vertex to get a better quality (but slower) vertex normal
	float delta = 0.001;
	float hdx   = hval - get_delta_h(npos + vec3(delta, 0.0, 0.0));
	float hdy   = hval - get_delta_h(npos + vec3(0.0, delta, 0.0));
	float hdz   = hval - get_delta_h(npos + vec3(0.0, 0.0, delta));
	normal      = normalize(normal + 0.5*(fg_NormalMatrix * vec3(hdx, hdy, hdz)));
#endif

#else
	normal       = normalize(fg_NormalMatrix * fg_Normal);
#endif
	vec4 epos    = fg_ModelViewMatrix * vec4(vertex, 1.0);
	world_space_pos = (inverse(fg_ViewMatrix) * epos).xyz;
	gl_Position  = fg_ProjectionMatrix * epos;
	gl_FrontColor= vec4(1.0); // always white - color will come from the texture
}

