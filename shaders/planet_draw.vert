uniform mat4 fg_ViewMatrix;
uniform float water_val = 0.0;
uniform float lava_val  = 0.0;
uniform float obj_radius;

out vec3 normal, world_space_pos, vertex;
out vec2 tc;

void main()
{
	tc     = fg_TexCoord;
	vertex = obj_radius*fg_Vertex.xyz; // fg_Vertex is typically normalized to a length of 1.0
#ifdef PROCEDURAL_DETAIL
	// since we know the incoming vertices represent a perfect sphere around (0,0,0),
	// we can calculate the vertex from the normal, allowing us to skip sending the normals to the GPU
	normal       = normalize(fg_NormalMatrix * fg_Vertex.xyz);
	vec3 spos    = terrain_scale*fg_Vertex.xyz;
	vec3 npos    = spos + vec3(noise_offset);
	float hval   = eval_terrain_noise(npos, 8);
	float cutoff = max(water_val, lava_val);
	float height = pow((1.0 - cutoff), 0.1)*max(0.0, 1.8*(hval-0.7)); // can go outside the [0,1] range
	float omcinv = 1.0/max(0.01, (1.0 - cutoff)); // avoid div-by-zero
	vertex      *= 1.0 + 0.01*(omcinv*(max(cutoff, height) - cutoff) - 0.5);
#else
	normal       = normalize(fg_NormalMatrix * fg_Normal);
#endif
	vec4 epos    = fg_ModelViewMatrix * vec4(vertex, 1.0);
	world_space_pos = (inverse(fg_ViewMatrix) * epos).xyz;
	gl_Position  = fg_ProjectionMatrix * epos;
	gl_FrontColor= vec4(1.0); // always white - color will come from the texture
}

