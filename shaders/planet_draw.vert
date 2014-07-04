uniform mat4 fg_ViewMatrix;
varying vec3 normal, world_space_pos, vertex;
varying vec2 tc;

#ifdef PROCEDURAL_DETAIL
uniform float terrain_scale, obj_radius, noise_offset;
uniform sampler3D cloud_noise_tex;

float eval_terrain_noise(in vec3 npos, const int num_octaves) {
	float val  = 0.0;
	float mag  = 1.0;
	float freq = 0.5; // lower freq for ridged noise

	for (int i = 0; i < num_octaves; ++i) { // similar to gen_cloud_alpha_time()
		float v = texture3D(cloud_noise_tex, freq*npos).r;
		v = 2.0*v - 1.0; // map [0,1] range to [-1,1]
		v = max(0.0, (0.75 - abs(v))); // ridged noise
		val  += v*mag;
		freq *= 1.92;
		mag  *= 0.5;
	}
	return val;
}
#endif

void main()
{
	tc           = fg_TexCoord;
	normal       = normalize(fg_NormalMatrix * fg_Normal);
	vec4 vertex2 = fg_Vertex;
#ifdef PROCEDURAL_DETAIL
	vec3 spos    = fg_Vertex.xyz*(terrain_scale/obj_radius);
	vec3 npos    = spos + vec3(noise_offset);
	float hval   = eval_terrain_noise(npos, 8);
	float height = max(0.0, 1.8*(hval-0.7)); // can go outside the [0,1] range
	vertex2.xyz += (0.02*height - 0.01)*obj_radius*fg_Normal;
#endif
	vertex       = vertex2.xyz;
	vec4 epos    = fg_ModelViewMatrix * vertex2;
	world_space_pos = (inverse(fg_ViewMatrix) * epos).xyz;
	gl_Position  = fg_ProjectionMatrix * epos;
	gl_FrontColor= vec4(1.0); // always white - color will come from the texture
}

