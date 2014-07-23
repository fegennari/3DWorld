uniform mat4 fg_ViewMatrix;
varying vec3 normal, world_space_pos, vertex;
varying vec2 tc;

uniform float water_val = 0.0;
uniform float lava_val  = 0.0;

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
	float cutoff = max(water_val, lava_val);
	float omcinv = 1.0/max(0.01, (1.0 - cutoff)); // avoid div-by-zero
	vertex2.xyz += 0.015*(omcinv*(max(cutoff, height) - cutoff) - 0.5)*obj_radius*fg_Normal;
#endif
	vertex       = vertex2.xyz;
	vec4 epos    = fg_ModelViewMatrix * vertex2;
	world_space_pos = (inverse(fg_ViewMatrix) * epos).xyz;
	gl_Position  = fg_ProjectionMatrix * epos;
	gl_FrontColor= vec4(1.0); // always white - color will come from the texture
}

