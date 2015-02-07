uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform mat4 fg_ViewMatrix;
uniform float tex_scale_s  = 1.0;
uniform float tex_scale_t  = 1.0;
uniform vec3 world_space_offset = vec3(0.0);
uniform vec3 sun_pos; // used for dynamic smoke shadows line clipping

in vec4 tex0_s, tex0_t;

out vec3 vpos, normal; // world space
out vec4 epos;
out vec3 eye_norm;
// tc comes from texture_gen.part.vert

void main()
{
	if (use_texgen == 1) {
		setup_texgen_st();
	}
	else if (use_texgen == 2) {
		tc = vec2(dot(fg_Vertex, tex0_s), dot(fg_Vertex, tex0_t));
	}
	else if (use_texgen == 3) {
		set_tc0_from_vert_id();
	}
	else {
		tc = fg_TexCoord * vec2(tex_scale_s, tex_scale_t);
	}
	fg_Color_vf = fg_Color;
	epos        = fg_ModelViewMatrix * fg_Vertex;
	gl_Position = fg_ProjectionMatrix * epos;

	if (use_fg_ViewMatrix) {
		eye_norm = normalize(fg_NormalMatrix * fg_Normal);
		normal   = normalize((transpose(fg_ViewMatrix) * vec4(eye_norm, 1)).xyz);
		vpos     = (inverse(fg_ViewMatrix) * epos).xyz; // world space
	}
	else {
		eye_norm = normalize(mat3(fg_ModelViewMatrix) * fg_Normal); // Note: avoids the fg_NormalMatrix upload
		vpos     = fg_Vertex.xyz + world_space_offset;
		normal   = normalize(fg_Normal);
	}
#ifdef USE_BUMP_MAP
	setup_tbn();
#endif
} 
