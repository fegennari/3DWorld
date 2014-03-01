uniform float half_dxy;
uniform float indir_vert_offset = 0.25;
// Note: these next two come from dynamic_lighting.part, which must always be included first
//uniform vec3 scene_llc, scene_scale; // scene bounds (world space)
uniform sampler3D smoke_and_indir_tex;
uniform vec3 const_indir_color = vec3(0.0);

varying vec3 vpos, normal; // world space

void add_indir_lighting(inout vec3 lit_color) {
#ifdef USE_LIGHT_COLORS
	vec3 diffuse = gl_Color.rgb;
#else // use light material properties
	vec3 diffuse = gl_FrontMaterial.diffuse.rgb;
#endif
	lit_color += diffuse * const_indir_color; // add constant indir
	
	if (indir_lighting) {
		vec3 spos = vpos + (indir_vert_offset*half_dxy)*normal; // move slightly away from the vertex
		spos = clamp((spos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec3 indir_light = texture3D(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
		//indir_light = pow(indir_light, vec3(0.45)); // gamma correction
		lit_color += diffuse * indir_light; // indirect lighting
	}
}
