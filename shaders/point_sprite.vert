uniform float point_scale = 1.0; // or point size in pixels, for constant size
in float point_size; // world space

#ifdef ENABLE_LIGHTING
uniform float half_dxy;
//uniform vec3 scene_llc, scene_scale; // scene bounds (world space), comes from dynamic_lighting.part
uniform sampler3D smoke_and_indir_tex;
uniform vec3 const_indir_color = vec3(0.0);

void add_indir_lighting(inout vec3 lit_color) {
	vec3 indir_color = const_indir_color; // add constant indir

	if (indir_lighting) {
		vec3 spos    = fg_Vertex.xyz + (0.25*half_dxy)*fg_Normal; // move slightly away from the vertex
		spos         = clamp((spos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
		indir_color += texture(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
	}
	lit_color += fg_Color.rgb * indir_color;
}
#endif


void main()
{
	vec4 epos    = fg_ModelViewMatrix * fg_Vertex;
	gl_Position  = fg_ProjectionMatrix * epos;
#ifdef CONSTANT_PT_SIZE
	gl_PointSize = point_scale;
#else
	gl_PointSize = clamp(point_scale*point_size/length(epos.xyz), 1.0, 64.0);
#endif

#ifdef ENABLE_LIGHTING
	vec3 normal = normalize(fg_NormalMatrix * fg_Normal); // eye space
	vec3 color  = vec3(0.0);
	if (enable_light0) {color += add_light_comp_pos_smap_light0(normal, epos).rgb;}
	if (enable_light1) {color += add_light_comp_pos_smap_light0(normal, epos).rgb;}
	add_indir_lighting(color);
	if (enable_dlights) {add_dlights(color.rgb, fg_Vertex.xyz, normalize(fg_Normal), vec3(1.0));} // dynamic lighting
	fg_Color_vf = vec4(color, fg_Color.a);
	//fg_Color_vf = apply_fog_epos(fg_Color_vf, epos);
#else
	fg_Color_vf = fg_Color;
#endif
} 
