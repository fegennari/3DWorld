uniform sampler3D smoke_and_indir_tex;
uniform vec3 const_indir_color = vec3(0,0,0);

varying vec3 spos; // world space

void add_indir_lighting(inout vec3 lit_color) {
	lit_color += gl_FrontMaterial.diffuse.rgb * const_indir_color; // add constant indir
	
	if (indir_lighting) {
		vec3 indir_light = texture3D(smoke_and_indir_tex, spos.zxy).rgb; // add indir light color from texture
		lit_color += gl_FrontMaterial.diffuse.rgb * indir_light; // indirect lighting
	}
}
