uniform float x_scene_size, y_scene_size, czmin, czmax; // scene bounds (world space)
uniform float step_delta;
uniform sampler2D tex0;
uniform sampler3D smoke_tex;
uniform sampler2D dlight_tex;
uniform usampler2D dlelm_tex, dlgb_tex;
uniform float min_alpha = 0.0;

// clipped eye position, clipped vertex position, starting vertex position
varying vec3 eye, vpos, spos, dlpos, normal; // world space

const float SMOKE_SCALE = 0.25;


float get_intensity_at(in vec3 pos, in vec4 lpos_r) {
	float radius = lpos_r.w;
	if (radius == 0.0) return 1.0; // no falloff
	if (abs(pos.z - lpos_r.z) > radius) return 0.0; // fast test
	float dist   = distance(pos, lpos_r.xyz);
	float rscale = max(0.0, (radius - dist))/radius;
	float mag = 0.25 + 0.75*max(0.0, dot(normal, normalize(lpos_r.xyz - pos))); // ambient + diffuse
	// FIXME: add specular
	return mag*rscale*rscale; // quadratic 1/r^2 attenuation
}

vec3 add_dlights(in vec3 pos, in vec3 off, in vec3 scale) {
	vec3 color  = vec3(0,0,0);
	uint gb_ix  = texture2D(dlgb_tex, pos.xy).r; // get grid bag element index range (uint32)
	uint st_ix  = (gb_ix & 0xFFFFU);
	uint end_ix = ((gb_ix >> 16U) & 0xFFFFU);
	const uint elem_tex_sz = 256;  // must agree with value in C++ code, or can use textureSize()
	const uint max_dlights = 1024; // must agree with value in C++ code, or can use textureSize()
	
	for (uint i = st_ix; i < end_ix; ++i) { // iterate over grid bag elements
		uint dl_ix  = texelFetch(dlelm_tex, ivec2((i%elem_tex_sz), (i/elem_tex_sz)), 0).r; // get dynamic light index (uint16)
		vec4 lpos_r = texelFetch(dlight_tex, ivec2(0, dl_ix), 0); // light center, radius
		vec4 lcolor = texelFetch(dlight_tex, ivec2(1, dl_ix), 0); // light color
		lpos_r.xyz *= scale; // convert from [0,1] back into world space
		lpos_r.xyz += off;
		lpos_r.w   *= x_scene_size;
		color += lcolor.rgb * (lcolor.a * get_intensity_at(dlpos, lpos_r));
		color  = clamp(color, 0.0, 1.0);
		if (color.rgb == vec3(1,1,1)) break; // saturated
	}
	return color;
}


// Note: This may seem like it can go into the vertex shader as well,
//       but we don't have the tex0 value there and can't determine the full init color
void main()
{
	vec3 off   = vec3(-x_scene_size, -y_scene_size, czmin);
	vec3 scale = vec3(2.0*x_scene_size, 2.0*y_scene_size, (czmax - czmin));
	vec3 lit_color  = gl_Color.rgb; // base color (with some lighting)
	
	if (do_lighting) {
		vec3 sp  = clamp((spos  - off)/scale, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec3 dlp = clamp((dlpos - off)/scale, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec3 indir_light = texture3D(smoke_tex, sp.zxy).rgb; // add indir light color from texture
		lit_color += gl_FrontMaterial.diffuse.rgb * indir_light; // indirect lighting
		lit_color += add_dlights(dlp, off, scale); // dynamic lighting
	}
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec4 color = vec4((texel.rgb * lit_color), (texel.a * gl_Color.a));
	if (keep_alpha && color.a <= min_alpha) discard;
	
	if (eye == vpos) {
		if (color.a <= min_alpha) discard;
		gl_FragColor = color;
		return;
	}
	vec3 pos   = (vpos - off)/scale;
	vec3 dir   = eye - vpos;
	vec3 delta = normalize(dir)*step_delta/scale;
	float nsteps = length(dir)/step_delta;
	int num_steps = 1 + min(1000, int(nsteps)); // round up
	float step_weight = fract(nsteps);
	
	// smoke volume iteration using 3D texture, pos to eye
	for (int i = 0; i < num_steps; ++i) {
		vec3 p = clamp(pos, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec4 tex_val = texture3D(smoke_tex, p.zxy); // rgba = {color.rgb, smoke}
		float smoke = SMOKE_SCALE*tex_val.a*step_weight;
		vec3 rgb_comp = (tex_val.rgb * gl_Fog.color.rgb);
		color = ((!keep_alpha && color.a == 0.0) ? vec4(rgb_comp, smoke) :
					mix(color, vec4(rgb_comp, (keep_alpha ? color.a : 1.0)), smoke));
		pos += delta*step_weight;
		step_weight = 1.0;
	}
	if (color.a <= min_alpha) discard;
	gl_FragColor = color;
}
