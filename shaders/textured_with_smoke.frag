uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta;
uniform sampler2D tex0;
uniform sampler3D smoke_tex;
uniform float min_alpha = 0.0;

uniform float light_atten = 0.0, refract_ix = 1.0;
uniform float cube_bb[6];

// clipped eye position, clipped vertex position, starting vertex position
varying vec3 eye, vpos, spos, normal, lpos0, vposl; // world space
varying vec3 eye_norm;
varying vec4 epos;

const float SMOKE_SCALE = 0.25;

// Note: dynamic point lights use reflection vector for specular, and specular doesn't move when the eye rotates
//       global directional lights use half vector for specular, which seems to be const per pixel, and specular doesn't move when the eye translates
#define ADD_LIGHT(i) lit_color += add_light_comp(n, i).rgb

vec3 add_light(in int ix, in vec3 off, in vec3 scale, in vec3 n, in vec3 source, in vec3 dest) {

	float lscale = 1.0;

	if (smoke_enabled && dynamic_smoke_shadows && source != dest) {
		vec3 dir      = dest - source;
		vec3 pos      = (source - off)/scale;
		vec3 delta    = normalize(dir)*step_delta/scale;
		float nsteps  = length(dir)/step_delta;
		int num_steps = 1 + min(100, int(nsteps)); // round up
		float step_weight = fract(nsteps);
	
		// smoke volume iteration using 3D texture, light0 to vposl
		for (int i = 0; i < num_steps; ++i) {
			vec3 p = clamp(pos, 0.0, 1.0); // should be in [0.0, 1.0] range
			float smoke = SMOKE_SCALE*texture3D(smoke_tex, p.zxy).a*step_weight;
			lscale = mix(lscale, 0.0, smoke);
			pos   += delta*step_weight;
			step_weight = 1.0;
		}
	}
	return lscale * add_light_comp_pos_smap(n, epos, ix).rgb;
}

// Note: This may seem like it can go into the vertex shader as well,
//       but we don't have the tex0 value there and can't determine the full init color
void main()
{
	vec3 off   = vec3(-x_scene_size, -y_scene_size, czmin);
	vec3 scale = vec3(2.0*x_scene_size, 2.0*y_scene_size, (czmax - czmin));
	vec3 lit_color = gl_Color.rgb; // base color (with some lighting)
	
	if (indir_lighting) {
		vec3 sp    = clamp((spos  - off)/scale, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec3 indir_light = texture3D(smoke_tex, sp.zxy).rgb; // add indir light color from texture
		lit_color += gl_FrontMaterial.diffuse.rgb * indir_light; // indirect lighting
	}
	if (direct_lighting) { // directional light sources with no attenuation
		vec3 n = normalize(eye_norm);
		if (enable_light0) lit_color += add_light(0, off, scale, n, lpos0, vposl);
		if (enable_light1) lit_color += add_light_comp_pos_smap(n, epos, 1).rgb;
		//if (enable_light0) ADD_LIGHT(0);
		//if (enable_light1) ADD_LIGHT(1);
		if (enable_light2) ADD_LIGHT(2);
		if (enable_light3) ADD_LIGHT(3);
		if (enable_light4) ADD_LIGHT(4);
		if (enable_light5) ADD_LIGHT(5);
		if (enable_light6) ADD_LIGHT(6);
		if (enable_light7) ADD_LIGHT(7);
	}
	if (enable_dlights) lit_color += add_dlights(vpos, normalize(normal), eye); // dynamic lighting
	vec4 texel  = texture2D(tex0, gl_TexCoord[0].st);
	float alpha = gl_Color.a;

	if (do_lt_atten /*&& light_atten > 0.0*/) { // account for light attenuating/reflecting semi-transparent materials
		vec3 view_vec = vpos - eye;
		vec3 far_pt   = vpos + 100.0*view_vec/length(view_vec); // move it far away
		pt_pair res   = clip_line(vpos, far_pt, cube_bb);
		float dist    = length(res.v1 - res.v2);
		float atten   = exp(-light_atten*dist);
		alpha += (1.0 - alpha)*(1.0 - atten);
	}
	if (do_lt_atten) {
		vec3 v_inc = eye - vpos;
		
		if (refract_ix != 1.0 && dot(normal, v_inc) > 0.0) { // entering ray in front surface
			float reflect_w = get_fresnel_reflection(normalize(v_inc), normalize(normal), 1.0, refract_ix);
			alpha = reflect_w + alpha*(1.0 - reflect_w); // don't have a reflection color/texture, so just modify alpha
		} // else exiting ray in back surface - ignore for now since we don't refract the ray
	}
	vec4 color = vec4((texel.rgb * lit_color), (texel.a * alpha));
	if (keep_alpha && color.a <= min_alpha) discard;
	vec3 eye_c  = eye;
	vec3 vpos_c = vpos;

	if (smoke_enabled) {
		pt_pair res = clip_line(vpos, eye, smoke_bb);
		eye_c  = res.v1;
		vpos_c = res.v2;
	}
	if (!smoke_enabled || eye_c == vpos_c) {
		if (color.a <= min_alpha) discard;
		if (!smoke_enabled) color = apply_fog(color); // apply standard fog
		gl_FragColor = color;
		return;
	}
	
	// smoke code
	vec3 dir      = eye_c - vpos_c;
	vec3 pos      = (vpos_c - off)/scale;
	vec3 delta    = normalize(dir)*step_delta/scale;
	float nsteps  = length(dir)/step_delta;
	int num_steps = 1 + min(1000, int(nsteps)); // round up
	float step_weight = fract(nsteps);
	
	// smoke volume iteration using 3D texture, pos to eye
	for (int i = 0; i < num_steps; ++i) {
		vec3 p        = clamp(pos, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec4 tex_val  = texture3D(smoke_tex, p.zxy); // rgba = {color.rgb, smoke}
		float smoke   = SMOKE_SCALE*tex_val.a*step_weight;
		vec3 rgb_comp = (tex_val.rgb * gl_Fog.color.rgb);
		color = ((!keep_alpha && color.a == 0.0) ? vec4(rgb_comp, smoke) :
					mix(color, vec4(rgb_comp, (keep_alpha ? color.a : 1.0)), smoke));
		pos  += delta*step_weight;
		step_weight = 1.0;
	}
	if (color.a <= min_alpha) discard;
	gl_FragColor = color;
}
