uniform float x_scene_size, y_scene_size, czmin, czmax; // scene bounds (world space)
uniform float step_delta;
uniform sampler2D tex0;
uniform sampler3D smoke_tex;
uniform float min_alpha = 0.0;

varying vec3 eye, vpos, orig_vertex; // world space

const float SMOKE_SCALE = 0.25;

const int MAX_LIGHTS = 256;
// store light_source as: center.xyz, radius, color.rgba
uniform int num_lights = 0;
uniform sampler2D dlights_tex;

const float CTHRESH = 0.02;
const float t_scale = 1.0/MAX_LIGHTS;

vec4 get_tex_val(int i, int which) {
	return texture2D(dlights_tex, vec2(0.5*which, i*t_scale));
}

float get_intensity_at(in vec3 pos, in vec3 off, in vec3 scale, in int i) {
	vec4 pos_r = get_tex_val(i,0);
	float radius = pos_r.w;
	if (radius == 0.0) return get_tex_val(i,1).a; // no falloff
	vec3 center = pos_r.xyz*scale + off;
	if (abs(pos.z - center.z) > radius) return 0.0; // fast test
	float dist = length(pos - center);
	if (dist > radius) return 0.0;
	float rscale = (radius - dist)/radius;
	return rscale*rscale*get_tex_val(i,1).a; // quadratic 1/r^2 attenuation
}

// Note: This may seem like it can go into the vertex shader as well,
//       but we don't have the tex0 value there and can't determine the full init color
void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec4 color0 = gl_Color;
	
	if (num_lights > 0) {
		vec3 llc = vec3(-x_scene_size, -y_scene_size, czmin);
		vec3 urc = vec3( x_scene_size,  y_scene_size, czmax);
		vec3 scale = (urc - llc);
		
		for (int i = 0; i < num_lights; ++i) {
			float cscale = get_intensity_at(orig_vertex, llc, scale, i);
			if (cscale < CTHRESH) continue;
			color0 += vec4(get_tex_val(i,1).rgb, 0.0)*cscale;
		}
	}
	vec4 color = vec4(texel.rgb * color0.rgb, texel.a * color0.a);
	if (keep_alpha && color.a <= min_alpha) discard;
	
	if (eye == vpos) {
		if (color.a <= min_alpha) discard;
		gl_FragColor = color;
		return;
	}
	vec3 off   = vec3(-x_scene_size, -y_scene_size, czmin);
	vec3 scale = vec3(2.0*x_scene_size, 2.0*y_scene_size, (czmax - czmin));
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
		color = ((!keep_alpha && color.a == 0.0) ? vec4(rgb_comp, smoke) : mix(color, vec4(rgb_comp, (keep_alpha ? color.a : 1.0)), smoke));
		pos += delta*step_weight;
		step_weight = 1.0;
	}
	if (color.a <= min_alpha) discard;
	gl_FragColor = color;
	
	//need to convert coordinate space and set in early termination case
	//gl_FragDepth = length(dir); // mean free path distance
}
