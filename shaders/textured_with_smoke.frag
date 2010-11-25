uniform float x_scene_size, y_scene_size, czmin, czmax; // scene bounds (world space)
uniform float step_delta;
uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform sampler2D tex0;
uniform sampler3D smoke_tex;

varying vec3 eye, vpos; // world space

const float SMOKE_SCALE = 0.25;

#define TEST_CLIP_T(reg, va, vb, vd, vc) {float t = ((va) - (vb))/(vd); if ((vc) > 0.0) {if (t > tmin) tmin = t;} else {if (t < tmax) tmax = t;}}

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec4 color = vec4(texel.rgb * gl_Color.rgb, texel.a * gl_Color.a);
	
	if (!smoke_enabled) {
		gl_FragColor = color;
		return;
	}
	
	// clip to scene bounds
	vec3 v1 = vpos;
	vec3 v2 = eye;
	float tmin = 0.0, tmax = 1.0;
	vec3 dv = v2 - v1;
	TEST_CLIP_T(0x01, smoke_bb[0], v1.x, dv.x,  dv.x); // -x plane
	TEST_CLIP_T(0x02, smoke_bb[1], v1.x, dv.x, -dv.x); // +x plane
	TEST_CLIP_T(0x04, smoke_bb[2], v1.y, dv.y,  dv.y); // -y plane
	TEST_CLIP_T(0x08, smoke_bb[3], v1.y, dv.y, -dv.y); // +y plane
	TEST_CLIP_T(0x10, smoke_bb[4], v1.z, dv.z,  dv.z); // -z plane
	TEST_CLIP_T(0x20, smoke_bb[5], v1.z, dv.z, -dv.z); // +z plane
	
	if (tmin >= tmax) { // clipped away
		gl_FragColor = color;
		return;
	}
	v2  = v1 + dv*tmax;
	v1 += dv*tmin;
	
	vec3 off   = vec3(-x_scene_size, -y_scene_size, czmin);
	vec3 scale = vec3(2.0*x_scene_size, 2.0*y_scene_size, (czmax - czmin));
	vec3 pos   = (v1 - off)/scale;
	vec3 dir   = v2 - v1;
	vec3 delta = normalize(dir)*step_delta/scale;
	float nsteps = length(dir)/step_delta;
	int num_steps = 1 + int(nsteps); // round up
	float step_weight = fract(nsteps);
	
	// smoke volume iteration using 3D texture, pos to eye
	for (int i = 0; i < num_steps; ++i) {
		vec3 p = clamp(pos, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec4 tex_val = texture3D(smoke_tex, p.zxy); // rgba = {color.rgb, smoke}
		float smoke = SMOKE_SCALE*tex_val.a*step_weight;
		color = mix(color, vec4((tex_val.rgb * gl_Fog.color.rgb), 1.0), smoke);
		pos += delta*step_weight;
		step_weight = 1.0;
	}
	gl_FragColor = color;
	// update depth if density is large?
}
