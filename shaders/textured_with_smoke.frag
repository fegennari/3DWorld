uniform float x_scene_size, y_scene_size, czmin, czmax; // scene bounds (world space)
uniform float step_delta;
uniform sampler2D tex0;
uniform sampler3D smoke_tex;

varying vec3 eye, vpos; // world space

const float SMOKE_MAX_CELL = 0.125;

#define TEST_CLIP_T(reg, va, vb, vd, vc) {float t = ((va) - (vb))/(vd); if ((vc) > 0.0) {if (t > tmin) tmin = t;} else {if (t < tmax) tmax = t;}}

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec4 color = vec4(texel.rgb * gl_Color.rgb, texel.a * gl_Color.a);
	
	// clip to scene bounds
	vec3 v1 = vpos;
	vec3 v2 = eye;
	float tmin = 0.0, tmax = 1.0;
	vec3 dv = v2 - v1;
	TEST_CLIP_T(0x01, -x_scene_size, v1.x, dv.x,  dv.x); // -x plane
	TEST_CLIP_T(0x02,  x_scene_size, v1.x, dv.x, -dv.x); // +x plane
	TEST_CLIP_T(0x04, -y_scene_size, v1.y, dv.y,  dv.y); // -y plane
	TEST_CLIP_T(0x08,  y_scene_size, v1.y, dv.y, -dv.y); // +y plane
	TEST_CLIP_T(0x10,  czmin,        v1.z, dv.z,  dv.z); // -z plane
	TEST_CLIP_T(0x20,  czmax,        v1.z, dv.z, -dv.z); // +z plane
	
	if (tmin >= tmax) { // clipped away
		gl_FragColor = color;
		return;
	}
	// FIXME: colors above house chimney
	v2  = v1 + dv*tmax;
	v1 += dv*tmin;
	
	vec3 off   = vec3(-x_scene_size, -y_scene_size, czmin);
	vec3 scale = vec3(2.0*x_scene_size, 2.0*y_scene_size, (czmax - czmin));
	vec3 pos   = (v1 - off)/scale;
	vec3 dir   = v2 - v1;
	//vec3 delta = normalize(dir)*step_delta/scale;
	int num_steps = 1 + min(1000, int(length(dir)/step_delta)); // round up
	vec3 delta = ((v2 - off)/scale - pos)/float(num_steps);
	
	// smoke volume iteration using 3D texture, pos to eye
	for (int i = 0; i < num_steps; ++i) {
		// FIXME: smoother steps
		vec3 p = clamp(pos, 0.0, 1.0); // should be in [0.0, 1.0] range
		vec4 tex_val = texture3D(smoke_tex, p.zxy); // rgba = {color.rgb, smoke}
		float smoke = 2.0*SMOKE_MAX_CELL*tex_val.a;
		color = mix(color, vec4((tex_val.rgb * gl_Fog.color.rgb), 1.0), smoke);
		pos += delta;
	}
	gl_FragColor = color;
	// update depth if density is large?
}
