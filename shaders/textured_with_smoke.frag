uniform float x_scene_size, y_scene_size, czmin, czmax; // scene bounds (world space)
uniform float step_delta;
uniform sampler2D tex0;
uniform sampler3D smoke_tex;

varying vec3 eye, vpos; // world space

const float SMOKE_MAX_CELL = 0.125;

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec4 color = vec4(texel.rgb * gl_Color.rgb, texel.a * gl_Color.a);
	
	float density = 0.0;
	float scolor  = 0.0;
	// FIXME: clip to scene bounds
	vec3 pos   = vpos;
	vec3 dir   = eye - pos;
	vec3 delta = normalize(dir)*step_delta;
	int num_steps = 1 + int(length(dir)/step_delta); // round up
	vec3 off   = vec3(-x_scene_size, -y_scene_size, czmin);
	vec3 scale = vec3(0.5/x_scene_size, 0.5/y_scene_size, 1.0/(czmax - czmin));
	
	// smoke volume iteration using 3D texture, pos to eye
	for (int i = 0; i < num_steps; ++i) {
		vec3 p = clamp(((pos - off) * scale), 0.0, 1.0); // should be in [0.0, 1.0] range
		vec4 val = texture3D(smoke_tex, p.zxy); // rgb = {smoke, val, luminance}
		float smoke = min(SMOKE_MAX_CELL, val.r);
			
		if (smoke > 0.0) {
			float val = 0.5*(val.g + val.b); // add in luminance from lights
			scolor    = ((density == 0.0) ? val : (smoke*val + (1.0 - smoke)*scolor));
			density  += smoke;
		}
		pos += delta;
	}
	float brightness = density*scolor;
	color.a = (1.0 - density)*color.a + density; // deal with transparent objects
	if (scolor < 1.0) color *= (1.0 - density)/(1.0 - brightness); // BLACK: attenuation
	float fog_coord = 15.0*brightness; // WHITE(GRAY): fog coord
	float fog = clamp((gl_Fog.end - fog_coord) * gl_Fog.scale, 0.0, 1.0);
	gl_FragColor = mix(gl_Fog.color, color, fog);
}
