#include <fresnel.part>

//uniform vec3 camera_pos; // from dynamic_lighting.part
uniform float water_depth, water_atten;
uniform vec3 uw_atten_max, uw_atten_scale;
uniform sampler2D reflection_tex;

const int MAX_SPLASHES = 16; // must agree with C++ code
uniform vec4 splashes[MAX_SPLASHES]; // {x, y, radius, intensity}

in vec3 vertex_ws; // world space
in vec4 epos, proj_pos;

// underwater attenuation code
void atten_color(inout vec3 color, in float atten) {
	color *= vec3(1.0) - min(uw_atten_max, uw_atten_scale*atten);
}

void main() {
	vec3 ws_normal = vec3(0.0, 0.0, 1.0); // +Z
	vec2 uv        = clamp((0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5)), 0.0, 1.0);
	
	// apply ripples when the player moves
	// TODO
	
	// apply lighting
	vec3 water_color = vec3(0.0);
	if (enable_dlights) {add_dlights(water_color.rgb, vertex_ws, epos, ws_normal, vec3(0.5));} // dynamic lighting
	
	// attenuate based on depth/distance, assuming a constant water depth and average light coming from above in +Z
	vec3 floor_pos    = vertex_ws - vec3(0.0, 0.0, water_depth);
	float t           = clamp(water_depth/max(0.01*water_depth, abs(camera_pos.z - floor_pos.z)), 0.0, 1.0);
	float camera_path = t*distance(camera_pos, floor_pos); // camera=>floor through water
	float light_path  = water_depth + camera_path; // light=>floor + camera=>floor
	float alpha       = min((0.1 + 1.25*water_atten*camera_path), 1.0);
	atten_color(water_color, 1.5*water_atten*light_path);

	// apply fresnel reflection
	vec3 normal      = normalize(fg_NormalMatrix * ws_normal);
	vec3 epos_n      = normalize(epos.xyz);
	float reflect_w  = get_fresnel_reflection(-epos_n, normal, 1.0, 1.333);
	vec4 reflect_tex = texture(reflection_tex, uv);
	fg_FragColor     = mix(vec4(water_color, alpha), reflect_tex, reflect_w);

	// debug visualization
	float splash_amplitude = 0.0;

	for (int i = 0; i < MAX_SPLASHES; ++i) {
		float height = splashes[i].w;
		if (height == 0.0) continue;
		float radius = splashes[i].z;
		float dist   = distance(vertex_ws.xy, splashes[i].xy);
		if (dist > radius) continue; // TODO: is this faster or slower?
		height           *= 1.0 - min(dist/radius, 1.0);
		splash_amplitude += height*sin(8.0*dist*water_atten);
	} // for i
	if (splash_amplitude != 0.0) {
		splash_amplitude = clamp(splash_amplitude, -1.0, 1.0);
		fg_FragColor = vec4(max(splash_amplitude, 0.0), 0.0, max(-splash_amplitude, 0.0), 1.0); // red-blue
	}
}
