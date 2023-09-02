#include <fresnel.part>

const float PI = 3.14159;

//uniform vec3 camera_pos; // from dynamic_lighting.part
uniform float water_depth, water_atten, time;
uniform vec3 uw_atten_max, uw_atten_scale;
uniform sampler2D reflection_tex;

const int MAX_SPLASHES = 32; // must agree with C++ code
uniform vec4 splashes[MAX_SPLASHES]; // {x, y, radius, intensity}

in vec3 vertex_ws; // world space
in vec4 epos, proj_pos;

// underwater attenuation code
void atten_color(inout vec3 color, in float atten) {
	color *= vec3(1.0) - min(uw_atten_max, uw_atten_scale*atten);
}

vec2 get_splash_amplitude() {
	float splash_amp = 0.0;
	float height_sum = 0.0;

	for (int i = 0; i < MAX_SPLASHES; ++i) {
		float height = splashes[i].w;
		if (height == 0.0) continue;
		float radius = splashes[i].z;
		float dist   = distance(vertex_ws.xy, splashes[i].xy);
		if (dist > radius) continue; // is this faster or slower? unclear
		// ideally we should check the room walls and allow them to block the wave, maybe with a ray march through a texture with wall=0 and space=1?
		height     *= 1.0 - min(dist/radius, 1.0);
		splash_amp += height*sin(8.0*PI*dist*water_atten - 4.0*time); // expand outward with time
		height_sum += abs(height);
	} // for i
	return vec2(splash_amp, height_sum);
}

void main() {
	vec3 ws_normal = vec3(0.0, 0.0, 1.0); // +Z
	vec2 uv        = clamp((0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5)), 0.0, 1.0);
	
	// apply ripples when the player moves
	vec2 splash = get_splash_amplitude();
	float foam  = min(2.0*splash.g, 1.0);
	vec2 delta  = clamp(10.0*vec2(dFdx(splash.r), dFdy(splash.r)), -1.0, 1.0); // take derivatves of signed splash amplitude
	uv          = clamp((uv + 0.08*delta), 0.0, 1.0);
	ws_normal   = normalize(ws_normal + vec3(delta, 0.0));
	// is there any way to do refraction here? we don't have any control over drawing of the floor
	// what about adjusting vertex heights in the vertex or tessellation shader?
	
	// apply lighting
	vec3 light_color = vec3(0.0);
	if (enable_dlights) {add_dlights(light_color, vertex_ws, epos, ws_normal, vec3(0.5));} // dynamic lighting
	// what about caustics? we could use a texture, but it's unclear how this should look under many ceiling lights
	
	// attenuate based on depth/distance, assuming a constant water depth and average light coming from above in +Z;
	// this isn't correct at the vertical walls because distance is lower since rays don't reach the floor, but it looks okay with low water depth
	vec3 floor_pos    = vertex_ws - vec3(0.0, 0.0, water_depth);
	float t           = clamp(water_depth/max(0.01*water_depth, abs(camera_pos.z - floor_pos.z)), 0.0, 1.0);
	float camera_path = t*distance(camera_pos, floor_pos); // camera=>floor through water
	float light_path  = water_depth + camera_path; // light=>floor + camera=>floor
	float alpha       = min((0.1 + 1.25*water_atten*camera_path), 1.0);
	vec4 uw_color     = vec4(light_color, alpha);
	atten_color(uw_color.rgb, 1.5*water_atten*light_path);
	vec4 water_color  = mix(uw_color, vec4(light_color, 1.0), foam);

	// apply fresnel reflection and compute final color
	vec3 normal      = normalize(fg_NormalMatrix * ws_normal);
	vec3 epos_n      = normalize(epos.xyz);
	float reflect_w  = get_fresnel_reflection(-epos_n, normal, 1.0, 1.333);
	vec4 reflect_tex = texture(reflection_tex, uv);
	fg_FragColor     = mix(water_color, reflect_tex, reflect_w);
}
