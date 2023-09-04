#include <fresnel.part>

const float PI = 3.14159;

//uniform vec3 camera_pos; // from dynamic_lighting.part
uniform float water_depth, water_atten, time;
uniform vec3 uw_atten_max, uw_atten_scale;
uniform sampler2D reflection_tex, frame_buffer;

struct splash_t {
	vec4 loc_rh; // {x, y, radius, height}
	vec4 bounds; // {x1, y1, x2, y2}
};

const int MAX_SPLASHES = 40; // must agree with C++ code
uniform splash_t splashes[MAX_SPLASHES];

in vec3 vertex_ws; // world space
in vec4 epos, proj_pos;

// underwater attenuation code
void atten_color(inout vec3 color, in float atten) {
	color *= vec3(1.0) - min(uw_atten_max, uw_atten_scale*atten); // apply scattering and absorption; faster and more tweak-able than using exp()
}

float calc_water_height(in float dist) {
	return sin(10.0*PI*dist*water_atten - 4.0*time); // expand outward with time
}
vec4 get_splash_amplitude(in vec3 vpos) { // returns {signed_amplitude, abs_mag_amplitude, amp_dx, amp_dy}
	float splash_amp = 0.0;
	float height_sum = 0.0;
	float amp_dx     = 0.0;
	float amp_dy     = 0.0;
	vec2 vpos_dx     = vpos.xy + dFdx(vpos.xy); // vpos shifted one pixel in screen space
	vec2 vpos_xy     = vpos.xy + dFdy(vpos.xy);

	for (int i = 0; i < MAX_SPLASHES; ++i) {
		float height = splashes[i].loc_rh.w; // amplitude at center, nominally < 1.0
		if (height == 0.0) continue;  // this splash is disabled, skip
		float radius = splashes[i].loc_rh.z; // max radius of effect of this splash
		vec2 loc_xy  = splashes[i].loc_rh.xy;
		float dist   = distance(vpos.xy, loc_xy);
		if (dist > radius) continue; // is this faster or slower? unclear
		// ideally we should check the room walls and allow them to block the wave, maybe with a ray march through a texture with wall=0 and space=1?
		vec4 bounds = splashes[i].bounds;
		if (vpos.x < bounds.x || vpos.y < bounds.y || vpos.x > bounds.z || vpos.y > bounds.w) continue; // check XY bounds
		dist        = max(dist, 0.01*radius); // avoid a singularity right at the center
		float dr    = dist/radius;
		height     *= 1.0 - min(dr*dr, 1.0); // height decreases sqrt-ish with distance; not quadratic, doesn't preserve water volume
		splash_amp += height*calc_water_height(dist);
		// compute XY derivatives; slower, but using dFdx()/dFdy() on splash_amp produces blocky artifacts
		amp_dx     += height*calc_water_height(distance(vpos_dx, loc_xy));
		amp_dy     += height*calc_water_height(distance(vpos_xy, loc_xy));
		height_sum += abs(height);
	} // for i
	return vec4(splash_amp, height_sum, (amp_dx - splash_amp), (amp_dy - splash_amp));
}

void main() {
	vec3 ws_normal = vec3(0.0, 0.0, 1.0); // world space, +Z
	vec2 uv        = clamp((0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5)), 0.0, 1.0); // base reflection texture coordinate
	
	// apply ripples when the player or other object moves or falls into the water
	vec4 splash = get_splash_amplitude(vertex_ws);
	float foam  = min(2.0*splash.y, 0.5); // or mud, depending on the fluid; reduces transparency near the splash center
	vec2 delta  = clamp(10.0*vec2(splash.z, splash.w), -1.0, 1.0);
	uv          = clamp((uv + 0.1*delta), 0.0, 1.0); // offset by ripple angle, but keep in the valid texture range
	ws_normal   = normalize(ws_normal + vec3(delta, 0.0)); // average with the incoming (+Z) normal
	
	// apply lighting
	vec3 light_color = vec3(0.0);
	if (enable_dlights) {add_dlights(light_color, vertex_ws, epos, ws_normal, vec3(0.5));} // add dynamic lighting; if disabled, water will be black
	
	// attenuate based on depth/distance, assuming a constant water depth and average light coming from above in +Z;
	// this isn't correct at the vertical walls because distance is lower since rays don't reach the floor, but it looks okay with low water depth
	vec3 floor_pos    = vertex_ws - vec3(0.0, 0.0, water_depth); // bottom surface of water which will be visible through/under the water
	float t           = clamp(water_depth/max(0.01*water_depth, abs(camera_pos.z - floor_pos.z)), 0.0, 1.0);
	float camera_path = t*distance(camera_pos, floor_pos); // camera=>floor through water
	float light_path  = water_depth + camera_path; // light=>floor + camera=>floor
	float alpha       = min((0.1 + 1.25*water_atten*camera_path), 1.0); // short path is transparent, long path is opaque
	vec4 uw_color     = vec4(light_color, alpha);
	atten_color(uw_color.rgb, 1.5*water_atten*light_path); // apply scattering and absorption
	vec4 water_color  = mix(uw_color, vec4(light_color, 1.0), foam);

	// apply fresnel reflection and compute final color
	vec3 normal      = normalize(fg_NormalMatrix * ws_normal); // eye space
	vec3 epos_n      = normalize(epos.xyz); // eye direction
	float reflect_w  = get_fresnel_reflection(-epos_n, normal, 1.0, 1.333);
	vec4 reflect_tex = texture(reflection_tex, uv);
	vec4 surf_color  = mix(water_color, reflect_tex, reflect_w);

	// handle refractions by mixing the frame buffer (with the floor and underwater objects) with the water using the alpha value
	vec3 refract_color = texture(frame_buffer, uv).rgb;
	// adjust intensity by the ratio of signed amplitude to height sum to get a signed value within an envelope that falls off at both ends:
	// at the low radius (to prevent sharp edges near the singularity at the splash center), and at the outer radius for a smooth transition
	float inten_adj    = 0.2*splash.x/max(splash.y, 0.05);
	refract_color     *= 1.0 + clamp(inten_adj, 0.0, 1.0); // add fake caustics based on splash amplitude
	fg_FragColor       = vec4(mix(refract_color, surf_color.rgb, surf_color.a), 1.0);
}
