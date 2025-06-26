#include <fresnel.part>

const float PI = 3.14159;

uniform float water_depth, water_atten, foam_scale, droplet_scale, time;
uniform float ripple_freq = 10.0;
uniform vec3 uw_atten_max, uw_atten_scale, closest_droplet;
uniform sampler2D reflection_tex, frame_buffer;

struct splash_t {
	vec4 loc_rh; // {x, y, radius, height}
	vec4 bounds; // {x1, y1, x2, y2}
};

const int MAX_SPLASHES = 40; // must agree with C++ code
uniform splash_t splashes[MAX_SPLASHES];

#ifdef USE_SPLASH_LINES
struct splash_line_t {
	vec4 lpts; // line points: {p1x, p1y, p2x, p2y}
	vec4 rh; // {radius, height, pad, pad}
	vec4 bounds; // {x1, y1, x2, y2}
};

const int MAX_SPLASHE_LINES = 20; // must agree with C++ code
uniform splash_line_t splash_lines[MAX_SPLASHE_LINES];

float point_line_seg_dist(vec2 p, vec2 p1, vec2 p2) {
    vec2 center = (p1 + p2) * 0.5;
    float len   = length(p2 - p1);
    vec2 dir    = (p2 - p1) / len;
    vec2 rel_p  = p - center;
    float dist1 = abs(dot(rel_p, vec2(dir.y, -dir.x)));
    float dist2 = abs(dot(rel_p, dir)) - 0.5*len;
    return max(dist1, dist2);
}
#endif // USE_SPLASH_LINES

in vec3 vpos; // world space
in vec4 epos, proj_pos;

// underwater attenuation code
void atten_color(inout vec3 color, in float atten) {
	color *= vec3(1.0) - min(uw_atten_max, uw_atten_scale*sqrt(atten)); // apply scattering and absorption; use sqrt() for gradual transition
}

float calc_water_height(in float dist) {
	return sin(ripple_freq*PI*dist*water_atten - 4.0*time); // expand outward with time
}
vec4 get_splash_amplitude() { // returns {signed_amplitude, abs_mag_amplitude, amp_dx, amp_dy}
	float splash_amp = 0.0;
	float height_sum = 0.0;
	float amp_dx     = 0.0;
	float amp_dy     = 0.0;
	vec2 vpos_dx     = vpos.xy + dFdx(vpos.xy); // vpos shifted one pixel in screen space
	vec2 vpos_dy     = vpos.xy + dFdy(vpos.xy);

	for (int i = 0; i < MAX_SPLASHES; ++i) {
		splash_t s   = splashes[i];
		float height = s.loc_rh.w; // amplitude at center, nominally < 1.0
		if (height == 0.0) continue; // this splash is disabled, skip
		float radius = s.loc_rh.z; // max radius of effect of this splash
		vec2 loc_xy  = s.loc_rh.xy;
		float dist   = distance(vpos.xy, loc_xy);
		if (dist > radius) continue; // is this faster or slower? unclear
		// ideally we should check the room walls and allow them to block the wave, maybe with a ray march through a texture with wall=0 and space=1?
		if (vpos.x < s.bounds.x || vpos.y < s.bounds.y || vpos.x > s.bounds.z || vpos.y > s.bounds.w) continue; // check XY bounds
		dist        = max(dist, 0.01*radius); // avoid a singularity right at the center
		float dr    = dist/radius;
		height     *= 1.0 - min(dr*dr, 1.0); // height decreases sqrt-ish with distance; not quadratic, doesn't preserve water volume
		splash_amp += height*calc_water_height(dist);
		// compute XY derivatives; slower, but using dFdx()/dFdy() on splash_amp produces blocky artifacts
		amp_dx     += height*calc_water_height(distance(vpos_dx, loc_xy));
		amp_dy     += height*calc_water_height(distance(vpos_dy, loc_xy));
		height_sum += abs(height);
	} // for i

#ifdef USE_SPLASH_LINES
	for (int i = 0; i < MAX_SPLASHE_LINES; ++i) {
		splash_line_t sl = splash_lines[i];
		float height = sl.rh.y; // amplitude at center, nominally < 1.0
		if (height == 0.0) continue; // this splash line is disabled, skip
		float radius = sl.rh.x; // max radius of effect of this splash
		float dist   = point_line_seg_dist(vpos.xy, sl.lpts.xy, sl.lpts.zw);
		if (dist > radius) continue; // is this faster or slower? unclear
		if (vpos.x < sl.bounds.x || vpos.y < sl.bounds.y || vpos.x > sl.bounds.z || vpos.y > sl.bounds.w) continue; // check XY bounds
		dist        = max(dist, 0.01*radius); // avoid a singularity right at the center
		float dr    = dist/radius;
		height     *= 1.0 - min(dr*dr, 1.0); // height decreases sqrt-ish with distance; not quadratic, doesn't preserve water volume
		splash_amp += height*calc_water_height(dist);
		// compute XY derivatives
		amp_dx     += height*calc_water_height(point_line_seg_dist(vpos_dx, sl.lpts.xy, sl.lpts.zw));
		amp_dy     += height*calc_water_height(point_line_seg_dist(vpos_dy, sl.lpts.xy, sl.lpts.zw));
		height_sum += abs(height);
	} // for i
#endif // USE_SPLASH_LINES
	return vec4(splash_amp, height_sum, (amp_dx - splash_amp), (amp_dy - splash_amp));
}

void main() {
	vec3 ws_normal = vec3(0.0, 0.0, 1.0); // world space, +Z
	vec2 uv        = clamp((0.5*proj_pos.xy/proj_pos.w + vec2(0.5)), 0.0, 1.0); // base reflection texture coordinate
	
	// apply ripples when the player or other object moves or falls into the water
	vec4 splash = get_splash_amplitude();
	float foam  = foam_scale*min(2.0*splash.y, 0.5); // or mud, depending on the fluid; reduces transparency near the splash center
	vec2 delta  = clamp(10.0*vec2(splash.z, splash.w), -1.0, 1.0);
	vec2 uv_adj = clamp((uv + 0.1*delta), 0.0, 1.0); // offset by ripple angle, but keep in the valid texture range
	ws_normal   = normalize(ws_normal + vec3(delta, 0.0)); // average with the incoming (+Z) normal
	
	// apply lighting
	vec3 light_color = vec3(0.0);
	//add_indir_lighting(light_color, 1.0); // indir lighting is disabled as it has little effect and is difficult to setup
	if (enable_dlights) {add_dlights(light_color, vpos, epos, ws_normal, vec3(0.5));} // add dynamic lighting; if disabled, water will be black
	
	// attenuate based on depth/distance, assuming a constant water depth and average light coming from above in +Z;
	float eye_to_floor = get_linear_depth_zval(uv); // floor, wall, object, etc.
	float eye_to_water = log_to_linear_depth(gl_FragCoord.z);
	float camera_path  = eye_to_floor - eye_to_water; // camera=>fragment through water
	float light_path   = water_depth + camera_path; // light=>floor + camera=>floor; may be approximate if water_depth isn't accurate
	float alpha        = min((0.3 + 0.6*sqrt(water_atten*camera_path)), 1.0); // short path is transparent, long path is opaque

	if (closest_droplet != vec3(0.0)) { // darken/more opaque around droplets
		alpha += max(0.0, 0.33*(droplet_scale - distance(vpos, closest_droplet))/droplet_scale);
	}
	vec4 uw_color      = vec4(light_color, alpha);
	atten_color(uw_color.rgb, 1.5*water_atten*light_path); // apply scattering and absorption
	vec4 water_color   = mix(uw_color, vec4(light_color, 1.0), foam);

	// apply fresnel reflection and compute final color
	vec3 normal      = normalize(fg_NormalMatrix * ws_normal); // eye space
	vec3 epos_n      = normalize(epos.xyz); // eye direction
	float reflect_w  = get_fresnel_reflection(-epos_n, normal, 1.0, 1.333);
	vec4 reflect_tex = texture(reflection_tex, uv_adj);
	vec4 surf_color  = mix(water_color, reflect_tex, reflect_w);

	// handle refractions by mixing the frame buffer (with the floor and underwater objects) with the water using the alpha value
	float darken_fact  = 0.8; // darken underwater objects to simulate wetness; is there a way to incrementally darken walls above the water line as well?
	vec3 refract_color = darken_fact*texture(frame_buffer, uv_adj).rgb;
	// adjust intensity by the ratio of signed amplitude to height sum to get a signed value within an envelope that falls off at both ends:
	// at the low radius (to prevent sharp edges near the singularity at the splash center), and at the outer radius for a smooth transition
	float inten_adj    = 0.2*splash.x/max(splash.y, 0.05);
	refract_color     *= 1.0 + clamp(inten_adj, 0.0, 1.0); // add fake caustics based on splash amplitude
	fg_FragColor       = vec4(mix(refract_color, surf_color.rgb, surf_color.a), 1.0);
}
