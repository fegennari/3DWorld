// http://glslsandbox.com/e#20376.4

// Planet rendering
// Written by Frank Gennari
// 10/6/14

// Note: This shader was ported from 3DWorld,
// which does the planet rendering entirely in the fragment shader;
// however, there are some problems with this approach:
// 1. Instead of passing the noise values in as 2D and 3D textures,
// they're computed in the shader, which makes it take forever to compile and run slowly
// 2. The large block of constants at the top is normally computed on the CPU
// and passed in as uniforms, so the values may work better and be more consistent
// with each other than the values used here
// 3. The original shader used both world and eye space for the transforms,
// and this shader probably has some transform/lighting bugs due to the porting
// 4. The original shader was in an entirely different coordinate scale from this one,
// and I haven't had time to tune all of the various constants
// 5. The original shader had variable bounds for some for loops,
// but they are hard-coded here, which reduces the quality/performance tradeoff
// 6. The original shader was multiple files + code generation,
// and bugs may have been introduced in the process of merging all of the code together

#ifdef GL_ES
precision mediump float;
#endif

// standard uniforms
uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;

// global constants (to be changed by the user)
const vec3 planet_pos     = vec3(0,0,0);
const float viewer_z      = 2.0;
const float planet_radius = 0.25;
const float moon_radius   = 0.06;
const float atmos_radius  = 0.26;
const float sun_radius    = 0.4; // far away
const float orbit         = 0.5;
const float terrain_scale = 1.0;
const float noise_scale   = 1.0;
const float atmosphere    = 1.0;
const float population    = 0.0;
const float water_val     = 0.2;
const float lava_val      = 0.0;
const float snow_thresh   = 0.9;
const float temperature   = 20.0;
const float spec_shine    = 50.0;
const vec3 cloud_freq     = vec3(2.0);
const vec3 atmos_density  = vec3(0.0, 2.0, 0.0); // {constant, linear, quadratic}
const vec3 atmos_inner_col = vec3(1.0, 1.0, 1.0); // white
const vec3 atmos_outer_col = vec3(0.0, 0.0, 1.0); // blue
const vec4 water_color    = vec4(0.2, 0.2, 0.8, 1.0); // blue-ish
const vec4 color_a        = vec4(0.05, 0.35, 0.05, 1.0); // dark green vegetation
const vec4 color_b        = vec4(0.60, 0.45, 0.25, 1.0); // brown dirt
const vec4 sun_color      = vec4(1.0, 1.0, 1.0, 1.0);
const vec4 specular_color = vec4(1.0, 1.0, 1.0, 1.0);
#define NUM_OCTAVES 8

// v1 is the line dir, p1 is the line start point, sphere center is at (0,0,0), sphere radius is 1.0
vec3 line_sphere_int(in vec3 v1, in vec3 p1) {
	float t = -dot(v1, p1); // v1 should be normalized
	vec3 v2 = p1 + v1*t;
	float v2len = length(v2); // sqrt of length of perpendicular to sphere center
	return v2 - v1*sqrt(1.0 - v2len*v2len); // distance along line to outside surface of sphere
}

vec3 ray_hit_dir(in vec2 screen_pos) {
	vec3 p1 = vec3(screen_pos, viewer_z);
	return normalize(line_sphere_int(vec3(0,0,1), p1));
}

vec3 draw_simple_planet(in vec3 color, in vec3 normal, in vec3 sun_dir, in float diffuse_scale);
void apply_atmosphere(inout vec3 cur_color, in vec3 sun_pos, in vec3 moon_pos, in vec3 light_dir,
		      in vec3 normal, in vec3 position, in vec3 camera_pos);
float simplex(in vec2 v);
float calc_sphere_shadow_atten(in vec3 pos, in vec3 lpos, in float lradius, in vec3 spos, in float sradius);
void adjust_normal_for_craters(inout vec3 norm, in vec3 vertex);
vec3 draw_planet(in vec3 vertex, in vec3 normal, in vec3 light_dir, in vec3 sun_pos, in vec3 ss_pos, in float ss_radius);

void main(void) {
	float rmin = min(resolution.x, resolution.y);
	vec2 position  = ((gl_FragCoord.xy - 0.5*resolution.xy) / rmin);

	vec3 sun_dir = normalize(vec3(2.0*mouse.x-1.0, 2.0*mouse.y-1.0, -0.5));
	vec3 sun_pos = 10.0*sun_dir;
	float angle  = 0.25*time;
	vec3 moon_pos = planet_pos + orbit*vec3(sin(angle), 0.0, cos(angle)); // orbits around y-axis
	vec2 ppos = position - planet_pos.xy;
	vec2 mpos = position - moon_pos.xy;
	bool in_planet = (length(ppos) < planet_radius);
	bool in_moon   = (length(mpos) < moon_radius);
	vec3 p_normal  = ray_hit_dir(ppos/planet_radius);
	vec3 real_ppos = planet_pos + p_normal*planet_radius;
	vec3 color;

	if (in_moon && (!in_planet || moon_pos.z < 0.0)) { // in moon (and in front of planet)
		vec3 normal = ray_hit_dir(mpos/moon_radius);
		vec3 pos = moon_pos + normal*moon_radius;
		float ds = calc_sphere_shadow_atten(pos, sun_pos, sun_radius, planet_pos, planet_radius);
		adjust_normal_for_craters(normal, normal); // slow
		color = draw_simple_planet(vec3(0.5,0.5,0.5), normal, sun_dir, ds);
	}
	else if (in_planet) {
		//float ds = calc_sphere_shadow_atten(real_ppos, sun_pos, sun_radius, moon_pos, moon_radius);
		//color = draw_simple_planet(vec3(0.5,0.5,1.0), p_normal, sun_dir, ds);
		color = draw_planet(real_ppos, p_normal, sun_dir, sun_pos, moon_pos, moon_radius);
	}
	else { // background stars
		float intensity = min(1.0, 4.0*(simplex(100.0*position) - 0.8));
		color = vec3(intensity);
	}
	if (length(ppos) < atmos_radius && (!in_moon || moon_pos.z > 0.0)) {
		vec3 camera_pos = vec3(position, viewer_z);
		apply_atmosphere(color, sun_pos, moon_pos, sun_dir, p_normal, real_ppos, camera_pos);
	}
	gl_FragColor = sun_color * vec4(color, 1.0);
}

vec3 draw_simple_planet(in vec3 color, in vec3 normal, in vec3 sun_dir, in float diffuse_scale) {
	float diffuse = 0.95*diffuse_scale*max(0.0, dot(normal, sun_dir));
	float ambient = 0.05;
	return color*(diffuse + ambient);
}

// ****************** SPHERE SHADOW ************************

float pt_line_dist(in vec3 P, in vec3 L1, in vec3 L2) {
	return length(cross((L2 - L1), (L1 - P)))/distance(L2, L1);
}

// analytical soft sphere shadow with spherical light source
float calc_sphere_shadow_atten(in vec3 pos, in vec3 lpos, in float lradius, in vec3 spos, in float sradius) {
	float atten = 1.0;
	float ldist = length(lpos - pos);

	if (ldist > length(lpos - spos)) { // behind the shadowing object
		const float PI = 3.14159;
		float d = pt_line_dist(spos, lpos, pos);
		float r = sradius;
		float R = lradius*length(spos - pos)/ldist;
		float tot_area = PI*R*R;

		if (d < abs(R - r)) { // fully overlapped
			atten *= 1.0 - PI*min(r,R)*min(r,R)/tot_area;
		}
		else if (d < (r + R)) { // partially overlapped
			float shadowed_area = r*r*acos((d*d+r*r-R*R)/(2.0*d*r)) + R*R*acos((d*d+R*R-r*r)/(2.0*d*R)) - 0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
			atten *= 1.0 - clamp(shadowed_area/tot_area, 0.0, 1.0);
		}
	}
	return atten;
}

// ****************** ATMOSPHERE ************************

float get_density_at(vec3 pos) {
	float dist = distance(pos, planet_pos);
	float v = 1.0 - clamp((dist - planet_radius)/(atmos_radius - planet_radius), 0.0, 1.0);
	return dot(atmos_density, vec3(1.0, v, v*v));
}

void apply_atmosphere(inout vec3 cur_color, in vec3 sun_pos, in vec3 moon_pos, in vec3 light_dir, in vec3 normal, in vec3 position, in vec3 camera_pos) {
	float ascale = min(4.0*(dot(normal, light_dir) + 0.25), 1.0);
	if (ascale <= 0.0) return;

	// alpha is calculated from distance between sphere intersection points
	float wpdist   = distance(position, planet_pos);
	vec3 ldir      = normalize(position - camera_pos);
	float dp       = dot(ldir, (position - planet_pos));
	float adist_sq = dp*dp - wpdist*wpdist + atmos_radius*atmos_radius;
	if (adist_sq <= 0.0) discard; // no sphere intersection
	float dist     = sqrt(adist_sq);
	float pdist_sq = dp*dp - wpdist*wpdist + planet_radius*planet_radius;
	if (pdist_sq > 0.0) {dist -= sqrt(pdist_sq);} // ray intersects planet, adjust distance
	vec3 pos       = position - ldir*(dp + 0.5*dist); // midpoint of ray in atmosphere
	float density  = get_density_at(pos)*dist/atmos_radius;
	float alpha    = ascale*clamp(3.5*density, 0.0, 1.0);
	float lt_atten = 1.0;

	if (sun_radius > 0.0 && moon_radius > 0.0) {
		lt_atten *= calc_sphere_shadow_atten(position, sun_pos, sun_radius, moon_pos, moon_radius);
	}
	vec3 color = vec3(0.05); // ambient only
	color += lt_atten*max(0.0, dot(normal, light_dir)); // diffuse
	vec3 scatter_color = mix(atmos_outer_col, atmos_inner_col, min(1.6*density, 1.0));
	cur_color = mix(cur_color, color*scatter_color, alpha);
}

// ****************** SIMPLEX NOISE ************************

vec3 mod289(in vec3 x) {
	return x - floor(x * 1.0 / 289.0) * 289.0;
}
vec3 permute(in vec3 x) {
	return mod289(((x * 34.0) + 1.0) * x);
}
vec4 mod289(in vec4 x) {
	return x - floor(x * 1.0 / 289.0) * 289.0;
}
vec4 permute(in vec4 x) {
	return mod289(((x * 34.0) + 1.0) * x);
}

float hash( in vec2 p ) {
	float h = dot(p,vec2(127.1,311.7));	
    return fract(sin(h)*43758.5453123);
}
float noise( in vec2 p ) {
    vec2 i = floor( p );
    vec2 f = fract( p );	
	vec2 u = f*f*(3.0-2.0*f);
    return -1.0+2.0*mix( mix( hash( i + vec2(0.0,0.0) ), hash( i + vec2(1.0,0.0) ), u.x),
                         mix( hash( i + vec2(0.0,1.0) ), hash( i + vec2(1.0,1.0) ), u.x), u.y);
}

float simplex(in vec2 v)
{
#if 0
	return noise(v);
#else
	vec4 C = vec4(
		 0.211324865405187,  // (3.0 -  sqrt(3.0)) / 6.0
		 0.366025403784439,  //  0.5 * (sqrt(3.0)  - 1.0)
		-0.577350269189626,	 // -1.0 + 2.0 * C.x
		 0.024390243902439); //  1.0 / 41.0

	// First corner
	vec2 i  = floor(v + dot(v, C.yy));
	vec2 x0 = v - i   + dot(i, C.xx);

	// Other corners
	vec2 i1 = (x0.x > x0.y) ? vec2(1, 0) : vec2(0, 1);
	vec4 x12 = x0.xyxy + C.xxzz;
	x12 = vec4(x12.xy - i1, x12.z, x12.w);

	// Permutations
	i = mod(i, vec2(289)); // Avoid truncation effects in permutation
	vec3 p = permute(permute(i.y + vec3(0, i1.y, 1)) + i.x + vec3(0, i1.x, 1));
	vec3 m = max(vec3(0.5) - vec3(dot(x0, x0), dot(x12.xy, x12.xy), dot(x12.zw, x12.zw)), vec3(0));
	m = m * m;
	m = m * m;

	// Gradients: 41 points uniformly over a line, mapped onto a diamond.
	// The ring size 17*17 = 289 is close to a multiple of 41 (41*7 = 287)
	vec3 x = 2.0 * fract(p * C.w) - 1.0;
	vec3 h = abs(x) - 0.5;
	vec3 ox = floor(x + 0.5);
	vec3 a0 = x - ox;

	// Normalise gradients implicitly by scaling m
	// Inlined for speed: m *= taylorInvSqrt(a0*a0 + h*h);
	m *= 1.79284291400159 - 0.85373472095314 * (a0 * a0 + h * h);

	// Compute final noise value at P
	vec3 g;
	g.x  = a0.x  * x0.x   + h.x  * x0.y;
	g.yz = a0.yz * x12.xz + h.yz * x12.yw;
	return 130.0 * dot(m, g);
#endif
}

float simplex(in vec3 v)
{
	vec2 C = vec2(1.0 / 6.0, 1.0 / 3.0);
	vec4 D = vec4(0.0, 0.5, 1.0, 2.0);

	// First corner
	vec3 i = floor(v + dot(v, vec3(C.y)));
	vec3 x0 = v - i + dot(i, vec3(C.x));

	// Other corners
	vec3 g = step(vec3(x0.y, x0.z, x0.x), x0);
	vec3 l = vec3(1.0 - g);
	vec3 i1 = min(g, vec3(l.z, l.x, l.y));
	vec3 i2 = max(g, vec3(l.z, l.x, l.y));

	//   x0 = x0 - 0.0 + 0.0 * C.xxx;
	//   x1 = x0 - i1  + 1.0 * C.xxx;
	//   x2 = x0 - i2  + 2.0 * C.xxx;
	//   x3 = x0 - 1.0 + 3.0 * C.xxx;
	vec3 x1 = x0 - i1 + C.x;
	vec3 x2 = x0 - i2 + C.y; // 2.0*C.x = 1/3 = C.y
	vec3 x3 = x0 - D.y;      // -1.0+3.0*C.x = -0.5 = -D.y

	// Permutations
	i = mod289(i); 
	vec4 p = vec4(permute(permute(permute(i.z + vec4(0.0, i1.z, i2.z, 1.0)) + i.y + vec4(0.0, i1.y, i2.y, 1.0)) + i.x + vec4(0.0, i1.x, i2.x, 1.0)));

	// Gradients: 7x7 points over a square, mapped onto an octahedron.
	// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
	float n_ = 0.142857142857; // 1.0/7.0
	vec3 ns = n_ * vec3(D.w, D.y, D.z) - vec3(D.x, D.z, D.x);

	vec4 j = vec4(p - 49.0 * floor(p * ns.z * ns.z));  //  mod(p,7*7)
	vec4 x_ = vec4(floor(j * ns.z));
	vec4 y_ = vec4(floor(j - 7.0 * x_));    // mod(j,N)
	vec4 x = vec4(x_ * ns.x + ns.y);
	vec4 y = vec4(y_ * ns.x + ns.y);
	vec4 h = vec4(1.0 - abs(x) - abs(y));
	vec4 b0 = vec4(x.x, x.y, y.x, y.y);
	vec4 b1 = vec4(x.z, x.w, y.z, y.w);

	// vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
	// vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
	vec4 s0 = vec4(floor(b0) * 2.0 + 1.0);
	vec4 s1 = vec4(floor(b1) * 2.0 + 1.0);
	vec4 sh = -step(h, vec4(0.0));

	vec4 a0 = vec4(b0.x, b0.z, b0.y, b0.w) + vec4(s0.x, s0.z, s0.y, s0.w) * vec4(sh.x, sh.x, sh.y, sh.y);
	vec4 a1 = vec4(b1.x, b1.z, b1.y, b1.w) + vec4(s1.x, s1.z, s1.y, s1.w) * vec4(sh.z, sh.z, sh.w, sh.w);
	vec3 p0 = vec3(a0.x, a0.y, h.x);
	vec3 p1 = vec3(a0.z, a0.w, h.y);
	vec3 p2 = vec3(a1.x, a1.y, h.z);
	vec3 p3 = vec3(a1.z, a1.w, h.w);

	// Normalise gradients
	//vec4 norm = taylorInvSqrt(vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3)));
	vec4 norm = 1.79284291400159 - 0.85373472095314 * vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3));
	p0 *= norm.x;
	p1 *= norm.y;
	p2 *= norm.z;
	p3 *= norm.w;

	// Mix final noise value
	vec4 m = max(0.6 - vec4(dot(x0, x0), dot(x1, x1), dot(x2, x2), dot(x3, x3)), vec4(0));
	m = m * m;
	return 42.0 * dot(m * m, vec4(dot(p0, x0), dot(p1, x1), dot(p2, x2), dot(p3, x3)));
}

// ****************** RAND GEN ************************

float rand_01  (float val) {return fract(sin(12.9898 * val) * 43758.5453);}
float rand_pm1 (float val) {return 2.0*(rand_01(val) - 0.5);}
vec3  rand_vec3(float val) {return vec3(rand_pm1(val), rand_pm1(val+1.0), rand_pm1(val+2.0));}

// ****************** CRATERS ************************

// add craters by modifying the normal - assumes the object center is at (0,0,0)
void adjust_normal_for_craters(inout vec3 norm, in vec3 vertex) {
	
	float v0 = 1.0; // using a variable here is slow
	vec3 dir = normalize(vertex); // world space normal

	for (int i = 0; i < 50; ++i) { // Note: inefficient
		vec3 center = rand_vec3(v0);
		vec3 dir2   = dir - normalize(center);
		float dist  = length(dir2);
		float rad1  = 0.15*(0.25 + 0.75*rand_01(v0+3.0));
		float rad2  = 1.5*rad1;
		v0         += 4.0;
		
		if (dist < rad2) { // at crater (parabola)
			vec3 cnorm = normalize(dir2/dist);
			float cwt;

			if (dist < rad1) { // inside crater
				cwt  = 0.75*dist/rad1; // higher power?
				cnorm = -cnorm;
			}
			else { // on rim of crater
				cwt  = 0.5*sqrt(1.0 - (dist - rad1)/(rad2 - rad1));
			}
			norm = normalize(mix(norm, cnorm, smoothstep(0.0, 1.0, cwt)));
		}
	}
}

// ****************** PERLIN CLOUDS 3D ************************

#define RIDGED_NOISE

float gen_cloud_alpha_time(in vec3 pos, in vec3 ftime)
{
	float alpha = 0.0;
	float freq  = 1.0;

	for (int i = 0; i < NUM_OCTAVES; ++i) {
		float v = simplex(noise_scale*(freq*pos + ftime));
#ifdef RIDGED_NOISE
		v = 2.0*v - 1.0; // map [0,1] range to [-1,1]
		v = 1.0 - abs(v); // ridged noise
		v = v*v;
#endif
		alpha += v/freq;
		freq  *= 2.0;
	}
	return 2.0*(0.5*alpha-0.4);
}

float fract_smooth(in float t) {
	return 2.0*abs(fract(0.5*t) - 0.5);
}

vec3 get_ftime() {
	float time2 = 0.01*time;
	return vec3(fract_smooth(time2), fract_smooth(0.95*time2), fract_smooth(0.85*time2));
}

float gen_cloud_alpha_non_norm(in vec3 pos) {
	return gen_cloud_alpha_time(pos, get_ftime());
}
float gen_cloud_alpha(in vec3 pos) {
	return clamp(gen_cloud_alpha_non_norm(pos), 0.0, 1.0);
}
float gen_cloud_alpha_static_non_norm(in vec3 pos) {
	return gen_cloud_alpha_time(pos, vec3(0.0));
}

// ****************** PROCEDURAL PLANET ************************

float eval_terrain_noise_base(in vec3 npos, const in float gain, const in float lacunarity) {
	float val  = 0.0;
	float mag  = 1.0;
	float freq = 0.5; // lower freq for ridged noise

	for (int i = 0; i < NUM_OCTAVES; ++i) { // similar to gen_cloud_alpha_time()
		float v = simplex(freq*npos);
		v = 2.0*v - 1.0; // map [0,1] range to [-1,1]
		v = max(0.0, (0.75 - abs(v))); // ridged noise
		val  += v*mag;
		freq *= lacunarity;
		mag  *= gain;
	}
	return val;
}

float eval_terrain_noise(in vec3 npos) {
	return 0.7*eval_terrain_noise_base(npos, 0.7, 2.0);
}
float eval_terrain_noise_normal(in vec3 npos) {
	return eval_terrain_noise_base(npos, 0.5, 1.92);
}

// ****************** PLANET DRAWING ************************

float calc_cloud_density(in vec3 lv) {
	float cloud_den = atmosphere*gen_cloud_alpha(lv);
	return pow(clamp(1.4*(cloud_den - 0.1), 0.0, 1.0), 0.7); // increase contrast/sharpen edges
}

vec3 calc_cloud_coord(in vec3 cloud_vertex, in vec3 normal) {
	vec3 lv = 1.2*cloud_freq*cloud_vertex;
#if 0 // too slow
	if (atmosphere > 0.6) { // add swirls
		float v_adj = 0.0;
		float v0    = 1.0;

		for (int i = 0; i < 30; ++i) { // slow when close to the planet
			float dist = (0.5 + 0.25*rand_01(v0+3.0))*length(normal - normalize(rand_vec3(v0)));
			v_adj     += max(0.0, (0.1 - dist))*sin(0.1/max(dist, 1.0));
			v0        += 4.0;
		}
		lv.z += 0.4*v_adj;
	}
#endif
	return lv;
}

vec3 draw_planet(in vec3 vertex, in vec3 normal, in vec3 light_dir, in vec3 sun_pos, in vec3 ss_pos, in float ss_radius) {

	float city_light = 0.0;
	float spec_mag   = 0.0;
	float coldness = pow(abs(normalize(vertex).z), 2.0); // 0 at equator and 1 at the poles
	vec3 spos    = vertex*(terrain_scale/planet_radius);
	vec3 bpos    = 32.0*spos;
	float hval   = eval_terrain_noise(spos);
	float height = max(0.0, 1.8*(hval-0.7)); // can go outside the [0,1] range
	float nscale = 0.0;
	float dnval  = 0.0;
	vec4 texel;

	if (height < water_val) {
		texel    = water_color;
		spec_mag = 1.0;
	}
	else {
		spec_mag = 0.0;
		nscale   = 1.0;
		float height_ws = (height - water_val)/(1.0 - water_val); // rescale to [0,1] above water

		if (water_val > 0.1 && atmosphere > 0.1) { // Earthlike planet
			vec4 gray = vec4(0.4, 0.4, 0.4, 1.0); // gray rock
			if      (height_ws < 0.1) {texel = color_b;} // low ground
			else if (height_ws < 0.4) {texel = mix(color_b, color_a, 3.3333*(height_ws - 0.1));}
			else if (height_ws < 0.5) {texel = color_a;} // medium ground
			else if (height_ws < 1.0) {texel = mix(color_a, gray, 2.0*(height_ws - 0.5));}
			else                      {texel = gray;} // high ground
		}
		else { // alien-like planet
			texel = mix(color_b, color_a, min(1.0, height_ws));
		}
		if (lava_val > 0.0) { // hot lava planet
			if (height < lava_val) {
				texel    = vec4(1,0,0,1); // red
				spec_mag = 0.75;
				nscale   = 0.2;
			}
			else if (height < lava_val + 0.07) {
				float val= (height - lava_val)/0.07;
				texel    = mix(vec4(1,0,0,1), texel, val); // close to lava line
				spec_mag = mix(0.75, spec_mag, val);
				nscale   = mix(0.2,  nscale,   val);
			}
		}
		else if (water_val > 0.0 && temperature < 30.0) { // handle water/ice/snow
			if (height < water_val + 0.07) { // close to water line (can have a little water even if water == 0)
				float val = (height - water_val)/0.07;
				texel     = mix(water_color, texel, val);
				spec_mag  = 1.0 - val;
				nscale    = val*val; // faster falloff
			}
			else if (height_ws > 1.0 && snow_thresh < 1.0) {
				dnval    = eval_terrain_noise_normal(bpos);
				float sv = 0.5 + 0.5*clamp(20.0*(1.0 - snow_thresh), 0.0, 1.0); // snow_thresh 1.0 => no snow, 0.95 => lots of snow
				float mv = dnval * sv * sqrt(height_ws - 1.0);
				spec_mag = 0.5*clamp((1.5*mv*mv - 0.25), 0.0, 1.0);
				texel    = mix(texel, vec4(1,1,1,1), spec_mag); // blend in some snow on peaks
			}
		}
	}
	if (/*snow_thresh < 1.0 &&*/ water_val > 0.2 && temperature < 30.0) { // add polar ice caps
		float icv = 0.7 + 0.01*temperature; // 1.0 @ T=30, 0.9 @ T=20, 0.7 @ T=0
		float val = (coldness - icv)/(1.0 - icv) + 1.0*(height - water_val);
		val       = clamp(3.0*val-1.0, 0.0, 1.0); // sharpen edges
		spec_mag  = mix(spec_mag, 0.7, val);
		texel     = mix(texel, vec4(1,1,1,1), val); // ice/snow
		nscale   *= mix(1.0, 0.25, val);
	}
	vec3 norm = normal;

	if (nscale > 0.0) { // compute normal + bump map
		// Note: using doubles/dvec3 has better precision/quality, but is much slower (what about making them precise?)
		float delta = 0.001;
		float hval0 = ((dnval == 0.0) ? eval_terrain_noise_normal(bpos) : dnval);
		float hdx   = hval0 - eval_terrain_noise_normal(bpos + vec3(delta, 0.0, 0.0));
		float hdy   = hval0 - eval_terrain_noise_normal(bpos + vec3(0.0, delta, 0.0));
		float hdz   = hval0 - eval_terrain_noise_normal(bpos + vec3(0.0, 0.0, delta));
		norm        = normalize(norm) + 1.0*nscale*vec3(hdx, hdy, hdz);
	}
	if (population > 0.0 && spec_mag < 0.5) {
		float thresh = 0.38*population - 0.42*simplex(4.5*spos) - 1.0;
		float freq   = 50.0;

		for (int i = 0; i < 4; ++i) {
			city_light = max(city_light, clamp(4.0*(simplex(freq*spos) + thresh), 0.0, 1.0));
			freq       *= 1.93;
		}
		city_light *= max((0.5 - spec_mag), 0.0) * population; // colonized and not over water/snow/ice
	}
	float atten0 = 1.0;
	float atten2 = 1.0;
	float sscale = atten0;

	if (sun_radius > 0.0 && ss_radius > 0.0) {
		sscale *= calc_sphere_shadow_atten(vertex, sun_pos, sun_radius, ss_pos, ss_radius);
	}
	norm          = normalize(norm); // renormalize
	vec3 ldir0    = light_dir;
	vec3 ldir2    = normalize(ss_pos - vertex);
	vec3 ldir20   = normalize(ss_pos - sun_pos);
	float lscale0 = (dot(norm, ldir0) > 0.0) ? 1.0 : 0.0;
	float lscale2 = (dot(norm, ldir2) > 0.0) ? 1.0 : 0.0;
	float dterm0 = max(dot(norm, ldir0), 0.0);
	float dterm2 = max(dot(norm, ldir2), 0.0);

	// add clouds
	float cloud_den    = 0.0;
	float cloud_shadow = 0.0;
	float cloud_diff   = 1.0;
	vec3 lv = calc_cloud_coord(vertex, normal);

	if (atmosphere > 0.0) {
		cloud_den = calc_cloud_density(lv);

		if (dterm0 > 0.0) {
			float cloud_alt = 0.01*planet_radius; // 1% of planet radius
			vec3 vertex_adj = planet_radius*normalize(vertex + planet_radius*ldir0*cloud_alt/dot(ldir0, vertex)); // approximate
			cloud_shadow = 0.75*calc_cloud_density(calc_cloud_coord(vertex_adj, normal));
			cloud_diff   = 0.8 + 0.2*(1.0 - cloud_shadow);
		}
	}
	float d0 = max(0.0, dot(norm, ldir0));
	float d2 = 0.25*max(0.0, dot(norm, ldir2)) * max(0.0, dot(ldir2, ldir20));
	float diffuse = 0.95*(d0 + d2);
	float ambient = 0.05;
	
	vec3 epos_norm = normalize(vertex);
	vec3 color     = (texel.rgb * (ambient + diffuse*(1.0 - cloud_shadow))); // add light cloud shadows
	vec3 half_vect = normalize(ldir0 - normal); // Eye + L = -eye_space_pos + L
	float specval  = pow(max(dot(norm, half_vect), 0.0), spec_shine);
	color         += ((water_val > 0.0) ? 1.0 : 0.0) * specular_color.rgb * specval * spec_mag * sscale;

	if (lava_val > 0.0) {
		float heat = max(0.0, (texel.r - texel.g - texel.b));

		if (heat > 0.0) { // lava patch
			heat *= gen_cloud_alpha(2.0*vertex);
			color += heat*vec3(1.0, 0.25*heat, 0.0); // add lava
		}
	}
	color += clamp(4.0*city_light*(0.2 - dterm0), 0.0, 1.0)*vec3(1.0, 0.8, 0.5);

	if (cloud_den > 0.0) { // add cloud color
		color = mix(color, vec3(ambient + cloud_diff*diffuse), cloud_den); // no clouds over high mountains?
	}
	return color;
}
