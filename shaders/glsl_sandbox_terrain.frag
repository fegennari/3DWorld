// 

// Terrain rendering
// Written by Frank Gennari
// 10/13/14

#ifdef GL_ES
precision mediump float;
#endif

// standard uniforms
uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;

// constants to be changed by the user
const float water_val     = 0.4;
const float snow_thresh   = 1.0;
const float terrain_scale = 1.0;
const float height_bias   = 0.0;
const float noise_scale   = 1.0;
const float cloud_scale   = 4.0;
const float cloud_density = 0.75;
const float cloud_height  = 0.04;
const vec2 cloud_speed    = vec2(0.03, 0.05);
const vec2 translate      = vec2(1.5, 1.0);
const vec4 water_color    = vec4(0.0, 0.1, 0.6, 1.0); // blue-ish
const vec4 color_a        = vec4(0.05, 0.35, 0.05, 1.0); // dark green vegetation
const vec4 color_b        = vec4(0.60, 0.45, 0.25, 1.0); // brown dirt
#define NUM_OCTAVES 8

vec3 draw_terrain(in vec2 vertex, in vec3 light_dir);

void main(void) {
	float rmin = min(resolution.x, resolution.y);
	vec2 position  = ((gl_FragCoord.xy - 0.5*resolution.xy) / rmin) + translate;
	//position.x    += 0.05*time; // scroll the terrain
	vec3 light_dir = normalize(vec3(2.0*mouse.x-1.0, 2.0*mouse.y-1.0, 1.0));
	gl_FragColor   = vec4(draw_terrain(position, light_dir), 1.0);
}

// ****************** SIMPLEX NOISE ************************

vec3 mod289(in vec3 x) {
	return x - floor(x * 1.0 / 289.0) * 289.0;
}
vec3 permute(in vec3 x) {
	return mod289(((x * 34.0) + 1.0) * x);
}

float hash( vec2 p ) {
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
#if 1 // faster but lower quality
	return noise(4.0*v);
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

// ****************** PERLIN CLOUDS 3D ************************

float gen_cloud_alpha_time(in vec2 pos, in vec2 ftime)
{
	float alpha = 0.0;
	float freq  = 1.0;

	for (int i = 0; i < NUM_OCTAVES; ++i) {
		float v = simplex(noise_scale*(freq*pos + ftime));
		//float v = noise(4.0*noise_scale*(freq*pos + ftime));
		alpha += v/freq;
		freq  *= 2.0;
	}
	return alpha;
}

float fract_smooth(in float t) {
	return 2.0*abs(fract(0.5*t) - 0.5);
}
vec2 get_ctime() {
	return time*cloud_speed;
}
float calc_cloud_density(in vec2 pos) {
	float density = gen_cloud_alpha_time(cloud_scale*pos, get_ctime());
	//return pow(clamp(1.4*(density - 0.1), 0.0, 1.0), 0.7); // increase contrast/sharpen edges
	return clamp(cloud_density*density, 0.0, 1.0);
}

// ****************** PROCEDURAL TERRAIN ************************

float eval_terrain_noise_base(in vec2 npos, const in float gain, const in float lacunarity) {
	float val  = 0.0;
	float mag  = 1.0;
	float freq = 0.5; // lower freq for ridged noise

	for (int i = 0; i < NUM_OCTAVES; ++i) { // similar to gen_cloud_alpha_time()
		float v = simplex(freq*npos);
		v = 0.9 - abs(v); // ridged noise
		val  += v*mag;
		freq *= lacunarity;
		mag  *= gain;
	}
	return val + height_bias;
}

float eval_terrain_noise(in vec2 npos) {
	return 0.7*eval_terrain_noise_base(npos, 0.7, 2.0);
}
float eval_terrain_noise_detail(in vec2 npos) {
	return eval_terrain_noise_base(npos, 0.5, 1.92);
}

vec3 draw_terrain(in vec2 vertex, in vec3 light_dir) {

	vec2 spos    = vertex*terrain_scale;
	float hval   = eval_terrain_noise(spos);
	float height = max(0.0, 1.8*(hval-0.7)); // can go outside the [0,1] range
	float nscale = 0.0;
	vec4 texel;

	if (height < water_val) {
		texel = water_color;
	}
	else { // Earthlike planet
		nscale = 1.0;
		float height_ws = (height - water_val)/(1.0 - water_val); // rescale to [0,1] above water
		vec4 gray = vec4(0.4, 0.4, 0.4, 1.0); // gray rock
		if      (height_ws < 0.1) {texel = color_b;} // low ground
		else if (height_ws < 0.4) {texel = mix(color_b, color_a, 3.3333*(height_ws - 0.1));}
		else if (height_ws < 0.5) {texel = color_a;} // medium ground
		else if (height_ws < 1.0) {texel = mix(color_a, gray, 2.0*(height_ws - 0.5));}
		else                      {texel = gray;} // high ground

		if (water_val > 0.0) { // handle water
			if (height < water_val + 0.07) { // close to water line (can have a little water even if water == 0)
				float val = (height - water_val)/0.07;
				texel     = mix(water_color, texel, val);
				nscale    = val*val; // faster falloff
			}
			else if (height_ws > 1.0) {
				float mv = eval_terrain_noise_detail(32.0*spos) * snow_thresh * sqrt(height_ws - 1.0);
				float mag= 0.5*clamp((1.5*mv*mv - 0.25), 0.0, 1.0);
				texel    = mix(texel, vec4(1,1,1,1), mag); // blend in some snow on peaks
			}
		}
	}

	// compute normal
	float delta = 0.001;
	float hdx   = hval - eval_terrain_noise(spos + vec2(delta, 0.0));
	float hdy   = hval - eval_terrain_noise(spos + vec2(0.0, delta));
	vec3 norm   = normalize(mix(vec3(0.0, 0.0, 1.0), normalize(vec3(hdx, hdy, 0.1)), nscale));

	// cloud shadows
	vec2 cloud_pos = vertex + cloud_height*light_dir.xy;
	float cloud_shadow = calc_cloud_density(cloud_pos);
	norm *= 1.0 - cloud_shadow; // attenuate diffuse color through the normal
	
	// lighting
	float diffuse = 0.9*max(dot(norm, light_dir), 0.0);
	float ambient = 0.1;
	vec3 color = texel.rgb * (ambient + diffuse);

	// clouds
	vec3 cloud_color = vec3(1.0, 1.0, 1.0);
	return mix(color, cloud_color, calc_cloud_density(vertex));
}
