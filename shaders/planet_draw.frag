// Note: Light 0 is the sun (A+D+S point light), light 1 is universe ambient (constant A), light 2 is planet reflection (D point light)
uniform float atmosphere = 1.0; // technically not needed for gas giants since assumed to be 1.0
uniform vec3 cloud_freq  = vec3(1.0);
uniform vec3 light_scale = vec3(1.0);
uniform vec4 emission    = vec4(0,0,0,1);
uniform vec3 sun_pos, ss_pos, rscale;
uniform float obj_radius, sun_radius, ss_radius, ring_ri, ring_ro, population;
uniform mat4 fg_ViewMatrix;
uniform sampler1D ring_tex;

#ifdef GAS_GIANT
uniform sampler1D tex0;
#else
uniform float water_val  = 0.0;
uniform float lava_val   = 0.0;
uniform float crater_val = 0.0;
uniform sampler2D tex0;
#ifdef PROCEDURAL_DETAIL
const int NORMAL_OCTAVES = 6;
uniform float snow_thresh, cold_scale, temperature, nmap_mag;
uniform vec4 water_color, color_a, color_b;
#endif // PROCEDURAL_DETAIL
#endif // not GAS_GIANT

in vec3 normal, world_space_pos, vertex;
in vec2 tc;


float calc_cloud_density(in vec3 lv) {
	float cloud_den = atmosphere*gen_cloud_alpha(lv);
	return pow(clamp(1.4*(cloud_den - 0.1), 0.0, 1.0), 0.7); // increase contrast/sharpen edges
}

vec3 calc_cloud_coord(in vec3 cloud_vertex) {
	vec3 lv = 1.2*cloud_freq*cloud_vertex;
#ifndef GAS_GIANT
	if (atmosphere > 0.6) { // add swirls
		float v_adj = 0.0;
		float v0    = 1.0;
		vec3 dir    = normalize(vertex); // world space normal

		for (int i = 0; i < 30; ++i) { // slow when close to the planet
			float dist = (0.5 + 0.25*rand_01(v0+3.0))*length(dir - normalize(rand_vec3(v0)));
			v_adj     += max(0.0, (0.1 - dist))*sin(0.1/max(dist, 1.0));
			v0        += 4.0;
		}
		lv.z += 0.4*v_adj;
	}
#endif // GAS_GIANT
	return lv;
}

void main()
{
	vec4 epos = fg_ModelViewMatrix * vec4(vertex, 1.0);
	vec3 norm        = normal;
	float city_light = 0.0;
	float spec_mag   = 0.0;

#ifdef GAS_GIANT
	float tc_adj = tc.t;
	float noise  = gen_cloud_alpha_static_non_norm(5.0*vertex);
	vec3 dir     = normalize(vertex); // world space normal
	tc_adj      += 0.2*sin(20.0*dir.z + 1.0*noise);
	//tc_adj      += 0.04*(noise - 0.5);
	float v0     = 1.0; // using a variable here is slow

	for (int i = 0; i < 50; ++i) { // Note: inefficient, but fast enough for a single gas giant render
		vec3 center = vec3(1.0, 1.0, 0.5)*rand_vec3(v0);
#ifdef ANIMATE_STORMS // slow but neat
		float angle = 50.0*time*rand_pm1(v0+2.5);
		float st    = sin(angle);
		float ct    = cos(angle);
		center.xy   = vec2((center.x*ct - center.y*st), (center.y*ct + center.x*st)); // rotation
#endif // ANIMATE_STORMS
		float dist  = 1.0*(0.25 + 0.75*rand_01(v0+3.0))*length(vec3(1.0, 1.0, 2.0)*(dir - normalize(center)));
		tc_adj     += 1.5*max(0.0, (0.1 - dist))*sin(0.1/max(dist, 0.01));
		v0         += 4.0;
	}
	vec4 texel = texture(tex0, tc_adj);
#else // not GAS_GIANT

#ifdef PROCEDURAL_DETAIL
	float coldness = cold_scale*pow(abs(normalize(vertex).z), 2.0); // 0 at equator and 1 at the poles
#ifdef ALL_WATER_ICE
	vec4 texel = water_color;
	if (coldness > 0.75) {texel = mix(texel, vec4(1,1,1,1), clamp(4.0*(coldness - 0.75f), 0.0, 1.0));} // ice/snow
	spec_mag = 1.0; // always specular
#else // not ALL_WATER_ICE

	vec3 spos    = vertex*(terrain_scale/obj_radius);
	vec3 npos    = spos + vec3(noise_offset);
	vec3 bpos    = 32.0*spos;
	float hval   = eval_terrain_noise(npos, 8);
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

		if (water_val > 0.2 && atmosphere > 0.1) { // Earthlike planet
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
				dnval    = eval_terrain_noise_normal(bpos, NORMAL_OCTAVES);
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
	nscale *= nmap_mag;

	if (nscale > 0.0) { // compute normal + bump map
		// Note: using doubles/dvec3 has better precision/quality, but is much slower (what about making them precise?)
		float delta = 0.001;
		float hval0 = ((dnval == 0.0) ? eval_terrain_noise_normal(bpos, NORMAL_OCTAVES) : dnval);
		float hdx   = hval0 - eval_terrain_noise_normal(bpos + vec3(delta, 0.0, 0.0), NORMAL_OCTAVES);
		float hdy   = hval0 - eval_terrain_noise_normal(bpos + vec3(0.0, delta, 0.0), NORMAL_OCTAVES);
		float hdz   = hval0 - eval_terrain_noise_normal(bpos + vec3(0.0, 0.0, delta), NORMAL_OCTAVES);
		norm = normalize(norm) + 1.0*nscale*(fg_NormalMatrix * vec3(hdx, hdy, hdz));
	}
	if (population > 0.0 && spec_mag < 0.5) {
		float thresh = 0.38*population - 0.42*texture(cloud_noise_tex, 4.5*spos).r - 1.0;
		float freq   = 50.0;

		for (int i = 0; i < 4; ++i) {
			city_light = max(city_light, clamp(4.0*(texture(cloud_noise_tex, freq*spos).r + thresh), 0.0, 1.0));
			freq       *= 1.93;
		}
		city_light *= max((0.5 - spec_mag), 0.0) * population; // colonized and not over water/snow/ice
	}
#endif // ALL_WATER_ICE

#else // not PROCEDURAL_DETAIL
	vec4 texel = texture(tex0, tc);
	spec_mag   = pow(texel.b, 4.0);
#endif // PROCEDURAL_DETAIL

#endif // GAS_GIANT

	float atten0 = light_scale[0] * calc_light_atten0(epos);
	float sscale = atten0;

	if (sun_radius > 0.0) {
		if (ss_radius > 0.0) {
			sscale *= calc_sphere_shadow_atten(world_space_pos, sun_pos, sun_radius, ss_pos, ss_radius);
		}
		if (has_rings) { // calculate shadows due to rings
			vec3 sun_local = (fg_ModelViewMatrixInverse * (fg_ViewMatrix * vec4(sun_pos, 1.0))).xyz;
			vec3 line_dir  = sun_local - vertex;
			float dist     = -vertex.z/line_dir.z; // Note: ring normal is always in z

			if (dist > 0.0) {
				vec3 ring_ipt = dist*line_dir + vertex;
				float rval    = (length(ring_ipt/rscale) - ring_ri)/(ring_ro - ring_ri);
				
				if (rval > 0.0 && rval < 1.0) {
					float dscale = length(vertex)/length(ring_ipt); // fake penumbra
					sscale      *= 1.0 - dscale*texture(ring_tex, rval).a;
				}
			}
		}
	}
	norm          = normalize(norm); // renormalize
	vec3 ldir0    = normalize(fg_LightSource[0].position.xyz - epos.xyz);
	vec3 ldir2    = normalize(fg_LightSource[2].position.xyz - epos.xyz);
	float lscale0 = (dot(norm, ldir0) > 0.0) ? 1.0 : 0.0;
	float lscale2 = (dot(norm, ldir2) > 0.0) ? 1.0 : 0.0;

#ifdef HAS_CRATERS
	// facing the sun or planet (reflected light), and not over water (blue)
	if ((lscale0 > 0.0 || lscale2 > 0.0) && (texel.b - texel.r - texel.g) < 0.0) { // smoother transition?
		adjust_normal_for_craters(norm, vertex); // add craters by modifying the normal
	}
#endif // HAS_CRATERS
	float dterm0 = max(dot(norm, ldir0), 0.0);

	// add clouds
	float cloud_den    = 0.0;
	float cloud_shadow = 0.0;
	float cloud_diff   = 1.0;
	vec3 lv = calc_cloud_coord(vertex);

	if (atmosphere > 0.0) {
		cloud_den = calc_cloud_density(lv);
#ifndef NO_CLOUD_SHADOWS
		if (dterm0 > 0.0) {
			float cloud_alt = 0.01*obj_radius; // 1% of planet radius
			vec3 obj_space_ldir = inverse(fg_NormalMatrix) * ldir0; // no normalization needed
			vec3 vertex_adj = obj_radius*normalize(vertex + obj_radius*obj_space_ldir*cloud_alt/dot(obj_space_ldir, vertex)); // approximate
			cloud_shadow = 0.75*calc_cloud_density(calc_cloud_coord(vertex_adj));
			cloud_diff   = 0.8 + 0.2*(1.0 - cloud_shadow);
		}
#else
		cloud_shadow = 0.25*cloud_den;
#endif
	}
	vec3 epos_norm = normalize(epos.xyz);
	vec3 ambient   = (fg_LightSource[0].ambient.rgb * atten0) + (fg_LightSource[1].ambient.rgb * light_scale[1]);
	vec3 diffuse   = (fg_LightSource[0].diffuse.rgb * dterm0 * lscale0 * sscale);
	
	if (light_scale[2] > 0.0) {
		float dterm2 = max(dot(norm, ldir2), 0.0);
		vec3 ldir20  = normalize(fg_LightSource[2].position.xyz - fg_LightSource[0].position.xyz);
		diffuse += (fg_LightSource[2].diffuse.rgb * dterm2 * lscale2 * light_scale[2] * calc_light_atten(epos, 2) * max(dot(ldir2, ldir20), 0.0));
	}
	vec3 color = (texel.rgb * (ambient + diffuse*(1.0 - cloud_shadow))); // add light cloud shadows

#ifndef GAS_GIANT
	vec3 half_vect = normalize(ldir0 - epos_norm); // Eye + L = -eye_space_pos + L
	float specval  = pow(max(dot(norm, half_vect), 0.0), get_shininess());
	color         += ((water_val > 0.0) ? 1.0 : 0.0) * fg_LightSource[0].specular.rgb*specular_color.rgb * specval * spec_mag * sscale;

	if (lava_val > 0.0) {
		float heat = max(0.0, (texel.r - texel.g - texel.b));

		if (heat > 0.0) { // lava patch
			heat *= gen_cloud_alpha(2.0*vertex);
			//color = mix(color, vec3(1.0, 0.25*heat, 0.0), heat); // add lava
			color += heat*vec3(1.0, 0.25*heat, 0.0); // add lava
		}
	}
#endif // not GAS_GIANT
	color += emission.rgb + clamp(4.0*city_light*(0.2 - dterm0), 0.0, 1.0)*vec3(1.0, 0.8, 0.5);

	if (cloud_den > 0.0) { // add cloud color
		color = mix(color, (ambient + cloud_diff*diffuse), cloud_den); // no clouds over high mountains?
	}
	if (cloud_den > 0.5) { // maybe add lightning
		float v0   = float(int(4.0E4*time))/4.0E4;
		float dist = (10.0 + 10.0*rand_01(v0+3.0))*length(normalize(vertex) - normalize(rand_vec3(v0)));
		float val  = max(0.0, (1.0 - dist));
		color     += 2.5*(cloud_den - 0.5)*val*vec3(0.6, 0.8, 1.0); // lightning color
	}
	fg_FragColor = gl_Color * vec4(color, 1.0);
}
