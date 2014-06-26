// Note: Light 0 is the sun (A+D+S point light), light 1 is universe ambient (constant A), light 2 is planet reflection (D point light)
uniform float atmosphere = 1.0; // technically not needed for gas giants since assumed to be 1.0
uniform vec3 cloud_freq  = vec3(1.0);
uniform vec3 light_scale = vec3(1.0);
uniform vec4 emission    = vec4(0,0,0,1);
uniform vec3 sun_pos, ss_pos, rscale;
uniform float sun_radius, ss_radius, ring_ri, ring_ro;
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
uniform float obj_radius, snow_thresh, cold_scale, noise_offset, terrain_scale, temperature;
uniform vec4 water_color, color_a, color_b;
#endif // PROCEDURAL_DETAIL
#endif // not GAS_GIANT

varying vec3 normal, world_space_pos, vertex;
varying vec2 tc;


void main()
{
	vec4 epos = fg_ModelViewMatrix * vec4(vertex, 1.0);
	if (dot(normal, epos.xyz) > 0.0) discard; // back facing
	vec3 norm      = normal;
	float spec_mag = 0.0;

#ifdef GAS_GIANT
	float tc_adj = tc.t + 0.04*(gen_cloud_alpha_static_non_norm(5.0*vertex) - 0.5);
	float v0 = 1.0; // using a variable here is slow
	vec3 dir = normalize(vertex); // world space normal

	for (int i = 0; i < 50; ++i) { // Note: inefficient, but fast enough for a single gas giant render
		vec3 center = vec3(1.0, 1.0, 0.5)*rand_vec3(v0);
#ifdef ANIMATE_STORMS // slow but neat
		float angle = 50.0*time*rand_pm1(v0+2.5);
		float st    = sin(angle);
		float ct    = cos(angle);
		center.xy   = vec2((center.x*ct - center.y*st), (center.y*ct + center.x*st)); // rotation
#endif // ANIMATE_STORMS
		float dist  = (0.25 + 0.75*rand_01(v0+3.0))*length(vec3(1.0, 1.0, 2.0)*(dir - normalize(center)));
		tc_adj     += 0.5*max(0.0, (0.1 - dist))*sin(0.1/max(dist, 0.01));
		v0         += 4.0;
	}
	vec4 texel   = texture1D(tex0, tc_adj);
#else // not GAS_GIANT

#ifdef PROCEDURAL_DETAIL
	float coldness = cold_scale*pow(abs(normalize(vertex).z), 2.0); // 0 at equator and 1 at the poles
#ifdef ALL_WATER_ICE
	vec4 texel = water_color;
	if (coldness > 0.75) {texel = mix(texel, vec4(1,1,1,1), clamp(4.0*(coldness - 0.75f), 0.0, 1.0));} // ice/snow
	spec_mag = 1.0; // always specular
#else // not ALL_WATER_ICE
	float hval = 0.0;
	float freq = 1.0;
	vec3 npos  = vertex*(terrain_scale/obj_radius) + vec3(noise_offset);

	for (int i = 0; i < 8; ++i) { // similar to gen_cloud_alpha_time()
		hval += texture3D(cloud_noise_tex, freq*npos).r/freq;
		freq *= 2.0;
	}
	float height = max(0.0, 1.8*(hval-0.7)); // can go outside the [0,1] range
	vec4 texel;

	if (height < water_val) {
		texel    = water_color;
		spec_mag = 1.0;
	}
	else {
		spec_mag = 0.0;
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
			if      (height < lava_val)        {texel = vec4(1,0,0,1);} // red
			else if (height < lava_val + 0.07) {texel = mix(vec4(1,0,0,1), texel, (height - lava_val)/0.07);} // close to lava line
		}
		else if (water_val > 0.0 && temperature < 30.0) { // handle water/ice/snow
			if (height < water_val + 0.07) { // close to water line (can have a little water even if water == 0)
				float val = (height - water_val)/0.07;
				texel     = mix(water_color, texel, val);
				spec_mag  = 1.0 - val;
			}
			else if (height_ws > 1.0 && snow_thresh < 1.0) {
				float freq = 1.0;
				float val  = 0.0;

				for (int i = 0; i < 4; ++i) { // scample high frequencies from the above iteration?
					val  += texture3D(cloud_noise_tex, 32.0*freq*npos).r/freq;
					freq *= 2.0;
				}
				float sv = 0.5 + 0.5*clamp(20.0*(1.0 - snow_thresh), 0.0, 1.0); // snow_thresh 1.0 => no snow, 0.95 => lots of snow
				float mv = val * sv * sqrt(height_ws - 1.0);
				spec_mag = clamp((1.5*mv*mv - 0.25), 0.0, 1.0);
				texel    = mix(texel, vec4(1,1,1,1), spec_mag); // blend in some snow on peaks
			}
		}
	}
	if (/*snow_thresh < 1.0 &&*/ water_val > 0.2 && temperature < 30.0) { // add polar ice caps
		float icv = 0.7 + 0.01*temperature; // 1.0 @ T=30, 0.9 @ T=20, 0.7 @ T=0
		float val = (coldness - icv)/(1.0 - icv) + 1.0*(height - water_val);
		val       = clamp(3*val-1, 0.0, 1.0); // sharpen edges
		spec_mag  = mix(spec_mag, 1.0, val);
		texel     = mix(texel, vec4(1,1,1,1), val); // ice/snow
	}
	norm = fg_NormalMatrix * vertex; // recompute
	// FIXME: perturb the normal by looking at derivative of normal at this point
#endif // ALL_WATER_ICE

#else
	vec4 texel = texture2D(tex0, tc);
	spec_mag   = pow(texel.b, 4.0);
#endif // PROCEDURAL_DETAIL

#endif // GAS_GIANT

	float atten0 = light_scale[0] * calc_light_atten0(epos);
	float atten2 = light_scale[2] * calc_light_atten(epos, 2);
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
					sscale      *= 1.0 - dscale*texture1D(ring_tex, rval).a;
				}
			}
		}
	}
	norm          = normalize(norm); // renormalize
	vec3 ldir0    = normalize(fg_LightSource[0].position.xyz - epos.xyz);
	vec3 ldir2    = normalize(fg_LightSource[2].position.xyz - epos.xyz);
	vec3 ldir20   = normalize(fg_LightSource[2].position.xyz - fg_LightSource[0].position.xyz);
	float lscale0 = (dot(norm, ldir0) > 0.0) ? 1.0 : 0.0;
	float lscale2 = (dot(norm, ldir2) > 0.0) ? 1.0 : 0.0;

#ifdef HAS_CRATERS
	// facing the sun or planet (reflected light), and not over water (blue)
	if ((lscale0 > 0.0 || lscale2 > 0.0) && (texel.b - texel.r - texel.g) < 0.0) { // smoother transition?
		adjust_normal_for_craters(norm, vertex); // add craters by modifying the normal
	}
#endif // HAS_CRATERS

	vec3 epos_norm = normalize(epos.xyz);
	vec3 ambient   = (fg_LightSource[0].ambient.rgb * atten0) + (fg_LightSource[1].ambient.rgb * light_scale[1]);
	vec3 diffuse   = (fg_LightSource[0].diffuse.rgb * max(dot(norm, ldir0), 0.0) * lscale0 * sscale) +
	                 (fg_LightSource[2].diffuse.rgb * max(dot(norm, ldir2), 0.0) * lscale2 * atten2 * max(dot(ldir2, ldir20), 0.0));
	vec3 color     = (texel.rgb * (ambient + diffuse));

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
#endif // GAS_GIANT
	if (atmosphere > 0.0) {
		float cloud_val = atmosphere*gen_cloud_alpha(cloud_freq*vertex);
		if (cloud_val > 0.0) {color = cloud_val*(ambient + diffuse) + (1.0 - cloud_val)*color;} // no clouds over high mountains?
	}
	fg_FragColor = gl_Color * vec4((color + emission.rgb), 1.0);
}
