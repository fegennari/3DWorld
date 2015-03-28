uniform sampler2D weights_tex, detail_tex, tex2, tex3, tex4, tex5, tex6, shadow_normal_tex, caustic_tex;
uniform float ts2, ts3, ts4, ts5, ts6; // texture scales
uniform float cs2, cs3, cs4, cs5, cs6; // color scales
uniform float wave_time, wave_amplitude; // water_plane_z comes from water_fog.part
uniform float water_atten    = 1.0;
uniform float normal_z_scale = 1.0;
uniform float spec_scale     = 1.0;
uniform float spec_offset    = 0.0;
uniform float caustics_weight= 1.0;
uniform float smap_atten_cutoff = 10.0;
uniform float smap_atten_slope  = 0.5;
uniform float htex_scale = 1.0;
uniform float triplanar_texture_scale = 1.0;
uniform vec3 uw_atten_max;
uniform vec3 uw_atten_scale;
uniform vec3 snow_cscale = vec3(1.0);

in vec4 vertex; // world space
in vec2 tc2; // for water caustics
//in vec2 tc; // comes from detail_normal_map.part.frag

// underwater attenuation code
void atten_color(inout vec4 color, in float dist) {
	color.rgb *= vec3(1.0) - min(uw_atten_max, mix(vec3(1.0), uw_atten_scale*dist*2.0, htex_scale)); // if not using height texture, treat as max depth
}

float integrate_water_dist(in vec3 targ_pos, in vec3 src_pos, in float water_z) {
	float t = clamp((water_z - targ_pos.z)/max(0.001, abs(src_pos.z - targ_pos.z)), 0.0, 1.0);
	return t*distance(src_pos, targ_pos);
}

vec4 get_light_specular_comp(in vec3 normal, in vec3 light_dir, in vec3 eye_pos, in float spec, in float shininess) {
	vec3 half_vect = normalize(light_dir - normalize(eye_pos)); // Eye + L = -eye_space_pos + L
	return vec4(spec, spec, spec, 1.0) * pow(clamp(dot(normal, half_vect), 0.0001, 1.0), shininess);
}

vec4 add_light_comp(in vec3 normal, in vec4 epos, in int i, in float ds_scale, in float a_scale, in float spec, in float shininess, in float bump_scale) {
	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 light_dir  = normalize(fg_LightSource[i].position.xyz);
	vec3 epos_final = epos.xyz;
#ifdef USE_NORMAL_MAP
	//bump_scale *= pow(texture(detail_tex, 11.0*tc + fract(vec2(0.0, 0.002*wave_time))).r, 3.0); // moving specular when rainy?
	ds_scale *= clamp(5.0*dot(normal, light_dir), 0.0, 1.0); // fix self-shadowing
	normal    = apply_bump_map(light_dir, epos_final, normal, bump_scale);
#if 0 // toksvig antialiasing
	float nmag = min(1.0, length(2.0*texture(detail_normal_tex, detail_normal_tex_scale*tc).rgb - 1.0));
	float ft   = nmag/(nmag + shininess*(1.0 - nmag));
	spec      *= (1.0 + ft*shininess)/(1.0 + shininess);
	shininess *= ft;
#endif
#endif // USE_NORMAL_MAP
	
	// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex
	float NdotL = dot(normal, light_dir);
	vec4 light  = fg_ModelViewMatrixInverse * fg_LightSource[i].position; // world space

	if (apply_cloud_shadows /*&& vertex.z > water_plane_z*//*&& vertex.z < cloud_plane_z*/) { // conditionals are faster but cause seams
		ds_scale *= 1.0 - get_cloud_plane_alpha(vertex.xyz, light);
	}
	
	// compute the ambient and diffuse lighting
	vec4 color = a_scale*fg_LightSource[i].ambient + ds_scale*max(dot(normal, light_dir), 0.0)*fg_LightSource[i].diffuse;
	if (enable_light0) {color += get_light_specular_comp(normal, light_dir, epos_final, ds_scale*spec, shininess);}
	
#ifdef HAS_WATER
	if (vertex.z < water_plane_z) { // underwater
#ifdef WATER_CAUSTICS
		if (i == 0) { // only for light0 (sun)
			// apply underwater caustics texture (Note: matches shallow water wave normal map, but not deep water wave normal map)
			float cweight = ds_scale*wave_amplitude*caustics_weight*min(8.0*(water_plane_z - vertex.z), 0.5);
			float ntime   = 2.0*abs(fract(0.005*wave_time) - 0.5);
			vec3  cval    = 4.0*mix(texture(caustic_tex, tc2).rgb, texture(caustic_tex, (tc2 + vec2(0.3, 0.6))).rgb, ntime);
			color.rgb    *= mix(vec3(1.0), cval, cweight);
		}
#endif // WATER_CAUSTICS
		// apply underwater attenuation
		// Note: ok if vertex is above the water, dist will come out as 0
		vec4 eye    = fg_ModelViewMatrixInverse[3]; // local tile space
		float depth = water_plane_z - vertex.z;
		float dist  = integrate_water_dist(vertex.xyz, eye.xyz, water_plane_z) + min(4.0*depth, integrate_water_dist(vertex.xyz, light.xyz, water_plane_z)); // clamp light pos dir
		atten_color(color, dist*water_atten);
	}
#endif // HAS_WATER
#ifdef GOD_RAYS
	if (apply_cloud_shadows && vertex.z > water_plane_z) {
		color.rgb += get_god_rays(vertex.xyz, fg_ModelViewMatrixInverse[3].xyz, light);
	}
#endif // GOD_RAYS
	return color;
}

vec3 add_texture(in sampler2D tex, in float tc_scale, in vec3 world_n) {
	// separate tc for diffuse texture, in case we want to sometimes mirror it to make tiling less periodic (though seems difficult and unnecessary)
#ifdef TRIPLANAR_TEXTURE
	return lookup_triplanar_texture_scaled(vertex.xyz, world_n, tex, tex, tex, triplanar_texture_scale*tc_scale).rgb;
#else
	return texture(tex, tc_scale*tc).rgb;
#endif
}

void main()
{
#ifdef REFLECTION_MODE
	if (vertex.z < water_plane_z) {discard;}
#endif
	vec4 shadow_normal  = texture(shadow_normal_tex, tc);
	float diffuse_scale = shadow_normal.w;
	float ambient_scale = 1.5*shadow_normal.z;
	vec2 nxy    = (2.0*shadow_normal.xy - 1.0);
	vec3 normal = vec3(nxy, normal_z_scale*(1.0 - sqrt(nxy.x*nxy.x + nxy.y*nxy.y))); // calculate n.z from n.x and n.y (we know it's always positive)
	vec3 world_n= normalize(normal);
	normal      = normalize(fg_NormalMatrix * normal); // eye space

	// sand, dirt, grass, rock, snow
	vec4 weights   = texture(weights_tex, tc);
	float weights4 = clamp((1.0 - weights.r - weights.g - weights.b - weights.a), 0.0, 1.0);
	weights  = smoothstep(0.0, 1.0, weights);
	weights4 = smoothstep(0.0, 1.0, weights4);

	vec3 texel0  = vec3(0.0);
	if (weights.r > 0) {texel0 += cs2*weights.r*add_texture(tex2, ts2, world_n);} // sand
	if (weights.g > 0) {texel0 += cs3*weights.g*add_texture(tex3, ts3, world_n);} // dirt
	if (weights.b > 0) {texel0 += cs4*weights.b*add_texture(tex4, ts4, world_n);} // grass
	if (weights.a > 0) {texel0 += cs5*weights.a*add_texture(tex5, ts5, world_n);} // rock
	if (weights4  > 0) {texel0 += cs6*weights4 *add_texture(tex6, ts6, world_n)*snow_cscale;} // snow
	vec3 texel1 = add_texture(detail_tex, 32.0, world_n).rgb; // detail texture 32x scale

	vec4 color  = vec4(0,0,0,1);
	vec4 epos   = fg_ModelViewMatrix * vertex;
	float vdist = length(epos.xyz);
	float bump_scale = 1.0 - weights.b; // bumps on everything but grass
	bump_scale *= clamp((2.5 - 0.1*vdist), 0.0, 1.0); // decrease scale with distance to reduce tiling artifacts on sand and snow
	float smap_scale = 0.0;
	if (use_shadow_map) {smap_scale = clamp(smap_atten_slope*(smap_atten_cutoff - vdist), 0.0, 1.0);}
	
	if (enable_light0) { // sun
		float spec      = spec_scale*(spec_offset + 0.2*weights.b + 0.25*weights4); // grass and snow
		float shininess = 80.0*spec_offset + 20.0*weights.b + 40.0*weights4;
		float dscale    = diffuse_scale;
		if (use_shadow_map) {dscale = min(dscale, mix(1.0, get_shadow_map_weight_light0(epos, normal), smap_scale));}
		color += add_light_comp(normal, epos, 0, dscale, ambient_scale, spec, shininess, bump_scale);
	}
	if (enable_light1) { // moon
		float dscale    = diffuse_scale;
		if (use_shadow_map) {dscale = min(dscale, mix(1.0, get_shadow_map_weight_light1(epos, normal), smap_scale));}
		color += add_light_comp(normal, epos, 1, dscale, ambient_scale, 0.0, 1.0, bump_scale);
	}
	if (enable_light2) {color += add_light_comp(normal, epos, 2, 1.0, 1.0, 0.0, 1.0, bump_scale) * calc_light_atten(epos, 2);} // lightning

	vec4 mesh_color = vec4((texel0.rgb * texel1.rgb * color.rgb), color.a);
	vec3 clipped_vert;
	float fog_coord = get_water_fog_coord(vertex, fg_ModelViewMatrixInverse[3].xyz, clipped_vert)*get_custom_fog_scale(clipped_vert.z);
	fg_FragColor    = apply_fog_ffc(mesh_color, fog_coord, fog_color);
	//fg_FragColor    = color; // untextured (white) for debugging
}
