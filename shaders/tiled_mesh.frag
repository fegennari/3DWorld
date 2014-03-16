uniform sampler2D weights_tex, detail_tex, tex2, tex3, tex4, tex5, tex6, shadow_normal_tex, noise_tex, caustic_tex, detail_normal_tex;
uniform float ts2, ts3, ts4, ts5, ts6; // texture scales
uniform float cs2, cs3, cs4, cs5, cs6; // color scales
uniform float water_plane_z, cloud_plane_z, wave_time, wave_amplitude;
uniform float water_atten    = 1.0;
uniform float normal_z_scale = 1.0;
uniform float spec_scale     = 1.0;
uniform float cloud_alpha    = 1.0;
uniform float caustics_weight= 1.0;
uniform vec3 cloud_offset    = vec3(0.0);
uniform vec3 uw_atten_max;
uniform vec3 uw_atten_scale;
uniform vec3 snow_cscale = vec3(1.0);

varying vec4 vertex; // world space
varying vec2 tc, tc2, tc3;

// underwater attenuation code
void atten_color(inout vec4 color, in float dist) {
	color.rgb *= vec3(1.0) - min(uw_atten_max, uw_atten_scale*dist*2.0);
}

float integrate_water_dist(in vec3 targ_pos, in vec3 src_pos, in float water_z) {
	float t = clamp((water_z - targ_pos.z)/max(0.001, abs(src_pos.z - targ_pos.z)), 0.0, 1.0);
	return t*distance(src_pos, targ_pos);
}

vec4 get_light_specular(in vec3 normal, in vec3 light_dir, in vec3 eye_pos, in float spec, in float shininess) {
	vec3 half_vect = normalize(light_dir - normalize(eye_pos)); // Eye + L = -eye_space_pos + L
	return vec4(spec, spec, spec, 1.0) * pow(max(dot(normal, half_vect), 0.0), shininess);
}

vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos, in vec3 normal, in float bump_scale) {
	vec3 tan  = normalize(cross(gl_ModelViewMatrix[0].xyz, normal));
	mat3 TBN  = transpose(mat3(tan, cross(normal, tan), normalize(normal))); // world space {Y, X, Z} for normal in +Z
	light_dir = TBN * light_dir;
	eye_pos   = TBN * eye_pos;
	return normalize(mix(vec3(0,0,1), (2.0*texture2D(detail_normal_tex, 0.25*tc2).rgb - 1.0), bump_scale)); // scaled detail texture tc
}

vec4 add_light_comp(in vec3 normal, in vec4 epos, in int i, in float ds_scale, in float a_scale, in float spec, in float shininess, in float bump_scale) {
	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 light_dir = normalize(gl_LightSource[i].position.xyz);
#ifdef USE_NORMAL_MAP
	ds_scale *= clamp(5.0*dot(normal, light_dir), 0.0, 1.0); // fix self-shadowing
	normal    = apply_bump_map(light_dir, epos, normal, bump_scale);
#endif
		
	// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex
	float NdotL = dot(normal, light_dir);
	vec4 light  = gl_ModelViewMatrixInverse * gl_LightSource[i].position; // world space

	if (apply_cloud_shadows /*&& vertex.z > water_plane_z*//*&& vertex.z < cloud_plane_z*/) {
		vec3 cpos = vertex.xyz + cloud_offset;
		float t = clamp((cloud_plane_z - cpos.z)/(light.z - cpos.z), 0.0, 1.0); // sky intersection position along vertex->light vector
		normal *= 1.0 - cloud_alpha*gen_cloud_alpha(cpos.xy + t*(light.xy - cpos.xy));
	}
	
	// compute the ambient and diffuse lighting
	vec4 color = a_scale*gl_LightSource[i].ambient + ds_scale*max(dot(normal, light_dir), 0.0)*gl_LightSource[i].diffuse;
	if (enable_light0) {color += get_light_specular(normal, light_dir, epos.xyz, ds_scale*spec, shininess);}
	
#ifdef HAS_WATER
	if (vertex.z < water_plane_z) { // underwater
#ifdef WATER_CAUSTICS
		if (i == 0) { // only for light0 (sun)
			// apply underwater caustics texture (Note: matches shallow water wave normal map, but not deep water wave normal map)
			float cweight = ds_scale*wave_amplitude*caustics_weight*min(8.0*(water_plane_z - vertex.z), 0.5);
			float ntime   = 2.0*abs(fract(0.005*wave_time) - 0.5);
			vec3  cval    = 4.0*mix(texture2D(caustic_tex, tc3).rgb, texture2D(caustic_tex, (tc3 + vec2(0.3, 0.6))).rgb, ntime);
			color.rgb    *= mix(vec3(1.0), cval, cweight);
		}
#endif
		// apply underwater attenuation
		// Note: ok if vertex is above the water, dist will come out as 0
		vec4 eye    = gl_ModelViewMatrixInverse[3]; // world space
		float depth = water_plane_z - vertex.z;
		float dist  = integrate_water_dist(vertex.xyz, eye.xyz, water_plane_z) + min(4.0*depth, integrate_water_dist(vertex.xyz, light.xyz, water_plane_z)); // clamp light pos dir
		atten_color(color, dist*water_atten);
	}
#endif
	return color;
}

void main()
{
#ifdef REFLECTION_MODE
	if (vertex.z < water_plane_z) {discard;}
#endif
	// sand, dirt, grass, rock, snow
	vec2 diff_tc = tc; // separate tc for diffuse texture, in case we want to sometimes mirror it to make tiling less periodic (though seems difficult and unnecessary)
	//diff_tc.s += 0.1*vertex.z; // we really need something like triplanar texturing here to deal with stretching on steep slopes
	vec4 weights = texture2D(weights_tex, tc);
	float weights4 = clamp((1.0 - weights.r - weights.g - weights.b - weights.a), 0.0, 1.0);
	vec3 texel0  = cs2*weights.r*texture2D(tex2, ts2*diff_tc).rgb + // sand
	               cs3*weights.g*texture2D(tex3, ts3*diff_tc).rgb + // dirt
				   cs4*weights.b*texture2D(tex4, ts4*diff_tc).rgb + // grass
				   cs5*weights.a*texture2D(tex5, ts5*diff_tc).rgb + // rock
				   cs6*weights4 *texture2D(tex6, ts6*diff_tc).rgb*snow_cscale;  // snow
	vec3 texel1  = texture2D(detail_tex, tc2).rgb; // detail texture

	vec4 shadow_normal  = texture2D(shadow_normal_tex, tc);
	float diffuse_scale = shadow_normal.w;
	float ambient_scale = 1.5*shadow_normal.z;
	float bump_scale    = 1.0 - weights.b; // bumps on everything but grass
	vec2 nxy    = (2.0*shadow_normal.xy - 1.0);
	vec3 normal = vec3(nxy, normal_z_scale*(1.0 - sqrt(nxy.x*nxy.x + nxy.y*nxy.y))); // calculate n.z from n.x and n.y (we know it's always positive)
	normal      = normalize(gl_NormalMatrix * normal); // eye space
	//normal     += 0.05*weights4*vec3(texture2D(noise_tex, 571.0*tc).r-0.5, texture2D(noise_tex, 714.0*tc).r-0.5, texture2D(noise_tex, 863.0*tc).r-0.5); // add noise
	vec4 color  = vec4(0,0,0,1);
	vec4 epos   = (gl_ModelViewMatrix * vertex);
	bump_scale *= clamp((2.5 - 0.1*length(epos.xyz)), 0.0, 1.0); // decrease scale with distance to reduce tiling artifacts on sand and snow
	
	if (enable_light0) { // sun
		float spec      = spec_scale*(0.2*weights.b + 0.25*weights4); // grass and snow
		float shininess = 20.0*weights.b + 40.0*weights4;
		color += add_light_comp(normal, epos, 0, diffuse_scale, ambient_scale, spec, shininess, bump_scale);
	}
	if (enable_light1) {color += add_light_comp(normal, epos, 1, diffuse_scale, ambient_scale, 0.0, 1.0, bump_scale);} // moon
	if (enable_light2) {color += add_light_comp(normal, epos, 2, 1.0, 1.0, 0.0, 1.0, bump_scale) * calc_light_atten(epos, 2);} // lightning
	gl_FragColor = apply_fog_scaled(vec4((texel0.rgb * texel1.rgb * color.rgb), color.a), vertex.z); // add fog
	//gl_FragColor = apply_fog(color); // untextured (white) for debugging
}
