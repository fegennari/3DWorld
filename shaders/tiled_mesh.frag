uniform sampler2D weights_tex, detail_tex, tex2, tex3, tex4, tex5, tex6, shadow_normal_tex, noise_tex, caustic_tex;
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
varying vec4 vertex; // world space

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

vec4 add_light_comp(in vec3 normal, in vec4 epos, in int i, in float shadow_weight, in float spec, in float shininess) {
	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(gl_LightSource[i].position.xyz);
		
	// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex
	float NdotL = dot(normal, lightDir);
	vec4 light  = gl_ModelViewMatrixInverse * gl_LightSource[i].position; // world space

	if (apply_cloud_shadows /*&& vertex.z > water_plane_z*//*&& vertex.z < cloud_plane_z*/) {
		vec3 cpos = vertex.xyz + cloud_offset;
		float t = (cloud_plane_z - cpos.z)/(light.z - cpos.z); // sky intersection position along vertex->light vector
		normal *= 1.0 - cloud_alpha*gen_cloud_alpha(cpos.xy + t*(light.xy - cpos.xy));
	}
	
	// compute the ambient and diffuse lighting
	vec4 diffuse = gl_LightSource[i].diffuse;
	vec4 ambient = gl_LightSource[i].ambient;
	vec4 color   = ambient + shadow_weight*max(dot(normal, lightDir), 0.0)*diffuse;
	if (enable_light0) {color += get_light_specular(normal, lightDir, epos.xyz, shadow_weight*spec, shininess);}
	
#ifdef HAS_WATER
	if (vertex.z < water_plane_z) { // underwater
#ifdef WATER_CAUSTICS
		if (i == 0) { // only for light0 (sun)
			// apply underwater caustics texture
			float cweight = shadow_weight*wave_amplitude*caustics_weight*min(8.0*(water_plane_z - vertex.z), 0.5);
			float ntime   = 2.0*abs(fract(0.005*wave_time) - 0.5);
			vec3  cval    = 4.0*mix(texture2D(caustic_tex, gl_TexCoord[2].st).rgb, texture2D(caustic_tex, (gl_TexCoord[2].st + vec2(0.3, 0.6))).rgb, ntime);
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
	// sand, dirt, grass, rock, snow
	vec2 tc = gl_TexCoord[0].st;
	vec2 diff_tc = tc; // separate tc for diffuse texture, in case we want to sometimes mirror it to make tiling less periodic (though seems difficult and unnecessary)
	vec4 weights = texture2D(weights_tex, tc);
	float weights4 = clamp((1.0 - weights.r - weights.g - weights.b - weights.a), 0.0, 1.0);
	vec3 texel0  = cs2*weights.r*texture2D(tex2, ts2*diff_tc).rgb +
	               cs3*weights.g*texture2D(tex3, ts3*diff_tc).rgb +
				   cs4*weights.b*texture2D(tex4, ts4*diff_tc).rgb +
				   cs5*weights.a*texture2D(tex5, ts5*diff_tc).rgb +
				   cs6*weights4 *texture2D(tex6, ts6*diff_tc).rgb;
	vec3 texel1  = texture2D(detail_tex, gl_TexCoord[1].st).rgb; // detail texture

	vec4 shadow_normal = texture2D(shadow_normal_tex, tc);
	vec3 normal = normalize(gl_NormalMatrix * ((2.0*shadow_normal.xyz - 1.0) * vec3(1.0, 1.0, normal_z_scale))); // eye space
	//normal += 0.05*weights4*vec3(texture2D(noise_tex, 571.0*tc).r-0.5, texture2D(noise_tex, 714.0*tc).r-0.5, texture2D(noise_tex, 863.0*tc).r-0.5);
	vec4 color = gl_LightModel.ambient;
	//texel0 = mix(vec3(0.05, 0.25, 0.05), texel0, shadow_normal.w*shadow_normal.w);
	vec4 epos = (gl_ModelViewMatrix * vertex);
	
	if (enable_light0) {
		float spec      = spec_scale*(0.2*weights.b + 0.25*weights4); // grass and snow
		float shininess = 20.0*weights.b + 40.0*weights4;
		color += add_light_comp(normal, epos, 0, shadow_normal.w, spec, shininess);
	}
	if (enable_light1) {color += add_light_comp(normal, epos, 1, shadow_normal.w, 0.0, 1.0);}
	if (enable_light2) {color += add_light_comp(normal, epos, 2, 1.0, 0.0, 1.0) * calc_light_atten(epos, 2);}
	gl_FragColor = apply_fog(vec4((texel0.rgb * texel1.rgb * color.rgb), color.a)); // add fog
}
