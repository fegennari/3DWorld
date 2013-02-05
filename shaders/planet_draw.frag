// Note: Light 0 is the sun (A+D+S point light), light 1 is universe ambient (constant A), light 2 is planet reflection (D point light)
uniform float atmosphere = 1.0;
uniform float water_val  = 0.0;
uniform float lava_val   = 0.0;
uniform vec3 light_scale = vec3(1,1,1);
uniform vec3 sun_pos, ss_pos, rscale;
uniform float sun_radius, ss_radius, ring_ri, ring_ro;
uniform mat4 world_space_mvm;
uniform sampler1D ring_tex;

#ifdef GAS_GIANT
uniform sampler1D tex0;
#else
uniform sampler2D tex0;
#endif

varying vec4 epos;
varying vec3 normal, world_space_pos, vertex;

#if 0

uniform float planet_radius, hmap_scale, snow_thresh, temperature, wr_scale;
uniform vec3 color_a, color_b;

vec3 eval_color(in float val, in float zval) {
	const float FREEZE_TEMP = 12.0;
	const float BOIL_TEMP   = 30.0;
	const vec3 white = vec3(1.0, 1.0, 1.0);
	const vec3 gray  = vec3(0.4, 0.4, 0.4);
	vec3 wic = ((temperature < FREEZE_TEMP) ? vec3(0.5, 0.5, 0.9) : vec3(0.2, 0.2, 0.8));
	float coldness = pow(abs(zval - 0.5)*2.0, 4.0); // zval=0.5 => equator, zval=0.0 => north pole, zval=1.0 => south pole
	vec3 color = vec3(0,0,0);

	if (val < water_val) { // underwater
		color = wic;
		if (coldness > 0.8) {color = mix(color, white, clamp((5.0*(coldness - 0.8f) + 2.0*(val - water_val)), 0.0, 1.0));} // ice
		return color;
	}
	const float water_adj = 0.07;
	float val_ws = ((water_val > 0.0) ? wr_scale*(val - water_val) : val); // rescale for water

	if (water_val > 0.2 && atmosphere > 0.1) { // Earthlike planet
		if (val_ws < 0.1) { // low ground
			color = color_b;
		}
		else if (val_ws < 0.4) {
			color = mix(color_b, color_a, 3.3333*(val_ws - 0.1));
		}
		else if (val_ws < 0.45) { // medium ground
			color = color_a;
		}
		else if (val_ws < 0.6) {
			color = mix(color_a, gray, 6.6667*(val_ws - 0.45));
		}
		else { // high ground
			color = gray;
		}
	}
	else { // alien-like planet
		color = mix(color_b, color_a, val_ws);
	}
	if (lava_val > 0.0) { // hot lava planet
		const float lava_adj = 0.07;
		const vec3 lavac = vec3(1.0, 0.0, 0.0); // red

		if (val < lava_val) {
			color = lavac; // move up?
		}
		else if (val < lava_val + lava_adj) { // close to lava line
			color = mix(lavac, color, (val - lava_val)/lava_adj);
		}
	}
	else if (temperature < BOIL_TEMP) { // handle water/ice/snow
		if (val < water_val + water_adj) { // close to water line (can have a little water even if water == 0)
			color = mix(wic, color, (val - water_val)/water_adj);
			if (coldness > 0.9) {color = mix(color, white, 10.0*(coldness - 0.9));} // ice
			return color;
		}
		float st = (1.0 - coldness*coldness)*snow_thresh;

		if (val > (st + 1.0E-6)) { // blend in some snow
			color = mix(color, white, (val - st)/(1.0 - st));
			return color;
		}
	}
	return color;
}

vec3 eval_color_cur_pos() {

	float cutoff = max(water_val, lava_val);
	float val    = cutoff + (1.0 - cutoff)*((length(vertex) - planet_radius)/(hmap_scale*planet_radius) + 0.5);
	val += 0.05*(gen_cloud_alpha_static(5.0*vertex) - 0.5);
	float zval   = 0.5*vertex.z/planet_radius + 0.5;
	return eval_color(val, zval);
}

#endif


void main()
{
#ifdef GAS_GIANT
	float tc     = gl_TexCoord[0].t + 0.04*(gen_cloud_alpha_static_non_norm(1.5*vertex) - 0.5);
	vec4 texel   = texture1D(tex0, tc);
#else
	vec4 texel   = texture2D(tex0, gl_TexCoord[0].st);
	//vec4 texel   = vec4(eval_color_cur_pos(), 1.0);
#endif
	vec3 norm    = normalize(normal); // renormalize
	float atten0 = light_scale[0] * calc_light_atten(epos, 0);
	float atten2 = light_scale[2] * calc_light_atten(epos, 2);

	if (sun_radius > 0.0) {
		if (ss_radius > 0.0) {
			atten0 *= calc_sphere_shadow_atten(world_space_pos, sun_pos, sun_radius, ss_pos, ss_radius);
		}
		if (has_rings) { // calculate shadows due to rings
			vec3 sun_local = (gl_ModelViewMatrixInverse * (world_space_mvm * vec4(sun_pos, 1.0))).xyz;
			vec3 line_dir  = sun_local - vertex;
			float dist     = -vertex.z/line_dir.z; // Note: ring normal is always in z

			if (dist > 0.0) {
				vec3 ring_ipt = dist*line_dir + vertex;
				float rval    = (length(ring_ipt/rscale) - ring_ri)/(ring_ro - ring_ri);
				
				if (rval > 0.0 && rval < 1.0) {
					float dscale = length(vertex)/length(ring_ipt); // fake penumbra
					atten0      *= 1.0 - dscale*texture1D(ring_tex, rval).a;
				}
			}
		}
	}
	vec3 ldir0     = normalize(gl_LightSource[0].position.xyz - epos.xyz);
	vec3 ldir2     = normalize(gl_LightSource[2].position.xyz - epos.xyz);
	vec3 ldir20    = normalize(gl_LightSource[2].position.xyz - gl_LightSource[0].position.xyz);
	vec3 epos_norm = normalize(epos.xyz);
	vec3 half_vect = normalize(ldir0 - epos_norm); // Eye + L = -eye_space_pos + L
	vec3 ambient   = (gl_LightSource[0].ambient.rgb * atten0) + (gl_LightSource[1].ambient.rgb * light_scale[1]);
	vec3 diffuse   = (gl_LightSource[0].diffuse.rgb * max(dot(norm, ldir0), 0.0) * atten0) +
	                 (gl_LightSource[2].diffuse.rgb * max(dot(norm, ldir2), 0.0) * atten2 * max(dot(ldir2, ldir20), 0.0));
	vec3 specular  = ((water_val > 0.0) ? 1.0 : 0.0) * gl_FrontLightProduct[0].specular.rgb * pow(max(dot(norm, half_vect), 0.0), gl_FrontMaterial.shininess) * pow(texel.b, 4.0) * atten0;
	vec3 color     = (texel.rgb * (ambient + diffuse)) + specular;

#ifndef GAS_GIANT
	if (atmosphere > 0.0) {
		float cloud_val = atmosphere*gen_cloud_alpha(vertex);
		if (cloud_val > 0.0) {color = cloud_val*(ambient + diffuse) + (1.0 - cloud_val)*color;} // no clouds over high mountains?
	}
	if (lava_val > 0.0) {
		float heat = max(0.0, (texel.r - texel.g - texel.b));

		if (heat > 0.0) { // lava patch
			heat *= gen_cloud_alpha(2.0*vertex);
			//color = mix(color, vec3(1.0, 0.25*heat, 0.0), heat); // add lava
			color += heat*vec3(1.0, 0.25*heat, 0.0); // add lava
		}
	}
#endif
	gl_FragColor = gl_Color * vec4((color + gl_FrontMaterial.emission.rgb), 1.0);
}
