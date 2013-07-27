// Note: Light 0 is the sun (A+D+S point light), light 1 is universe ambient (constant A), light 2 is planet reflection (D point light)
uniform float atmosphere = 1.0; // technically not needed for gas giants since assumed to be 1.0
uniform vec3 cloud_freq  = vec3(1.0);
uniform vec3 light_scale = vec3(1.0);
uniform vec3 sun_pos, ss_pos, rscale;
uniform float sun_radius, ss_radius, ring_ri, ring_ro;
uniform mat4 world_space_mvm;
uniform sampler1D ring_tex;

#ifdef GAS_GIANT
uniform sampler1D tex0;
#else
uniform float water_val  = 0.0;
uniform float lava_val   = 0.0;
uniform float crater_val = 0.0;
uniform sampler2D tex0;
#endif

varying vec3 normal, world_space_pos, vertex;


void main()
{
	vec4 epos = gl_ModelViewMatrix * vec4(vertex, 1.0);
	if (dot(normal, epos.xyz) > 0.0) discard; // back facing

#ifdef GAS_GIANT
	float tc = gl_TexCoord[0].t + 0.04*(gen_cloud_alpha_static_non_norm(5.0*vertex) - 0.5);
	float v0 = 1.0; // using a variable here is slow
	vec3 dir = normalize(vertex); // world space normal

	for (int i = 0; i < 50; ++i) { // Note: inefficient, but fast enough for a single gas giant render
		vec3 center = vec3(1.0, 1.0, 0.5)*rand_vec3(v0);
#ifdef ANIMATE_STORMS // slow but neat
		float angle = 50.0*time*rand_pm1(v0+2.5);
		float st    = sin(angle);
		float ct    = cos(angle);
		center.xy   = vec2((center.x*ct - center.y*st), (center.y*ct + center.x*st)); // rotation
#endif
		float dist  = (0.25 + 0.75*rand_01(v0+3.0))*length(vec3(1.0, 1.0, 2.0)*(dir - normalize(center)));
		tc         += 0.5*max(0.0, (0.1 - dist))*sin(0.1/max(dist, 0.01));
		v0         += 4.0;
	}
	vec4 texel   = texture1D(tex0, tc);
#else
	vec4 texel   = texture2D(tex0, gl_TexCoord[0].st);
#endif
	float atten0 = light_scale[0] * calc_light_atten(epos, 0);
	float atten2 = light_scale[2] * calc_light_atten(epos, 2);
	float sscale = atten0;

	if (sun_radius > 0.0) {
		if (ss_radius > 0.0) {
			sscale *= calc_sphere_shadow_atten(world_space_pos, sun_pos, sun_radius, ss_pos, ss_radius);
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
					sscale      *= 1.0 - dscale*texture1D(ring_tex, rval).a;
				}
			}
		}
	}
	vec3 norm      = normalize(normal); // renormalize
	vec3 ldir0     = normalize(gl_LightSource[0].position.xyz - epos.xyz);
	vec3 ldir2     = normalize(gl_LightSource[2].position.xyz - epos.xyz);
	vec3 ldir20    = normalize(gl_LightSource[2].position.xyz - gl_LightSource[0].position.xyz);
	float lscale0  = (dot(norm, ldir0) > 0.0) ? 1.0 : 0.0;
	float lscale2  = (dot(norm, ldir2) > 0.0) ? 1.0 : 0.0;

#ifdef HAS_CRATERS
	// facing the sun or planet (reflected light), and not over water (blue)
	if ((lscale0 > 0.0 || lscale2 > 0.0) && (texel.b - texel.r - texel.g) < 0.0) {
		adjust_normal_for_craters(norm, vertex); // add craters by modifying the normal
	}
#endif

	vec3 epos_norm = normalize(epos.xyz);
	vec3 half_vect = normalize(ldir0 - epos_norm); // Eye + L = -eye_space_pos + L
	vec3 ambient   = (gl_LightSource[0].ambient.rgb * atten0) + (gl_LightSource[1].ambient.rgb * light_scale[1]);
	vec3 diffuse   = (gl_LightSource[0].diffuse.rgb * max(dot(norm, ldir0), 0.0) * lscale0 * sscale) +
	                 (gl_LightSource[2].diffuse.rgb * max(dot(norm, ldir2), 0.0) * lscale2 * atten2 * max(dot(ldir2, ldir20), 0.0));
	vec3 color     = (texel.rgb * (ambient + diffuse));

#ifndef GAS_GIANT
	float specval  = pow(max(dot(norm, half_vect), 0.0), gl_FrontMaterial.shininess);
	color         += ((water_val > 0.0) ? 1.0 : 0.0) * gl_FrontLightProduct[0].specular.rgb * specval * pow(texel.b, 4.0) * sscale;

	if (lava_val > 0.0) {
		float heat = max(0.0, (texel.r - texel.g - texel.b));

		if (heat > 0.0) { // lava patch
			heat *= gen_cloud_alpha(2.0*vertex);
			//color = mix(color, vec3(1.0, 0.25*heat, 0.0), heat); // add lava
			color += heat*vec3(1.0, 0.25*heat, 0.0); // add lava
		}
	}
#endif
	if (atmosphere > 0.0) {
		float cloud_val = atmosphere*gen_cloud_alpha(cloud_freq*vertex);
		if (cloud_val > 0.0) {color = cloud_val*(ambient + diffuse) + (1.0 - cloud_val)*color;} // no clouds over high mountains?
	}
	gl_FragColor = gl_Color * vec4((color + gl_FrontMaterial.emission.rgb), 1.0);
}
