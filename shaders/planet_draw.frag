// Note: Light 0 is the sun (A+D+S point light), light 1 is universe ambient (constant A), light 2 is planet reflection (D point light)
uniform float atmosphere = 1.0; // technically not needed for gas giants since assumed to be 1.0
uniform vec3 light_scale = vec3(1,1,1);
uniform vec3 sun_pos, ss_pos, rscale;
uniform float sun_radius, ss_radius, ring_ri, ring_ro;
uniform mat4 world_space_mvm;
uniform sampler1D ring_tex;

#ifdef GAS_GIANT
uniform sampler1D tex0;
#else
uniform float water_val = 0.0;
uniform float lava_val  = 0.0;
uniform sampler2D tex0;
#endif

varying vec4 epos;
varying vec3 normal, world_normal, world_space_pos, vertex;

float rand_01 (float val) {return fract(sin(12.9898 * val) * 43758.5453);}
float rand_pm1(float val) {return 2.0*(rand_01(val) - 0.5);}


void main()
{
#ifdef GAS_GIANT
	float tc = gl_TexCoord[0].t + 0.04*(gen_cloud_alpha_static_non_norm(1.5*vertex) - 0.5);
	float v0 = 1.0; // using a variable here is slow

	for (int i = 0; i < 50; ++i) { // Note: inefficient, but fast enough for a single gas giant render
		float v1   = rand_pm1(v0+1.0);
		float v2   = rand_pm1(v0+2.0);
		float v3   = rand_pm1(v0+3.0);
		float v4   = rand_01 (v0+4.0);
		v0        += 5.0;
		float dist = (0.25 + 0.75*v4)*length(vec3(1.0, 1.0, 2.0)*(world_normal - normalize(vec3(v1, v2, 0.5*v3))));
		tc        += 0.5*max(0.0, (0.1 - dist))*sin(0.1/dist);
	}
	vec4 texel   = texture1D(tex0, tc);
#else
	vec4 texel   = texture2D(tex0, gl_TexCoord[0].st);
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
	vec3 color     = (texel.rgb * (ambient + diffuse));

#ifndef GAS_GIANT
	float specval  = pow(max(dot(norm, half_vect), 0.0), gl_FrontMaterial.shininess);
	color         += ((water_val > 0.0) ? 1.0 : 0.0) * gl_FrontLightProduct[0].specular.rgb * specval * pow(texel.b, 4.0) * atten0;

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
