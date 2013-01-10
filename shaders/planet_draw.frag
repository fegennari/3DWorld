// Note: Light 0 is the sun (A+D+S point light), light 1 is universe ambient (constant A), light 2 is planet reflection (D point light)
uniform float atmosphere = 1.0;
uniform float spec_scale = 1.0;
uniform vec3 light_scale = vec3(1,1,1);
uniform vec3 sun_pos, ss_pos, rscale;
uniform float sun_radius, ss_radius, ring_ri, ring_ro;
uniform mat4 world_space_mvm;
uniform sampler2D tex0;
uniform sampler1D ring_tex;
varying vec4 epos;
varying vec3 normal, world_space_pos, vertex;


void main()
{
	vec4 texel   = texture2D(tex0, gl_TexCoord[0].st);
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
	vec3 specular  = spec_scale * gl_FrontLightProduct[0].specular.rgb * pow(max(dot(norm, half_vect), 0.0), gl_FrontMaterial.shininess) * pow(texel.b, 4.0) * atten0;
	vec3 color     = (texel.rgb * (ambient + diffuse)) + specular;
	float cloud_val= atmosphere*gen_cloud_alpha(vertex);
	if (cloud_val > 0.0) {color = cloud_val*(ambient + diffuse) + (1.0 - cloud_val)*color;} // no clouds over high mountains?
	gl_FragColor   = apply_fog(gl_Color * vec4((color + gl_FrontMaterial.emission.rgb), 1.0));
}
