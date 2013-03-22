uniform sampler2D tex0, tex1, tex2, tex3, tex4, tex5, tex6;
uniform float ts2, ts3, ts4, ts5, ts6; // texture scales
uniform float cs2, cs3, cs4, cs5, cs6; // color scales
uniform float water_plane_z, cloud_plane_z;
uniform float water_atten    = 1.0;
uniform float normal_z_scale = 1.0;
uniform float spec_scale     = 1.0;
uniform float cloud_alpha    = 1.0;
uniform vec3 cloud_offset    = vec3(0,0,0);
uniform sampler2D shadow_normal_tex, noise_tex;
varying vec3 vertex; // world space
varying vec4 epos;

// underwater attenuation code
void atten_color(inout vec4 color, in float dist) {
	color.rgb *= vec3(1,1,1) - min(vec3(0.98, 0.97, 0.95), vec3(1.5, 0.9, 0.5)*dist*2.0);
}

float integrate_water_dist(in vec3 targ_pos, in vec3 src_pos, in float water_z) {
	float t = clamp((water_z - targ_pos.z)/max(1.0E-6, abs(src_pos.z - targ_pos.z)), 0.0, 1.0);
	return t*distance(src_pos, targ_pos);
}

vec4 get_light_specular(in vec3 normal, in vec3 light_dir, in vec3 eye_pos, in float spec, in float shininess) {
	vec3 half_vect = normalize(light_dir - normalize(eye_pos)); // Eye + L = -eye_space_pos + L
	return vec4(spec, spec, spec, 1.0) * pow(max(dot(normal, half_vect), 0.0), shininess);
}

vec4 add_light_comp(in vec3 normal, in int i, in float shadow_weight, in float spec, in float shininess) {
	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(gl_LightSource[i].position.xyz);
		
	// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex.
	float NdotL = dot(normal, lightDir);
	vec4 light  = gl_ModelViewMatrixInverse * gl_LightSource[i].position; // world space

	if (apply_cloud_shadows /*&& vertex.z > water_plane_z*//*&& vertex.z < cloud_plane_z*/) {
		vec3 cpos = vertex + cloud_offset;
		float t = (cloud_plane_z - cpos.z)/(light.z - cpos.z); // sky intersection position along vertex->light vector
		normal *= 1.0 - cloud_alpha*gen_cloud_alpha(cpos.xy + t*(light.xy - cpos.xy));
	}
	
	// compute the ambient and diffuse lighting
	vec4 diffuse = gl_LightSource[i].diffuse;
	vec4 ambient = gl_LightSource[i].ambient;
	vec4 color   = ambient + shadow_weight*max(dot(normal, lightDir), 0.0)*diffuse;
	if (enable_light0) {color += get_light_specular(normal, lightDir, epos.xyz, shadow_weight*spec, shininess);}
	
	// apply underwater attenuation
	// Note: ok if vertex is above the water, dist will come out as 0
	vec4 eye   = gl_ModelViewMatrixInverse[3]; // world space
	float dist = integrate_water_dist(vertex, eye.xyz, water_plane_z) + integrate_water_dist(vertex, light.xyz, water_plane_z);
	atten_color(color, dist*water_atten);
	return color;
}

float calc_light0_caustics() { // use vertex and water_plane_z
	return 1.0;
}

void main()
{
	// sand, dirt, grass, rock, snow
	vec2 tc = gl_TexCoord[0].st;
	vec4 weights = texture2D(tex0, tc);
	float weights4 = clamp((1.0 - weights.r - weights.g - weights.b - weights.a), 0.0, 1.0);
	vec3 texel0  = cs2*weights.r*texture2D(tex2, ts2*tc).rgb +
	               cs3*weights.g*texture2D(tex3, ts3*tc).rgb +
				   cs4*weights.b*texture2D(tex4, ts4*tc).rgb +
				   cs5*weights.a*texture2D(tex5, ts5*tc).rgb +
				   cs6*weights4 *texture2D(tex6, ts6*tc).rgb;
	vec3 texel1  = texture2D(tex1, gl_TexCoord[1].st).rgb; // detail texture

	vec4 shadow_normal = texture2D(shadow_normal_tex, tc);
	vec3 normal = normalize(gl_NormalMatrix * ((2.0*shadow_normal.xyz - 1.0) * vec3(1.0, 1.0, normal_z_scale))); // eye space
	//normal += 0.05*weights4*vec3(texture2D(noise_tex, 571.0*tc).r-0.5, texture2D(noise_tex, 714.0*tc).r-0.5, texture2D(noise_tex, 863.0*tc).r-0.5);
	vec4 color  = gl_LightModel.ambient;
	//texel0 = mix(vec3(0.05, 0.25, 0.05), texel0, shadow_normal.w*shadow_normal.w);
	
	if (enable_light0) {
		float spec      = spec_scale*(0.2*weights.b + 0.25*weights4); // grass and snow
		float shininess = 20.0*weights.b + 40.0*weights4;
		color += add_light_comp(normal, 0, shadow_normal.w*calc_light0_caustics(), spec, shininess);
	}
	if (enable_light1) {color += add_light_comp(normal, 1, shadow_normal.w, 0.0, 1.0);}
	if (enable_light2) {color += add_light_comp(normal, 2, 1.0, 0.0, 1.0) * calc_light_atten(epos, 2);}
	gl_FragColor = apply_fog(vec4((texel0.rgb * texel1.rgb * color.rgb), color.a)); // add fog
}
