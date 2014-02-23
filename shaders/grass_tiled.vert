uniform float dist_const  = 10.0;
uniform float dist_slope  = 0.5;
uniform float cloud_alpha = 1.0;
uniform float x1, y1, dx_inv, dy_inv;
uniform sampler2D height_tex, shadow_normal_tex, weight_tex, noise_tex;
uniform float cloud_plane_z;
uniform vec3 cloud_offset = vec3(0.0);

attribute vec2 local_translate;

varying vec2 tc;

vec4 add_light_comp(in vec3 vertex, in vec3 normal, in vec4 epos, in int i, in float ds_scale, in float a_scale) {
	if (apply_cloud_shadows) {
		vec4 light = gl_ModelViewMatrixInverse * gl_LightSource[i].position; // world space
		vec3 cpos  = vertex + cloud_offset;
		float t    = (cloud_plane_z - cpos.z)/(light.z - cpos.z); // sky intersection position along vertex->light vector
		ds_scale  *= 1.0 - cloud_alpha*gen_cloud_alpha(cpos.xy + t*(light.xy - cpos.xy));
	}
	return add_light_comp_pos_scaled(normal, epos, i, ds_scale, a_scale);
}

void main()
{
	tc          = get_grass_tc();
	vec4 vertex = gl_Vertex;
	vertex.xy  += local_translate;
	float z_val = texture2D(height_tex, vec2((vertex.x - x1)*dx_inv, (vertex.y - y1)*dy_inv)).r;
	float ascale= 1.0;
#ifdef DEC_HEIGHT_WHEN_FAR
	float dist  = length((gl_ModelViewMatrix * (vertex + vec4(0, 0, z_val, 0))).xyz);
	float zscale= 1.0 - clamp(dist_slope*(dist - dist_const), 0.0, 1.0);
	vertex.z   *= zscale; // decrease height far away from camera
	ascale      = min(1.0, 10.0*zscale);
#endif
	vertex.z   += z_val;
	vec2 tc2    = vec2(vertex.x*dx_inv, vertex.y*dy_inv); // same as (x2 - x1 - 1.0*DX_VAL)
	if (enable_grass_wind) {vertex.xyz += get_grass_wind_delta(vertex.xyz, tc.s);}

	vec4 epos   = gl_ModelViewMatrix  * vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	gl_FogFragCoord = length(epos.xyz);
	float grass_weight = texture2D(weight_tex, tc2).b;
	//grass_weight = ((grass_weight < 0.2) ? 0.0 : grass_weight);
	float noise_weight = texture2D(noise_tex, 10.0*vec2(gl_Color.r, gl_Color.g)).r; // "hash" the color
	
	// calculate lighting
	vec4 shadow_normal  = texture2D(shadow_normal_tex, tc2);
	float diffuse_scale = shadow_normal.w;
	float ambient_scale = 1.5*shadow_normal.z * (1.5 - 0.75*tc.s); // decreased ambient at base, increased ambient at tip
	vec2 nxy    = (2.0*shadow_normal.xy - 1.0);
	vec3 normal = vec3(nxy, (1.0 - sqrt(nxy.x*nxy.x + nxy.y*nxy.y))); // calculate n.z from n.x and n.y (we know it's always positive)
	normal      = normalize(gl_NormalMatrix * normal); // eye space
	vec4 color  = gl_FrontMaterial.emission;
	//if (grass_weight < noise_weight) {
	if (enable_light0) {color += add_light_comp(vertex.xyz, normal, epos, 0, diffuse_scale, ambient_scale);}
	if (enable_light1) {color += add_light_comp(vertex.xyz, normal, epos, 1, diffuse_scale, ambient_scale);}
	if (enable_light2) {color += add_pt_light_comp(normal, epos, 2);}
	//}
	color.a *= ascale*((grass_weight < noise_weight) ? 0.0 : 1.0); // skip some grass blades by making them transparent
	gl_FrontColor = color;
} 
