uniform float dist_const  = 10.0;
uniform float dist_slope  = 0.5;
uniform float cloud_alpha = 1.0;
uniform float x1, y1, dx_inv, dy_inv;
uniform sampler2D height_tex, shadow_normal_tex, weight_tex, noise_tex;
uniform float cloud_plane_z;
uniform vec3 cloud_offset = vec3(0.0);
uniform vec2 xlate = vec2(0.0);

attribute vec2 local_translate;

varying vec2 tc;

float calc_light_scale(in vec3 vertex, in vec4 light_pos) {
	if (!apply_cloud_shadows) {return 1.0;}
	vec4 light = fg_ModelViewMatrixInverse * light_pos; // world space
	vec3 cpos  = vertex + cloud_offset;
	float t    = (cloud_plane_z - cpos.z)/(light.z - cpos.z); // sky intersection position along vertex->light vector
	return 1.0 - cloud_alpha*gen_cloud_alpha(cpos.xy + t*(light.xy - cpos.xy));
}

void main()
{
	tc          = get_grass_tc();
	vec4 vertex = fg_Vertex;
	vertex.xy  += local_translate;
	float z_val = texture2D(height_tex, vec2((vertex.x - x1)*dx_inv, (vertex.y - y1)*dy_inv)).r;
	float ascale= 1.0;
#ifdef DEC_HEIGHT_WHEN_FAR
	float dist  = length((fg_ModelViewMatrix * (vertex + vec4(xlate, z_val, 0))).xyz);
	float zscale= 1.0 - clamp(dist_slope*(dist - dist_const), 0.0, 1.0);
	vertex.z   *= zscale; // decrease height far away from camera
	ascale      = min(1.0, 10.0*zscale);
#endif
	vertex.z   += z_val;
	vec2 tc2    = vec2(vertex.x*dx_inv, vertex.y*dy_inv); // same as (x2 - x1 - 1.0*DX_VAL)
	if (enable_grass_wind) {vertex.xyz += get_grass_wind_delta(vertex.xyz, tc.s);}

	vec4 epos   = fg_ModelViewMatrix  * (vertex + vec4(xlate, 0, 0));
	gl_Position = fg_ProjectionMatrix * epos;
	gl_FogFragCoord = length(epos.xyz);
	float grass_weight = texture2D(weight_tex, tc2).b;
	//grass_weight = ((grass_weight < 0.2) ? 0.0 : grass_weight);
	float noise_weight = texture2D(noise_tex, 10.0*vec2(fg_Color.r, fg_Color.g)).r; // "hash" the color
	
	// calculate lighting
	vec4 shadow_normal  = texture2D(shadow_normal_tex, tc2);
	float diffuse_scale = shadow_normal.w;
	float ambient_scale = 1.5*shadow_normal.z * (1.5 - 0.75*tc.s); // decreased ambient at base, increased ambient at tip
	vec2 nxy    = (2.0*shadow_normal.xy - 1.0);
	vec3 normal = vec3(nxy, (1.0 - sqrt(nxy.x*nxy.x + nxy.y*nxy.y))); // calculate n.z from n.x and n.y (we know it's always positive)
	normal      = normalize(fg_NormalMatrix * normal); // eye space
	vec3 color  = vec3(0.0);
	//if (grass_weight < noise_weight) {
	if (enable_light0) {color += add_light_comp_pos_scaled_smap_light0(normal, epos, diffuse_scale*calc_light_scale(vertex.xyz, fg_LightSource[0].position), ambient_scale).rgb;}
	if (enable_light1) {color += add_light_comp_pos_scaled_smap_light1(normal, epos, diffuse_scale*calc_light_scale(vertex.xyz, fg_LightSource[1].position), ambient_scale).rgb;}
	if (enable_light2) {color += add_pt_light_comp(normal, epos, 2).rgb;}
	//}
	float alpha = fg_Color.a;
	alpha *= ascale*((grass_weight < noise_weight) ? 0.0 : 1.0); // skip some grass blades by making them transparent
	gl_FrontColor = vec4(color, alpha);
} 
