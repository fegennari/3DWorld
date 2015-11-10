uniform float dist_const = 10.0;
uniform float dist_slope = 0.5;
uniform float x1, y1, dx_inv, dy_inv;
uniform sampler2D height_tex, normal_tex, shadow_tex, weight_tex, noise_tex;
uniform vec2 xlate = vec2(0.0);

in vec2 local_translate;

out vec2 tc;

void main()
{
	tc          = get_grass_tc();
	vec4 vertex = fg_Vertex;
	vertex.xy  += local_translate;
	float z_val = texture(height_tex, vec2((vertex.x - x1)*dx_inv, (vertex.y - y1)*dy_inv)).r;
	float ascale= 1.0;
#ifdef DEC_HEIGHT_WHEN_FAR
	float dist  = length((fg_ModelViewMatrix * (vertex + vec4(xlate, z_val, 0))).xyz);
	float ds_val= clamp(dist_slope*(dist - dist_const), 0.0, 1.0);
	vertex.z   -= 0.9*ds_val*height; // move below mesh far away from camera
	ascale      = min(1.0, 10.0*(1.0 - ds_val)); // decrease alpha very far away from camera
#endif
	vertex.z   += z_val;
	vec2 tc2    = vec2(vertex.x*dx_inv, vertex.y*dy_inv); // same as (x2 - x1 - 1.0*DX_VAL)
	if (enable_grass_wind) {vertex.xyz += get_grass_wind_delta(vertex.xyz, tc.s);}

	vec4 epos   = fg_ModelViewMatrix  * (vertex + vec4(xlate, 0, 0));
	gl_Position = fg_ProjectionMatrix * epos;
	gl_FogFragCoord = length(epos.xyz);
	vec4 weights = texture(weight_tex, tc2);
	float grass_weight = weights.b; // grass weight in weights {sand, dirt, grass, rock, [snow]}
	//grass_weight = ((grass_weight < 0.2) ? 0.0 : grass_weight);
	float noise_weight = texture(noise_tex, 10.0*vec2(fg_Color.r, fg_Color.g)).r; // "hash" the color
	
	// calculate lighting
	vec3 shadow  = texture(shadow_tex, tc2).rgb; // {mesh_shadow, tree_shadow, ambient_occlusion}
	float ambient_scale = 1.5*shadow.b * (1.5 - 0.75*tc.s); // decreased ambient at base, increased ambient at tip
	vec3 eye_norm = normalize(fg_NormalMatrix * (2.0*texture(normal_tex, tc2).xyz - vec3(1.0))); // eye space
	vec4 ad_color = mix(gl_Color, vec4(1.0, 0.7, 0.4, 1.0), weights.r); // mix in yellow-brown grass color to match sand
	float diffuse_scale = min(shadow.r, shadow.g); // min of mesh and tree shadow
	vec3 color    = do_shadowed_lighting(vertex, epos, eye_norm, ad_color, ambient_scale, diffuse_scale);
	float alpha   = fg_Color.a * ascale * ((grass_weight < noise_weight) ? 0.0 : 1.0); // skip some grass blades by making them transparent
	fg_Color_vf   = vec4(color, alpha);
} 
