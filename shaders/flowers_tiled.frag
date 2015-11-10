uniform float min_alpha  = 0.0;
uniform float dist_const = 10.0;
uniform float dist_slope = 0.5;
uniform float dx_inv, dy_inv;
uniform sampler2D tex0, shadow_tex;

in vec4 vertex, epos;
in vec3 eye_norm;
in vec2 tc;

void main()
{
	vec4 texel = texture(tex0, tc);
	if (texel.a <= min_alpha) discard;

	float ds_val = dist_slope*(length(epos.xyz) - dist_const);
	if (ds_val >= 1.0) discard; // too far away
	texel.a *= (1.0 - max(ds_val, 0.0)); // decrease transparency when far

	vec2 tc2     = vec2(vertex.x*dx_inv, vertex.y*dy_inv); // same as (x2 - x1 - 1.0*DX_VAL)
	vec3 shadow  = texture(shadow_tex, tc2).rgb; // {mesh_shadow, tree_shadow, ambient_occlusion}
	float diffuse_scale = min(shadow.r, shadow.g); // min of mesh and tree shadow
	vec3 color   = do_shadowed_lighting(vertex, epos, eye_norm, gl_Color, shadow.b, diffuse_scale);
	fg_FragColor = apply_fog_epos(texel*vec4(color, gl_Color.a*texel.a), epos);
}
