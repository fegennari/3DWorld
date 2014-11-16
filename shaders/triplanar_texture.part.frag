uniform float tex_scale = 1.0;

vec3 get_blend_weights(in vec3 normal) {
	vec3 blend_weights = abs(normal); // Tighten up the blending zone:
	blend_weights = (blend_weights - 0.2) * 7;
	blend_weights = max(blend_weights, 0); // Force weights to sum to 1.0 (very important!)
	blend_weights /= (blend_weights.x + blend_weights.y + blend_weights.z);
	return blend_weights;
}

vec4 lookup_triplanar_texture(in vec3 pos, in vec3 normal, in sampler2D tex_x, in sampler2D tex_y, in sampler2D tex_z) {
	// http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
	// Determine the blend weights for the 3 planar projections.
	vec3 blend_weights = get_blend_weights(normal);

	// Now determine a color value for each of the 3 projections, blend them, and store blended results
	// Compute the UV coords for each of the 3 planar projections.
	vec2 coord_x = pos.yz * tex_scale;
	vec2 coord_y = pos.zx * tex_scale;
	vec2 coord_z = pos.xy * tex_scale;

	// Sample color maps for each projection, at those UV coords.
	vec4 col_x = texture(tex_x, coord_x);
	vec4 col_y = texture(tex_y, coord_y);
	vec4 col_z = texture(tex_z, coord_z);

	// Finally, blend the results of the 3 planar projections.
	return (col_x*blend_weights.xxxx + col_y*blend_weights.yyyy + col_z*blend_weights.zzzz);
}

vec4 lookup_triplanar_texture_2sz(in vec3 pos, in vec3 normal, in sampler2D tex_x, in sampler2D tex_y, in sampler2D tex_z1, in sampler2D tex_z2) { 
	vec3 blend_weights = get_blend_weights(normal); 
	vec2 coord_x = pos.yz * tex_scale;
	vec2 coord_y = pos.zx * tex_scale;
	vec2 coord_z = pos.xy * tex_scale;
	vec4 col_x   = texture(tex_x, coord_x);
	vec4 col_y   = texture(tex_y, coord_y);
	vec4 col_z1  = texture(tex_z1, coord_z);
	vec4 col_z2  = texture(tex_z2, coord_z);
	vec4 col_z   = mix(col_z1, col_z2, clamp(1000*normal.z, 0.0, 1.0));
	return (col_x*blend_weights.xxxx + col_y*blend_weights.yyyy + col_z*blend_weights.zzzz);
}

vec3 lookup_triplanar_texture_bump(in vec3 pos, in vec3 normal, in sampler2D bump_tex) {
	vec3 blend_weights = get_blend_weights(normal);
	vec2 bump_tx = texture(bump_tex, pos.yz).xy - 0.5;
	vec2 bump_ty = texture(bump_tex, pos.zx).xy - 0.5;
	vec2 bump_tz = texture(bump_tex, pos.xy).xy - 0.5;
	vec3 bump_x  = vec3(0.0, bump_tx.x, bump_tx.y);
	vec3 bump_y  = vec3(bump_ty.y, 0.0, bump_ty.x);
	vec3 bump_z  = vec3(bump_tz.x, bump_tz.y, 0.0);
	return normalize(bump_x*blend_weights.xxx + bump_y*blend_weights.yyy + bump_z*blend_weights.zzz);
}
