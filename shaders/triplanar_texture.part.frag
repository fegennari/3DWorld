uniform float tex_scale = 1.0;

vec4 lookup_triplanar_texture(in vec3 pos, in vec3 normal, in sampler2D tex_x, in sampler2D tex_y, in sampler2D tex_z) {
	// http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
	// Determine the blend weights for the 3 planar projections.  
	vec3 blend_weights = abs(normal);   // Tighten up the blending zone:  
	blend_weights = (blend_weights - 0.2) * 7;  
	blend_weights = max(blend_weights, 0);      // Force weights to sum to 1.0 (very important!)  
	blend_weights /= (blend_weights.x + blend_weights.y + blend_weights.z ).xxx;

	// Now determine a color value for each of the 3 projections, blend them, and store blended results
	// Compute the UV coords for each of the 3 planar projections.  
	vec2 coord_x = pos.yz * tex_scale;  
	vec2 coord_y = pos.zx * tex_scale;  
	vec2 coord_z = pos.xy * tex_scale;  

	// Sample color maps for each projection, at those UV coords.  
	vec4 col_x = texture2D(tex_x, coord_x);
	vec4 col_y = texture2D(tex_y, coord_y);
	vec4 col_z = texture2D(tex_z, coord_z);

	// Finally, blend the results of the 3 planar projections.  
	return (col_x*blend_weights.xxxx + col_y*blend_weights.yyyy + col_z*blend_weights.zzzz);
}
