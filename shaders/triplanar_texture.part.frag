uniform float tex_scale = 1.0;
uniform sampler2D tex0;

vec4 lookup_triplanar_texture(in vec3 pos, in vec3 normal) {
	// http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
	// Determine the blend weights for the 3 planar projections.  
	vec3 blend_weights = abs(normal);   // Tighten up the blending zone:  
	blend_weights = (blend_weights - 0.2) * 7;  
	blend_weights = max(blend_weights, 0);      // Force weights to sum to 1.0 (very important!)  
	blend_weights /= (blend_weights.x + blend_weights.y + blend_weights.z ).xxx;

	// Now determine a color value for each of the 3 projections, blend them, and store blended results
	// Compute the UV coords for each of the 3 planar projections.  
	vec2 coord1 = pos.yz * tex_scale;  
	vec2 coord2 = pos.zx * tex_scale;  
	vec2 coord3 = pos.xy * tex_scale;  

	// Sample color maps for each projection, at those UV coords.  
	vec4 col1 = texture2D(tex0, coord1);
	vec4 col2 = texture2D(tex0, coord2);
	vec4 col3 = texture2D(tex0, coord3);

	// Finally, blend the results of the 3 planar projections.  
	return (col1*blend_weights.xxxx + col2*blend_weights.yyyy + col3*blend_weights.zzzz);
}
