#ifdef USE_BUMP_MAP
uniform sampler2D bump_map;

varying vec3 tangent_v, binorm_v;
varying vec3 ts_pos; // tangent space pos

// Note: we use tex coord 0 for the bump map (assuming it's the same tex coord for the regular texture)
vec3 apply_bump_map(inout vec3 eye_norm, inout vec3 light_dir)
{
	vec3 new_ldir;
	new_ldir.x = dot(light_dir, tangent_v);
	new_ldir.y = dot(light_dir, binorm_v);
	new_ldir.z = dot(light_dir, eye_norm);
	light_dir  = new_ldir;
	vec3 bump  = normalize(texture2D(bump_map, gl_TexCoord[0].st).xyz * 2.0 - 1.0);
	return bump;
}
#endif
