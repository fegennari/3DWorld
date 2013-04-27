uniform vec3 camera_pos, planet_pos, sun_pos, ss_pos;
uniform float planet_radius, atmos_radius, sun_radius, ss_radius;
uniform float atmosphere = 1.0;
uniform vec3 light_scale = vec3(1.0);
varying vec4 epos;
varying vec3 normal, world_space_pos;

void main()
{
	// alpha is calculated from distance between sphere intersection points
	float wpdist   = length(world_space_pos - planet_pos);
	float dp       = dot(normalize(world_space_pos - camera_pos), (world_space_pos - planet_pos));
	float adist_sq = dp*dp - wpdist*wpdist + atmos_radius*atmos_radius;
	if (adist_sq <= 0.0) discard; // no sphere intersection
	float dist     = sqrt(adist_sq);
	float pdist_sq = dp*dp - wpdist*wpdist + planet_radius*planet_radius;
	if (pdist_sq > 0.0) {dist -= sqrt(pdist_sq);} // ray intersects planet, adjust distance
	float density  = dist/atmos_radius;
	float alpha    = atmosphere*clamp(4.0*density, 0.0, 1.0);
	float lt_atten = 1.0;

	if (sun_radius > 0.0 && ss_radius > 0.0) {
		lt_atten *= calc_sphere_shadow_atten(world_space_pos, sun_pos, sun_radius, ss_pos, ss_radius);
	}
	// Note: since only moons have a light2 set (from planet reflections), and moons have no atmosphere, light2 is not used here
	vec4 color  = gl_FrontMaterial.emission;
	color.rgb  += lt_atten*light_scale[0]*add_pt_light_comp(normalize(normal), epos, 0).rgb; // sun ADS
	color.rgb  += light_scale[1]*(gl_Color * gl_LightSource[1].ambient).rgb; // ambient only
	float rg_comp = atmosphere*min(1.6*density, 1.0);
	vec3 scatter_color = vec3(rg_comp, rg_comp, 1.0); // precomputed texture lookup or something else better?
	color *= vec4(scatter_color, alpha);
	gl_FragColor = color;
}
