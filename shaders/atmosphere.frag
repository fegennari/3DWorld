uniform vec3 camera_pos, planet_pos, sun_pos, ss_pos;
uniform float planet_radius, atmos_radius, sun_radius, ss_radius;
uniform vec3 light_scale   = vec3(1.0);
uniform vec3 atmos_density = vec3(1.0, 0.0, 0.0); // {constant, linear, quadratic}
uniform vec3 inner_color, outer_color;

in vec4 epos;
in vec3 normal, world_space_pos;

float get_density_at(vec3 pos) {
	float dist = distance(pos, planet_pos);
	float v = 1.0 - clamp((dist - planet_radius)/(atmos_radius - planet_radius), 0.0, 1.0);
	return dot(atmos_density, vec3(1.0, v, v*v));
}

void main()
{
	vec3 norm_norm = normalize(normal);
	vec3 light_dir = normalize(fg_LightSource[0].position.xyz - epos.xyz);
	float ascale   = min(4.0*(dot(norm_norm, light_dir) + 0.25), 1.0);
	if (ascale <= 0.0) discard;

	// alpha is calculated from distance between sphere intersection points
	float wpdist   = distance(world_space_pos, planet_pos);
	vec3 ldir      = normalize(world_space_pos - camera_pos);
	float dp       = dot(ldir, (world_space_pos - planet_pos));
	float adist_sq = dp*dp - wpdist*wpdist + atmos_radius*atmos_radius;
	if (adist_sq <= 0.0) discard; // no sphere intersection
	float dist     = sqrt(adist_sq);
	float pdist_sq = dp*dp - wpdist*wpdist + planet_radius*planet_radius;
	if (pdist_sq > 0.0) {dist -= sqrt(pdist_sq);} // ray intersects planet, adjust distance
	vec3 pos       = world_space_pos - ldir*(dp + 0.5*dist); // midpoint of ray in atmosphere
	float density  = get_density_at(pos)*dist/atmos_radius;
	float alpha    = ascale*clamp(4.0*density, 0.0, 1.0);
	float lt_atten = 1.0;

	if (sun_radius > 0.0 && ss_radius > 0.0) {
		lt_atten *= calc_sphere_shadow_atten(world_space_pos, sun_pos, sun_radius, ss_pos, ss_radius);
	}
	// Note: since only moons have a light2 set (from planet reflections), and moons have no atmosphere, light2 is not used here
	vec3 color = vec3(0.0);
	color += lt_atten*light_scale[0]*add_pt_light_comp(norm_norm, epos, 0).rgb; // sun ADS
	color += light_scale[1]*(gl_Color * fg_LightSource[1].ambient).rgb; // ambient only
	vec3 scatter_color = mix(outer_color, inner_color, min(1.6*density, 1.0)); // precomputed texture lookup
	fg_FragColor = vec4(color*scatter_color, alpha);
}
