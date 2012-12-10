uniform vec3 camera_pos, planet_pos;
uniform float cloud_radius;
uniform float atmosphere = 1.0;
uniform vec2 light_scale = vec2(1,1);
varying vec4 epos;
varying vec3 normal, world_space_pos;

void main()
{
	// alpha is calculated from distance between sphere intersection points
	float pdist   = length(world_space_pos - planet_pos);
	float dp      = dot(normalize(world_space_pos - camera_pos), (world_space_pos - planet_pos));
	float dist_sq = dp*dp - pdist*pdist + cloud_radius*cloud_radius;
	if (dist_sq <= 0.0) discard;
	float density = atmosphere*sqrt(dist_sq)/cloud_radius;
	float alpha = min(4.0*density, 1.0);
	vec4 color  = gl_FrontMaterial.emission;
	color.rgb  += light_scale[0]*add_pt_light_comp(normalize(normal), epos, 0).rgb; // ambient, diffuse, and specular
	color.rgb  += light_scale[1]*(gl_Color * gl_LightSource[1].ambient).rgb; // ambient only
	float rg_comp = min(1.6*density, 1.0);
	vec3 scatter_color = vec3(rg_comp, rg_comp, 1.0); // precomputed texture lookup or something else better?
	color *= vec4(scatter_color, alpha);
	gl_FragColor = apply_fog(color);
}
