uniform float alpha_scale = 1.0;
uniform sampler2D tex0;

in vec3 world_space_pos;
in vec4 epos;

vec3 get_sphere_point_sprite_normal() {

	vec3 normal;
	normal.xy = 2.0*gl_PointCoord - vec2(1.0); 
	float mag = dot(normal.xy, normal.xy);
	if (mag > 1.0) discard; // outside circle/sphere radius
	normal.z  = sqrt(1.0 - mag);
	vec3 zv = normalize(-epos.xyz);
	vec3 yv = normalize(cross(vec3(1,0,0), zv));
	vec3 xv = normalize(cross(zv, yv));
	return normalize(xv*normal.x + yv*normal.y + zv*normal.z);
}

vec2 get_sphere_uv(in vec3 normal) {
	return vec2(atan(normal.z, normal.x)/(3.14159*2.0), acos(-normal.y)/3.14159);
}

void main() {
	float alpha = min(1.0, (2.0/(1.0 + alpha_scale*length(epos.xyz)) - 0.5)); // attenuate when far from the camera
	if (alpha <= 0.0) {discard;}
#ifdef DRAW_AS_SPHERES
	vec3 normal = get_sphere_point_sprite_normal();
	//gl_FragDepth = ???;
#else
	vec3 normal = normalize(-epos.xyz); // facing the camera
#endif
	vec3 color   = vec3(0.0);
	float atten0 = 1.0;
#ifdef ENABLE_SHADOWS
	atten0 = calc_shadow_atten(world_space_pos);
#endif
	color += atten0*add_pt_light_comp(normal, epos, 0).rgb; // sun_diffuse
	color += add_pt_light_comp(normal, epos, 1).rgb; // galaxy_ambient
#ifdef DRAW_AS_SPHERES
	//color *= texture(tex0, gl_PointCoord).rgb;
	color *= texture(tex0, get_sphere_uv(normal)).rgb;
#endif
	fg_FragColor = vec4(color, alpha * gl_Color.a); // use gl_Color alpha directly
}

