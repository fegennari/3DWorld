uniform float alpha_scale = 1.0;

varying vec3 world_space_pos;
varying vec4 epos;

void main()
{
	float alpha = min(1.0, (2.0/(1.0 + alpha_scale*length(epos.xyz)) - 0.5)); // attenuate when far from the camera
	if (alpha <= 0.0) {discard;}
#ifdef DRAW_AS_SPHERES
	vec3 normal;
	normal.xy = 2.0*gl_PointCoord - vec2(1.0); 
	float mag = dot(normal.xy, normal.xy);
	if (mag > 1.0) discard; // outside circle/sphere radius
	normal.z  = sqrt(1.0 - mag);
#else
	vec3 normal = normalize(-epos.xyz); // facing the camera
#endif
	vec4 color     = gl_FrontMaterial.emission;
	float atten[2] = {1.0, 1.0};
#ifdef ENABLE_SHADOWS
	atten[0] = calc_shadow_atten(world_space_pos);
#endif

	for (int i = 0; i < 2; ++i) { // sun_diffuse, galaxy_ambient
		color += atten[i]*add_pt_light_comp(normal, epos, i);
	}
	gl_FragColor = vec4(color.rgb, alpha * gl_Color.a); // use diffuse alpha directly;
}

