uniform sampler3D noise_tex;
uniform float noise_scale = 1.0;
uniform vec4 color1i, color2i, color3i, color1o, color2o, color3o;
uniform vec3 view_dir;
uniform float radius = 1.0;
uniform float offset = 0.0;
uniform float alpha_scale= 1.0;
uniform float alpha_bias = -0.4; // intended to be changed for grayscale mode
uniform float dist_bias  = 0.0;
uniform vec3 rscale      = vec3(1.0);
uniform vec4 color_mult  = vec4(1.0);

in vec3 normal, vertex; // local object space

// Note: the nebula center is always assumed to be at 0,0,0 in local object space
void main() {

	vec4 color   = color1o;
	float ascale = alpha_scale;
	ascale *= clamp((1.5*abs(dot(normal, view_dir)) - 0.5), 0.0, 1.0); // attenuate billboards not facing the camera

	if (!line_mode) {
		vec3 view_vect = vertex/(rscale*radius);
		float dist = length(view_vect); // 0 at center, 1 at edge
		if (dist > 1.0) discard;
#ifdef USE_CLOUD_MODE
		ascale    *= pow(max(1.0, (1.0 - view_vect.z)), 2.0); // clouds denser near the bottom
#endif
		ascale    *= clamp((1.0 - dist), 0.0, 1.0); // attenuate near edges to create a spherical shape
		color      = mix(color1i, color, dist);
		if (color.a*ascale < 0.002) discard;
		vec3  pos  = vertex + vec3(offset);
		vec4  val  = vec4(0.0);
		float freq = 1.0;

		for (int i = 0; i < NUM_OCTAVES; ++i) {
			vec4 v = texture(noise_tex, noise_scale*(freq*pos));
#ifdef RIDGED_NOISE
			v = 2.0*v - 1.0; // map [0,1] range to [-1,1]
			v = 1.0 - abs(v); // ridged noise
			v = v*v;
#endif
			val  += v/freq;
			freq *= 2.0;
		}
		val = clamp(2.0*(0.5*val + alpha_bias + dist_bias*dist), 0.0, 1.0);
		
		if (noise_ncomp == 1) { // grayscale
			color.a *= val.r;
		}
		else { // RGBA
			color.a *= val.a;
			color = mix(color, mix(color2i, color2o, dist), abs(val.r - val.g));
			color = mix(color, mix(color3i, color3o, dist), abs(val.r - val.b));
		}
	} // end !line_mode
	color.a *= min(ascale, 1.0);
	if (color.a < 0.002) discard;
	//color.a *= min(1.0, 10.0*length((fg_ModelViewMatrix * vec4(vertex, 1.0)).xyz)/radius); // atten when very close to a plane (based on calculated epos)

#ifdef ENABLE_LIGHTING
	vec4 epos  = fg_ModelViewMatrix * vec4(vertex, 1.0);
#ifdef USE_CLOUD_MODE
	float dp   = dot(normalize(epos.xyz - fg_LightSource[0].position.xyz).xyz, normalize(epos.xyz));
	color.rgb *= 1.0 + 0.2*dp + 0.7*pow(max(-dp, 0.0), 30.0);
#else // !USE_CLOUD_MODE
	color.rgb *= (fg_LightSource[0].diffuse.rgb + fg_LightSource[1].ambient.rgb) * calc_light_atten(epos, 0); // sun_diffuse + galaxy_ambient
#endif // USE_CLOUD_MODE
#endif // ENABLE_LIGHTING
	fg_FragColor = color_mult * color;
}
