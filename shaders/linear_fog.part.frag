uniform float fog_scale = 0.0;
uniform float fog_fade_val = 8.0;

// linear/quadratic fog
vec4 apply_fog_ffc(in vec4 color, in float ffc) {
	float fog = clamp((gl_Fog.end - ffc) * gl_Fog.scale, 0.0, 1.0);
#ifdef USE_QUADRATIC_FOG
	fog = 1.0 - (1.0-fog)*(1.0-fog); // quadratic term
#endif
	float fog_mix_val = mix(1.0, fog, fog_scale);
	//vec4 fin_color = vec4(mix(gl_Fog.color.rgb, color.rgb, fog_mix_val), color.a);
	vec4 fin_color = mix(gl_Fog.color, color, fog_mix_val); // fog_mix_val=0 => gl_Fog.color, fog_mix_val=1 => color
#ifdef FOG_FADE_TO_TRANSPARENT
	fin_color.a *= min(fog_fade_val*fog_mix_val, 1.0);
#endif
	return fin_color;
}

vec4 apply_fog(in vec4 color) {
	return apply_fog_ffc(color, gl_FogFragCoord);
}

vec4 apply_fog_epos(in vec4 color, in vec4 epos) {
	return apply_fog_ffc(color, length(epos.xyz));
}
