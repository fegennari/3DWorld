uniform float fog_scale = 0.0;

// linear fog
vec4 apply_fog_ffc(in vec4 color, in float ffc) {
	float fog = clamp((gl_Fog.end - ffc) * gl_Fog.scale, 0.0, 1.0);
	//return vec4(mix(gl_Fog.color.rgb, color.rgb, mix(1.0, fog, fog_scale)), color.a);
	return mix(gl_Fog.color, color, mix(1.0, fog, fog_scale));
}

vec4 apply_fog(in vec4 color) {
	return apply_fog_ffc(color, gl_FogFragCoord);
}

vec4 apply_fog_epos(in vec4 color, in vec4 epos) {
	return apply_fog_ffc(color, length(epos.xyz));
}
