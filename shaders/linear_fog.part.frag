uniform float fog_scale = 0.0;

// linear fog
vec4 apply_fog(in vec4 color) {
	float fog = clamp((gl_Fog.end - gl_FogFragCoord) * gl_Fog.scale, 0.0, 1.0);
	return mix(gl_Fog.color, color, mix(1.0, fog, fog_scale));
}
