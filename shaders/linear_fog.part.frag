uniform float fog_scale = 0.0;
uniform float fog_fade_val = 8.0;

// linear/quadratic fog
vec4 apply_fog_ffc(in vec4 color, in float ffc, in vec4 fog_color) {
	float fog = clamp((gl_Fog.end - ffc) * gl_Fog.scale, 0.0, 1.0);
#ifdef USE_QUADRATIC_FOG
	fog = 1.0 - (1.0-fog)*(1.0-fog); // quadratic term
#endif
	float fog_mix_val = mix(1.0, fog, fog_scale);
	//vec4 fin_color = vec4(mix(fog_color.rgb, color.rgb, fog_mix_val), color.a);
	vec4 fin_color = mix(fog_color, color, fog_mix_val); // fog_mix_val=0 => fog_color, fog_mix_val=1 => color
#ifdef FOG_FADE_TO_TRANSPARENT
	fin_color.a *= min(fog_fade_val*fog_mix_val, 1.0);
#endif
	return fin_color;
}

vec4 apply_fog(in vec4 color) {
	return apply_fog_ffc(color, gl_FogFragCoord, gl_Fog.color);
}

vec4 apply_fog_epos(in vec4 color, in vec4 epos) {
	return apply_fog_ffc(color, length(epos.xyz), gl_Fog.color);
}

vec4 apply_fog_colored(in vec4  color,    // original color of the pixel
                       in float distance, // camera to point distance
                       in vec3  ray_dir,  // camera to point vector
                       in vec3  sun_dir)  // sun light direction
{
#ifdef USE_QUADRATIC_FOG // hack to make this only apply to tiled terrain mode, which uses quadratic fog
    float sun_amount = ((sun_dir.z > 0.0) ? max(dot(ray_dir, sun_dir), 0.0) : 0.0);
	vec4 sun_color   = vec4(max(gl_Fog.color.r, gl_Fog.color.b), max(gl_Fog.color.g, gl_Fog.color.b), 0.5*(gl_Fog.color.r + gl_Fog.color.g), gl_Fog.color.a);
    vec4 fog_color   = mix(gl_Fog.color, sun_color, pow(sun_amount, 8.0));
	return apply_fog_ffc(color, distance, fog_color);
#else
	return apply_fog_ffc(color, distance, gl_Fog.color);
#endif
}

vec4 apply_fog_colored(in vec4 color, in vec4 epos) {
	return apply_fog_colored(color, gl_FogFragCoord, -normalize(epos.xyz), normalize(epos.xyz - gl_LightSource[0].position.xyz));
}

vec4 apply_fog_colored_epos(in vec4 color, in vec4 epos) {
	return apply_fog_colored(color, length(epos.xyz), -normalize(epos.xyz), normalize(epos.xyz - gl_LightSource[0].position.xyz));
}
