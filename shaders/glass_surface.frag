uniform vec3 light_color = vec3(1.0);
uniform sampler2D reflection_tex;

in vec4 epos, proj_pos;

void main() {
#ifdef ENABLE_REFLECTION
	vec3 normal     = normalize(fg_NormalMatrix * vec3(0.0, 0.0, 1.0)); // +Z in eye space
	vec3 epos_n     = normalize(epos.xyz); // eye direction
	float reflect_w = get_fresnel_reflection(-epos_n, normal, 1.0, 1.333);
	vec2 uv = clamp((0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5)), 0.0, 1.0);
	fg_FragColor = mix(gl_Color, texture(reflection_tex, uv), reflect_w);
#else
	fg_FragColor = gl_Color;
#endif
	fg_FragColor.xyz *= light_color;
}
