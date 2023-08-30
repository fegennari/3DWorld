#include <fresnel.part>

uniform sampler2D reflection_tex;

in vec3 vertex_from_vs;
in vec4 epos, proj_pos;

void main() {
	vec3 ws_normal = vec3(0.0, 0.0, 1.0); // +Z
	vec3 normal = normalize(fg_NormalMatrix[2]); // eye space +z in world space
	vec3 epos_n = normalize(epos.xyz);
	vec2 uv     = clamp((0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5)), 0.0, 1.0);
	// TODO: apply ripples when the player moves
	vec4 reflect_tex = texture(reflection_tex, uv);
	vec3 water_color = vec3(0.0);
	if (enable_dlights) {add_dlights(water_color.rgb, vertex_from_vs, epos, ws_normal, gl_Color.rgb);} // dynamic lighting
	// TODO: attenuate based on depth/distance
	float reflect_w  = get_fresnel_reflection(-epos_n, normal, 1.0, 1.333);
	fg_FragColor     = mix(vec4(water_color, gl_Color.a), reflect_tex, reflect_w);
}
