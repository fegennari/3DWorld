#include <fresnel.part>

//uniform vec3 camera_pos; // from dynamic_lighting.part
uniform float water_depth, water_atten;
uniform vec3 uw_atten_max, uw_atten_scale;
uniform sampler2D reflection_tex;

in vec3 vertex_ws; // world space
in vec4 epos, proj_pos;

// underwater attenuation code
void atten_color(inout vec3 color, in float atten) {
	color *= vec3(1.0) - min(uw_atten_max, uw_atten_scale*atten*2.0);
}

void main() {
	vec3 ws_normal = vec3(0.0, 0.0, 1.0); // +Z
	vec3 normal = normalize(fg_NormalMatrix[2]); // eye space +z in world space
	vec3 epos_n = normalize(epos.xyz);
	vec2 uv     = clamp((0.5*proj_pos.xy/proj_pos.w + vec2(0.5, 0.5)), 0.0, 1.0);
	// TODO: apply ripples when the player moves
	vec4 reflect_tex = texture(reflection_tex, uv);
	vec3 water_color = vec3(0.0);
	if (enable_dlights) {add_dlights(water_color.rgb, vertex_ws, epos, ws_normal, gl_Color.rgb);} // dynamic lighting
	
	// attenuate based on depth/distance, assuming a constant water depth and average light coming from above in +Z
	vec3 floor_pos = vertex_ws - vec3(0.0, 0.0, water_depth);
	float t     = clamp(water_depth/max(0.001, abs(camera_pos.z - floor_pos.z)), 0.0, 1.0);
	float dist  = water_depth + t*distance(camera_pos, floor_pos); // light=>floor + camera=>floor through water
	float atten = water_atten*dist;
	//float alpha = gl_Color.a;
	float alpha = min(atten, 1.0);
	atten_color(water_color, atten);

	float reflect_w  = get_fresnel_reflection(-epos_n, normal, 1.0, 1.333);
	fg_FragColor     = mix(vec4(water_color, alpha), reflect_tex, reflect_w);
}
