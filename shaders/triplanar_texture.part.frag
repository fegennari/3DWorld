uniform sampler2D tex0;

vec4 lookup_triplanar_texture(in vec3 normal) {
	vec2 tc = vec2(0.0, 0.0); // FIXME
	return texture2D(tex0, tc);
}
