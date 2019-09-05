uniform samplerCube tex0;
//uniform float min_alpha = 0.0;

in vec3 dir;

void main() {
	vec4 texel = texture(tex0, dir);
	//if (texel.a <= min_alpha) discard; // no alpha test for cube maps?
	fg_FragColor = gl_Color*texel;
}
