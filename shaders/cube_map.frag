uniform samplerCube tex0;
//uniform float min_alpha = 0.0;

in vec3 dir;

void main() {
	// all the cube map textures I have use +y for up, but 3DWorld has +z for up, so we need to do some vector swizzles here
	vec4 texel = texture(tex0, -dir.yzx);
	//if (texel.a <= min_alpha) discard; // no alpha test for cube maps?
	fg_FragColor = gl_Color*texel;
}
