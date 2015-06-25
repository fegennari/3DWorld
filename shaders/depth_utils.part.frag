uniform float znear, zfar;
uniform sampler2D depth_tex;

float get_linear_depth_01(in vec2 pos) {
	float d = texture(depth_tex, pos).r;
	return (2.0 * znear) / (zfar + znear - d * (zfar - znear)); // [0,1] range
}

float get_linear_depth_zval(in vec2 pos) {
	float d = texture(depth_tex, pos).r;
	return 2.0 * zfar * znear / (zfar + znear - (zfar - znear)*(2*d - 1)); // actual z-value
}
