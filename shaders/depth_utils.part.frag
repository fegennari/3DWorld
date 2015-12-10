uniform float znear, zfar;
uniform vec2 xy_step;
uniform sampler2D depth_tex;

float log_to_linear_depth(in float d) {
	return 2.0 * zfar * znear / (zfar + znear - (zfar - znear)*(2*d - 1)); // actual z-value
}

float get_linear_depth_01(in vec2 pos) {
	float d = texture(depth_tex, pos).r;
	return (2.0 * znear) / (zfar + znear - d * (zfar - znear)); // [0,1] range
}

float get_linear_depth_zval(in vec2 pos) {
	return log_to_linear_depth(texture(depth_tex, pos).r);
}

float get_depth_at_fragment() {
	return texture(depth_tex, gl_FragCoord.xy*xy_step).r;
}

#if 0
// taken from https://mynameismjp.wordpress.com/page/6/
vec3 view_space_pos_from_depth(in vec2 pos) {
	float z = get_linear_depth_zval(pos); // get the depth value for this pixel
	float x = pos.x*2.0 - 1.0;
	float y = (1.0 - pos.y)*2.0 - 1.0;
	vec4 proj_pos = vec4(x, y, z, 1.0); // get x/w and y/w from the viewport position
	vec4 vpos = proj_pos * inverse(fg_ProjectionMatrix); // transform by the inverse projection matrix
	return vpos.xyz / vpos.w; // divide by w to get the view-space position
}
#endif

