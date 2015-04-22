uniform vec3 sun_pos, camera;
uniform sampler2D tex0;

in vec2 tc;

void main()
{
	// Note: constant, could be moved into vertex shader, but not performance critical
	float alpha  = clamp(1.1*get_cloud_plane_alpha(camera, vec4(sun_pos, 1.0)), 0.0, 1.0);
	fg_FragColor = gl_Color*vec4((1.0 - alpha)*texture(tex0, tc).rrr, 1.0); // since this uses additive blending, we modify RGB instead of A
}

