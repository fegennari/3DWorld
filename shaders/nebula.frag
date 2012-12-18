uniform vec3 view_dir;
uniform float radius;
uniform float alpha_scale = 1.0;
varying vec3 normal, vertex; // local object space

// Note: the nebula center is always assumed to be at 0,0,0 in local object space
void main()
{
	float alpha = gl_Color.a * alpha_scale;

	if (!line_mode) {
		float dist = length(vertex)/radius; // 0 at center, 1 at edge
		if (dist > 1.0) discard;
		alpha *= clamp(sqrt(1.0 - dist), 0.0, 1.0); // attenuate near edges to create a spherical shape
		alpha *= gen_cloud_alpha(vertex);
	}
	alpha *= abs(dot(normal, view_dir)); // attenuate billboards not facing the camera
	gl_FragColor = vec4(gl_Color.rgb, alpha);
}
