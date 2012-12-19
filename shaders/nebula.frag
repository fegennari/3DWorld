uniform vec4 secondary_color;
uniform vec3 view_dir;
uniform float radius;
uniform float offset = 0.0;
varying vec3 normal, vertex; // local object space

// Note: the nebula center is always assumed to be at 0,0,0 in local object space
void main()
{
	vec4 color = gl_Color;

	if (!line_mode) {
		float dist = length(vertex)/radius; // 0 at center, 1 at edge
		if (dist > 1.0) discard;
		color.a *= gen_cloud_alpha(vertex + vec3(offset, offset, offset));
		float a2 = gen_cloud_alpha(vertex - vec3(offset, offset, offset) + vec3(0.1, 0.4, 0.7));
		color    = a2*secondary_color + (1.0 - a2)*color;
		color.a *= clamp(sqrt(1.0 - dist), 0.0, 1.0); // attenuate near edges to create a spherical shape
	}
	color.a *= clamp((1.5*abs(dot(normal, view_dir)) - 0.5), 0.0, 1.0); // attenuate billboards not facing the camera
	gl_FragColor = color;
}
