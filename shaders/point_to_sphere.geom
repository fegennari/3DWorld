// input: point.xyz, radius, color, ndiv
// output: textured sphere as triangles
uniform int ndiv;
varying in float radius;

void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	vec4 pos = gl_ModelViewMatrix * gl_PositionIn[0]; // center of the sphere

	for (int s = 0; s < ndiv; ++s) {
		for (int t = 0; t < ndiv; ++t) {
			vec4 pts[4];
			vec3 n[4];
			vec2 tc[2];
			// FIXME: Write
			output_quad(pts, n, tc);
		}
	}
}
