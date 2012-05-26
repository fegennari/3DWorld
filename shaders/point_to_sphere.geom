// input: point.xyz, radius, color, ndiv
// output: textured sphere as triangles
uniform int ndiv;
varying in float radius;

void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	vec4 pos = gl_ModelViewMatrix * gl_PositionIn[0]; // center of the sphere
	float ndiv_inv = 1.0/ndiv;

	for (int s = 0; s < ndiv; ++s) {
		for (int t = 0; t <= ndiv; ++t) {
			for (int i = 0; i < 2; ++i) {
				vec3 delta = ?; // FIXME
				gl_Position = pos + radius*vec4(delta, 0.0);
				normal      = normalize(delta);
				gl_TexCoord[0].st = vec2((1.0 - (s+i)*ndiv_inv), (1.0 - t*ndiv_inv));
				EmitVertex();
			}
		}
		EndPrimitive();
	}
}
