// input: line as two points {point.xyz, radius, color [,next_dir.xyz]}, ndiv
// output: textured cylinder/cone as triangles
uniform int ndiv = 12;
varying out vec3 normal;

void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	vec4 pos1 = gl_ModelViewMatrix * gl_PositionIn[0];
	vec4 pos2 = gl_ModelViewMatrix * gl_PositionIn[1];

	vec3 v12 = pos2 - pos1;
	vec3 vt = v12;
	vec3 vab[2];
	((abs(v12.x) < abs(v12.y)) ? ((abs(v12.x) < abs(v12.z)) ? vt.x : vt.z) : ((abs(v12.y) < abs(v12.z)) ? vt.y : vt.z)) += 0.5; // add to min dim
	vab[0] = normalize(cross(vt,  v12   )); // vab[0] is orthogonal to v12
	vab[1] = normalize(cross(v12, vab[0])); // vab[1] is orthogonal to v12 and vab[0]

	float r[2] = {r1, r2}; // FIXME: where do these come from?
	vec3 ce[2] = {pos0, pos1};
	float PI  = 3.14159;
	float css = 2.0*PI/float(ndiv);
	float sin_ds = sin(css);
	float cos_ds = cos(css);
	float sin_s = 0.0;
	float cos_s = 1.0;
	// sin(x + y) = sin(x)*cos(y) + cos(x)*sin(y)
	// cos(x + y) = cos(x)*cos(y) - sin(x)*sin(y)

	for (int s = 0; s < ndiv; ++s) {
		float sin_vals[2] = {sin_s, (sin_s*cos_ds + cos_s*sin_ds)}; // cur, next
		float cos_vals[2] = {cos_s, (cos_s*cos_ds - sin_s*sin_ds)}; // cur, next
		vec4 pts[4];
		vec3 n[4];
		vec2 tc[2];

		for (int i = 0; i < 2; ++i) { // two ends
			for (int j = 0; j < 2; ++j) { // prev and cur edges
				int ix = 2*i + j;
				vec3 delta = vab[0]*sin_vals[j] + vab[1]*cos_vals[j];
				pts[ix]   = ce[i] + r[i]*delta;
				n  [ix]   = normalize(delta);
				tc [ix].s = float(s+j)/float(ndiv);
				tc [ix].t = float(i);
			}
		}
		output_quad(pts, n, tc);
		sin_s = sin_vals[1];
		cos_s = cos_vals[1];
	}
}