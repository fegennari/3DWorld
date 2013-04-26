// input: line as two points {point.xyz, radius, color [,next_dir.xyz]}, ndiv
// output: textured cylinder/cone as triangles
varying vec3 normal; // world space
varying vec4 epos;
varying vec3 eye_norm;
varying vec3 vpos, spos, lpos0, vposl; // world space, not used

void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	vec4 pos1 = gl_PositionIn[0]; // world space
	vec4 pos2 = gl_PositionIn[1];

	vec3 v12 = (pos2 - pos1).xyz;
	vec3 vt = v12;
	vec3 vab[2];
	int min_dim = ((abs(v12.x) < abs(v12.y)) ? ((abs(v12.x) < abs(v12.z)) ? 0 : 2) : ((abs(v12.y) < abs(v12.z)) ? 1 : 2));
	vt[min_dim] += 0.5; // add to min dim
	vab[0] = normalize(cross(vt,  v12   )); // vab[0] is orthogonal to v12
	vab[1] = normalize(cross(v12, vab[0])); // vab[1] is orthogonal to v12 and vab[0]

	vec2 r = vec2(gl_TexCoordIn[0][7].r, gl_TexCoordIn[1][7].r); // radius packed into tex coord 7
	vec4 ce[2]; ce[0] = pos1; ce[1] = pos2;
	float PI  = 3.14159;
	float css = 2.0*PI/float(ndiv);
	float sin_ds = sin(css);
	float cos_ds = cos(css);
	float sin_s = 0.0;
	float cos_s = 1.0;
	// sin(x + y) = sin(x)*cos(y) + cos(x)*sin(y)
	// cos(x + y) = cos(x)*cos(y) - sin(x)*sin(y)

	for (int s = 0; s <= ndiv; ++s) {
		float sin_val = (sin_s*cos_ds + cos_s*sin_ds);
		float cos_val = (cos_s*cos_ds - sin_s*sin_ds);

		for (int i = 0; i < 2; ++i) { // two ends
			vec3 delta  = vab[0]*sin_s + vab[1]*cos_s; // world space
			epos        = gl_ModelViewMatrix  * (ce[i] + r[i]*vec4(delta, 0.0)); // eye space
			gl_FogFragCoord = length(epos.xyz); // set standard fog coord
			gl_Position = gl_ProjectionMatrix * epos;
			normal      = normalize(delta); // world space
			eye_norm    = normalize(gl_NormalMatrix * normal); // eye space
			gl_TexCoord[0].st = vec2(float(s)/float(ndiv), float(i));
			EmitVertex();
		}
		sin_s = sin_val;
		cos_s = cos_val;
	}
	EndPrimitive();
}