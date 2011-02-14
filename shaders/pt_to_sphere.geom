varying vec3 normal;

void add_vert(in vec3 p, in vec3 n, in vec2 t)
{
	gl_Position       = vec4(p, 1.0);
	normal            = n;
	gl_TexCoord[0].st = t;
	EmitVertex();
}

void create_quad(in vec3 p[4], in vec3 n[4], in vec2 t[4])
{
	add_vert(p[0], n[0], t[0]);
	add_vert(p[1], n[1], t[1]);
	add_vert(p[2], n[2], t[2]);
	EndPrimitive();

	add_vert(p[1], n[1], t[1]);
	add_vert(p[2], n[2], t[2]);
	add_vert(p[3], n[3], t[3]);
	EndPrimitive();
}

void main()
{
	gl_FrontColor = gl_FrontColorIn[0];
	vec3 scale = gl_TexCoordIn[0][0].stq; // FIXME: hack
	vec4 pos = gl_PositionIn[0];

	// create sphere points
#define ndiv 4
	vec3 p[ndiv][ndiv+1];
	vec3 n[ndiv][ndiv+1];
	vec2 t[ndiv][ndiv+1];
	const float PI  = 3.1415926535;
	float ndiv_inv  = 1.0/ndiv;
	float cs_scale  = PI*ndiv_inv;
	float cs_scale2 = 2.0*cs_scale;
	float sin_dt    = sin(cs_scale);
	float cos_dt    = cos(cs_scale);

	for (uint s = 0; s < ndiv; ++s) { // build points and normals table
		float theta = s*cs_scale2;
		float tvc   = cos(theta); // theta1
		float tvs   = sin(theta); // theta2
		float sin_t = 0.0;
		float cos_t = 1.0;
		uint soff   = (ndiv+1)*s;

		for (uint t = 0; t <= ndiv; ++t) { // Note: x and y are swapped because theta is out of phase by 90 degrees to match tex coords
			float sv = sin_t;
			float cv = cos_t;
			n[s][t]  = vec3(sv*tvs, sv*tvc, cv); // R*sin(phi)*sin(theta), R*sin(phi)*cos(theta), R*cos(phi)
			sin_t    = sv*cos_dt + cv*sin_dt;
			cos_t    = cv*cos_dt - sv*sin_dt;
			p[s][t]  = (gl_ModelViewProjectionMatrix * (vec4(n[s][t]*scale, 0.0) + pos)).xyz;
		}
	}

	// create quads
	for (uint s = 0; s < ndiv; ++s) {
		uint sn  = (s+1)%ndiv;
		uint snt = min((s+1), ndiv);

		for (uint t = 0; t < ndiv; ++t) {
			uint ix = s*(ndiv+1)+t;
			uint tn = t+1;
			vec3 pts[4], norm[4];
			vec2 tex[4];
			pts [0] = p[s][t]; pts [1] = p[sn][t]; pts [2] = p[sn][tn]; pts [3] = p[s][tn];
			norm[0] = n[s][t]; norm[1] = n[sn][t]; norm[2] = n[sn][tn]; norm[3] = n[s][tn];

			for (uint i = 0; i < 4; ++i) {
				tex[i] = vec2(0.5,0.5);//vec2((1.0f - (((i&1)^(i>>1))!=0 ? snt : s)*ndiv_inv), (1.0f - ((i>>1)!=0 ? tn : t)*ndiv_inv));
			}
			create_quad(pts, norm, tex);
		}
	}
}