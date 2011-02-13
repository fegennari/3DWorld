vec4 get_color(in vec4 base_color, in vec3 n)
{
	vec4 color = base_color;
	if (enable_light0) color += add_light_comp(n, 0);
	if (enable_light1) color += add_light_comp(n, 1);
	return color;
}

void add_vert(in vec4 p, in vec4 c, in vec2 t)
{
	gl_Position       = p;
	gl_FrontColor     = c;
	gl_TexCoord[0].st = t;
	EmitVertex();
}

void create_quad(in vec4 p[4], in vec4 c[4], in vec2 t[4])
{
	add_vert(p[0], c[0], t[0]);
	add_vert(p[1], c[1], t[1]);
	add_vert(p[2], c[2], t[2]);
	EndPrimitive();

	add_vert(p[1], c[1], t[1]);
	add_vert(p[2], c[2], t[2]);
	add_vert(p[3], c[3], t[3]);
	EndPrimitive();
}

void main()
{
	const uint ndiv = 6;
	vec4 base_color = gl_FrontColorIn[0] * gl_LightModel.ambient;
	vec3 scale = gl_TexCoordIn[0][0].stq; // FIXME: hack
	vec4 pos = gl_PositionIn[0];

	vec4 p[4], c[4];
	vec2 t[4];
	vec4 dx = vec4(0,0,1,1); // 01 00 11 10
	vec4 dz = vec4(1,0,1,0);

	for (int i = 0; i < 4; ++i) {
		vec3 n = vec3(0,0,1);
		p[i] = gl_ModelViewProjectionMatrix * (pos + vec4(scale.x*(2*dx[i]-1), 0, scale.z*(2*dz[i]-1), 0));
		c[i] = get_color(base_color, n);
		t[i] = vec2(float(i>>1), float(!(i&1)));
	}
	create_quad(p, c, t);
}