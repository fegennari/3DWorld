uniform float normal_z = 1.0;

varying vec3 normal;
varying vec4 epos, proj_pos;

void main()
{
	setup_texgen0();
	tc2    = gl_MultiTexCoord0;
	normal = gl_NormalMatrix * vec3(0.0, 0.0, normal_z); // not normalized
	epos   = gl_ModelViewMatrix * gl_Vertex;
	proj_pos = ftransform();
	gl_Position = proj_pos;
}
