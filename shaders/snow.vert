varying vec4 epos;
varying vec3 eye, dlpos, normal; // world space
uniform mat4 world_space_mvm;
varying vec3 dl_normal;

void main()
{
	setup_texgen0();
	vec3 n = gl_NormalMatrix * gl_Normal;
	dl_normal = normalize(n);
	normal = (no_normalize ? n : dl_normal);
	epos   = gl_ModelViewMatrix * gl_Vertex;
	dlpos  = gl_Vertex.xyz;
	eye    = gl_ModelViewMatrixInverse[3].xyz; // world space
	gl_Position = ftransform();
	gl_FogFragCoord = length(epos.xyz);
}
