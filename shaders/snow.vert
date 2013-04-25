uniform mat4 world_space_mvm;
varying vec4 epos;
varying vec3 dlpos, dl_normal; // world space
varying vec3 normal;

void main()
{
	setup_texgen0();
	dl_normal = normalize(gl_Normal);
	vec3 n = gl_NormalMatrix * gl_Normal;
	normal = (no_normalize ? n : normalize(n));
	epos   = gl_ModelViewMatrix * gl_Vertex;
	dlpos  = gl_Vertex.xyz;
	gl_Position = ftransform();
	gl_FogFragCoord = length(epos.xyz);
}
