uniform mat4 world_space_mvm;

varying vec4 epos;
varying vec3 eye, dlpos, normal; // world space
varying vec3 dl_normal;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	epos      = gl_ModelViewMatrix * gl_Vertex;
	normal    = normalize(gl_NormalMatrix * gl_Normal);
	dl_normal = normalize((transpose(world_space_mvm) * vec4(normal, 1)).xyz);
	mat4 mvm_inv = inverse(world_space_mvm);
	dlpos     = (mvm_inv * epos).xyz;
	eye       = mvm_inv[3].xyz;
	gl_Position = ftransform();
	set_fog();
}
