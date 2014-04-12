uniform mat4 world_space_mvm;

varying vec4 epos;
varying vec3 eye, dlpos, normal; // world space
varying vec3 dl_normal;
varying vec2 tc;

void main()
{
	tc        = fg_TexCoord;
	epos      = gl_ModelViewMatrix * fg_Vertex;
	normal    = normalize(gl_NormalMatrix * fg_Normal);
	dl_normal = normalize((transpose(world_space_mvm) * vec4(normal, 1)).xyz);
	mat4 mvm_inv = inverse(world_space_mvm);
	dlpos     = (mvm_inv * epos).xyz;
	eye       = mvm_inv[3].xyz;
	gl_Position     = gl_ProjectionMatrix * epos;
	gl_FogFragCoord = length(epos.xyz);
}
