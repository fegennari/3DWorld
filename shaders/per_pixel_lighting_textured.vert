uniform mat4 fg_ViewMatrix;

varying vec4 epos;
varying vec3 eye, dlpos, normal; // world space
varying vec3 dl_normal;
varying vec2 tc;

void main()
{
	tc        = fg_TexCoord;
	epos      = fg_ModelViewMatrix * fg_Vertex;
	normal    = normalize(fg_NormalMatrix * fg_Normal);
	dl_normal = normalize((transpose(fg_ViewMatrix) * vec4(normal, 1)).xyz);
	mat4 vm_inv = inverse(fg_ViewMatrix);
	dlpos     = (vm_inv * epos).xyz;
	eye       = vm_inv[3].xyz;
	gl_Position     = fg_ProjectionMatrix * epos;
	gl_FogFragCoord = length(epos.xyz);
}
