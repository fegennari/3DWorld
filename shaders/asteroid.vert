varying vec3 vpos, normal, world_normal;
varying vec4 epos;

void main()
{
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
	world_normal  = gl_Normal;
	normal = normalize(gl_NormalMatrix * gl_Normal);
	epos   = gl_ModelViewMatrix * gl_Vertex;
	vpos   = gl_Vertex.xyz;
} 
