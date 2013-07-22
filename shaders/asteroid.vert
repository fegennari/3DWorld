varying vec3 vpos, normal, world_normal;
varying vec4 epos;
attribute mat4 inst_xform_matrix;

void main()
{
#ifdef USE_CUSTOM_XFORM
	gl_Position   = (gl_ProjectionMatrix * inst_xform_matrix) * gl_Vertex;
	normal        = normalize(transpose(inverse(mat3(inst_xform_matrix))) * gl_Normal); // for lighting (FIXME: incorrect?)
#else
	gl_Position   = ftransform();
	normal        = normalize(gl_NormalMatrix * gl_Normal); // for lighting
#endif
	gl_FrontColor = gl_Color;
	world_normal  = gl_Normal; // for triplanar texturing
	epos          = gl_ModelViewMatrix * gl_Vertex;
	vpos          = gl_Vertex.xyz;
} 
