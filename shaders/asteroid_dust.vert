uniform mat4 world_space_mvm;
uniform float sphere_size = 1.0;

varying vec3 world_space_pos;
varying vec4 epos;

void main()
{
	epos = gl_ModelViewMatrix * gl_Vertex;
#ifdef ENABLE_SHADOWS
	world_space_pos = (inverse(world_space_mvm) * epos).xyz;
#endif
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
#ifdef DRAW_AS_SPHERES
	float radius  = sphere_size*(0.5 + fract(223*gl_Vertex.x + 247*gl_Vertex.y + 262*gl_Vertex.z)); // random radius 0.5-1.5 * sphere_size
	gl_PointSize  = clamp(radius/length(epos.xyz), 1.0, 64.0);
#endif
}
