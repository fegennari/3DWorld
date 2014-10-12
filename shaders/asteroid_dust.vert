uniform mat4 fg_ViewMatrix;
uniform float sphere_size = 1.0;
uniform vec4 color = vec4(1.0);

out vec3 world_space_pos;
out vec4 epos;

void main()
{
	epos = fg_ModelViewMatrix * fg_Vertex;
#ifdef ENABLE_SHADOWS
	world_space_pos = (inverse(fg_ViewMatrix) * epos).xyz;
#endif
	gl_Position   = fg_ProjectionMatrix * epos;
	gl_FrontColor = color;
#ifdef DRAW_AS_SPHERES
	float radius  = sphere_size*(0.5 + fract(223*fg_Vertex.x + 247*fg_Vertex.y + 262*fg_Vertex.z)); // random radius 0.5-1.5 * sphere_size
	float pt_sz   = radius/length(epos.xyz);
	gl_PointSize  = clamp(pt_sz, 1.0, 64.0);
	gl_FrontColor.a *= min(1.0, 3.0*pt_sz); // attenuate very small points
#else
	gl_PointSize  = 1.0;
#endif
}
