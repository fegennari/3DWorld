uniform mat4 fg_ViewMatrix;
uniform float sphere_size = 1.0;
uniform vec4 color = vec4(1.0);

varying vec3 world_space_pos;
varying vec4 epos;

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
	gl_PointSize  = clamp(radius/length(epos.xyz), 1.0, 64.0);
#endif
}
