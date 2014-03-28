uniform vec4 xlate_scale = vec4(0,0,0,1); // encoded as {tx, ty, tz, scale}

void main()
{
	vec4 vertex   = (vec4(xlate_scale.xyz, 0.0) + (vec4(xlate_scale.www, 1.0) * gl_Vertex));
	gl_Position   = gl_ModelViewProjectionMatrix * vertex;
	gl_FrontColor = vec4(0,0,0,1); // black/unused
}
