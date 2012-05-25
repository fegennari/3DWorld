// input: point.xyz, normal.xyz, tangent.xyz, color
// output: textured quad as two triangles

void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	vec4 pos = gl_PositionIn[0]; // center of the quad
	vec3 normal  = gl_TexCoordIn[0][6].xyz;
	vec3 tangent = gl_TexCoordIn[0][7].xyz; // not normalized - provides the size of the quad (half-width)
	vec3 binorm  = cross(tangent, normal);
	vec4 pts[4];
	pts[0] = gl_ModelViewProjectionMatrix * (pos + vec4(tangent, 0.0));
	pts[1] = gl_ModelViewProjectionMatrix * (pos + vec4(tangent, 0.0) + vec4(binorm, 0.0));
	pts[2] = gl_ModelViewProjectionMatrix * (pos);
	pts[3] = gl_ModelViewProjectionMatrix * (pos + vec4(binorm,  0.0));
	output_textured_quad(pts);
}