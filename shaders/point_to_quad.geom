// input: point.xyz, normal.xyz, tangent.xyz, color
// output: textured quad as two triangles

void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	vec4 pos = gl_PositionIn[0]; // center of the quad
	vec3 normal  = gl_TexCoordIn[0][6].xyz;
	vec4 tangent = vec4(gl_TexCoordIn[0][7].xyz, 0.0); // not normalized - provides the size of the quad (half-width)
	vec4 binorm  = vec4(cross(tangent.xyz, normal), 0.0);
	vec4 pts[4];
	pts[0] = gl_ModelViewProjectionMatrix * (pos);
	pts[1] = gl_ModelViewProjectionMatrix * (pos + binorm);
	pts[2] = gl_ModelViewProjectionMatrix * (pos + tangent);
	pts[3] = gl_ModelViewProjectionMatrix * (pos + tangent + binorm);
	output_textured_quad(pts);
}