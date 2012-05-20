// input: point.xyz, normal.xyz, tangent.xyz, color
// output: textured quad as two triangles
varying in vec3 normal, tangent; // tangent not normalized - provides the size of the quad (half-width)

void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	//normal_out = normalize(normal); // not needed?
	vec4 pos = gl_ModelViewMatrix * gl_PositionIn[0]; // center of the quad
	vec4 pts[4];
	vec3 binorm = cross(normal, tangent);
	pts[0] = gl_ProjectionMatrix * (pos - vec4(tangent, 0.0) + vec4(binorm, 0.0));
	pts[1] = gl_ProjectionMatrix * (pos - vec4(tangent, 0.0) - vec4(binorm, 0.0));
	pts[2] = gl_ProjectionMatrix * (pos + vec4(tangent, 0.0) + vec4(binorm, 0.0));
	pts[3] = gl_ProjectionMatrix * (pos + vec4(tangent, 0.0) - vec4(binorm, 0.0));
	output_textured_quad(pts);
}