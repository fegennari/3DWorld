layout (points) in;
layout (triangle_strip, max_vertices=4) out;

// input: point.xyz, normal.xyz, tangent.xyz, color
// output: textured quad as two triangles
in VertexData {
    vec3 normal;
    vec3 tangent;
	vec4 color;
} VertexIn[1]; // from VS

void main()
{
	vec4 pos     = gl_in[0].gl_Position; // center of the quad
	vec3 normal  = vertexIn[0].normal;
	vec4 tangent = vec4(vertexIn[0].tangent, 0.0); // not normalized - provides the size of the quad (half-width)
	vec4 binorm  = vec4(cross(tangent.xyz, normal), 0.0);
	vec4 pts[4];
	pts[0] = fg_ModelViewProjectionMatrix * (pos);
	pts[1] = fg_ModelViewProjectionMatrix * (pos + binorm);
	pts[2] = fg_ModelViewProjectionMatrix * (pos + tangent);
	pts[3] = fg_ModelViewProjectionMatrix * (pos + tangent + binorm);
	output_textured_quad(pts, vertexIn[0].color);
}
