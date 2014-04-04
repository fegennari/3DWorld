layout (triangles) in;
layout (triangle_strip, max_vertices=3) out;

varying vec2 out tc;

void main()
{
	gl_FrontColor = gl_in[0].gl_FrontColor; // all colors are the same
	vec4 p0 = gl_in[0].gl_Position;
	vec4 p1 = gl_in[1].gl_Position;
	vec4 p2 = gl_in[2].gl_Position;
	
	gl_Position = p0;
	tc          = vec2(0.9, 0.1);
	EmitVertex();
	
	gl_Position = p1;
	tc          = vec2(0.9, 0.9);
	EmitVertex();
	
	// top vertex wind movement
	vec3  wind  = normalize(vec3(wind_x, wind_y, 0.0));
	float width = length(p1 - p0);
	vec3 normal = normalize(cross((p1.xyz - p0.xyz), (p2.xyz - p0.xyz)));
	float delta = width * abs(dot(wind, normal)) * get_wind_delta(p2.xyz, 1.0);
	vec4 offset = vec4(delta*wind_x, delta*wind_y, 0.0, 0.0);
	
	gl_Position = p2 + offset;
	tc          = vec2(0.1, 0.5);
	EmitVertex();
}
