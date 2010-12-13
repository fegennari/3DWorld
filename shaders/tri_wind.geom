void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	vec4 p0 = gl_PositionIn[0];
	vec4 p1 = gl_PositionIn[1];
	vec4 p2 = gl_PositionIn[2];
	
	gl_Position = p0;
	gl_TexCoord[0].st = vec2(0.9, 0.1);
	EmitVertex();
	
	gl_Position = p1;
	gl_TexCoord[0].st = vec2(0.9, 0.9);
	EmitVertex();
	
	// top vertex wind movement
	vec3  wind  = normalize(vec3(wind_x, wind_y, 0.0));
	float width = length(p1 - p0);
	vec3 normal = normalize(cross((p1.xyz - p0.xyz), (p2.xyz - p0.xyz)));
	float delta = width * abs(dot(wind, normal)) * get_wind_delta(p2.xyz, 1.0);
	vec4 offset = vec4(delta*wind_x, delta*wind_y, 0.0, 0.0);
	
	gl_Position = p2 + offset;
	gl_TexCoord[0].st = vec2(0.1, 0.5);
	EmitVertex();
}