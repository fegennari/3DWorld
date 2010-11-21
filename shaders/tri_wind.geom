#version 120
#extension GL_EXT_geometry_shader4 : enable

uniform float time = 0.0;
uniform float wind_x = 0.0;
uniform float wind_y = 0.0;
uniform sampler2D tex_noise;

void main()
{
	vec3 wind = normalize(vec3(wind_x, wind_y, 0.0));

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
	float width = length(p1 - p0);
	//float mag = sin(time);
	float tx = (0.2*time + 1.0*p2.x)*wind_x;
	float ty = (0.2*time + 1.0*p2.y)*wind_y;
	float mag = 1.0 - texture2D(tex_noise, vec2(tx, ty)).a;
	vec3 normal = normalize(cross((p1.xyz - p0.xyz), (p2.xyz - p0.xyz)));
	float delta = 10.0*width*mag * abs(dot(wind, normal));
	vec4 offset = vec4(delta*wind_x, delta*wind_y, 0.0, 0.0);
	
	gl_Position = p2 + offset;
	gl_TexCoord[0].st = vec2(0.1, 0.5);
	EmitVertex();
}