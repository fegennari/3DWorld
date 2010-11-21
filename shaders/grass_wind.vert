uniform float time = 0.0;
uniform float dist = 1.0;
uniform float wind_x = 0.0;
uniform float wind_y = 0.0;
uniform sampler2D tex_noise;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space (not normalized)
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	gl_FrontColor = color;
	vec4 vertex = gl_Vertex;
	
	if (gl_TexCoord[0].s < 0.5) { // top vertex
		vec3 wind = normalize(vec3(wind_x, wind_y, 0.0));
		//float mag = sin(time);
		float tx = (time + 1.2*vertex.x)*wind_x;
		float ty = (time + 1.2*vertex.y)*wind_y;
		float mag = texture2D(tex_noise, vec2(tx, ty)).a;
		float delta = dist * mag * abs(dot(wind, normalize(normal)));
		vertex += vec4(delta*wind_x, delta*wind_y, 0.0, 0.0);
	}
	gl_Position = gl_ModelViewProjectionMatrix * vertex;
} 
