uniform float time = 0.0;
uniform float height = 1.0;
uniform float wind_x = 0.0;
uniform float wind_y = 0.0;
uniform sampler2D tex_noise;

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0;
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space (not normalized)
	vec4 color  = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	gl_FrontColor = color;
	vec4 vertex = gl_Vertex;
	
	if (gl_TexCoord[0].s < 0.5) { // top vertex
		// Note: grass motion amplitude should depend on dot(wind, gl_Normal), but the normal is incorrect
		float tx = (time*wind_x + 1.2*vertex.x);
		float ty = (time*wind_y + 1.2*vertex.y);
		float mag = texture2D(tex_noise, vec2(tx, ty)).a;
		float flexibility = 0.75*gl_Color.g + 0.25; // burned grass is less flexible
		float delta = 1.5 * flexibility * height * mag;
		
		// apply x/y delta but maintain the existing height
		vec3 v = normalize(vec3(delta*wind_x, delta*wind_y, height)) * height;
		v.z -= height;
		vertex.xyz += v;
	}
	gl_Position = gl_ModelViewProjectionMatrix * vertex;
} 
