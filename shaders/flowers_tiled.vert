uniform vec2 xlate = vec2(0.0);
uniform float x1, y1, dx_inv, dy_inv;
uniform sampler2D height_tex;

out vec4 vertex, epos;
out vec3 eye_norm;
out vec2 tc;

void main()
{
	tc            = fg_TexCoord;
	vec3 gwdelta  = get_grass_wind_delta(fg_Vertex.xyz, 0.5); // wind_amt=0.5
	eye_norm      = fg_NormalMatrix * normalize(normalize(fg_Normal) + gwdelta/height); // eye space, height comes from wind.part
	vertex        = fg_Vertex + vec4(gwdelta, 0.0);
	vertex.z     += texture2D(height_tex, vec2((fg_Vertex.x - x1)*dx_inv, (fg_Vertex.y - y1)*dy_inv)).r;
	epos          = fg_ModelViewMatrix * (vertex + vec4(xlate, 0.0, 0.0));
	gl_Position   = fg_ProjectionMatrix * epos;
	gl_FrontColor = fg_Color;
} 
