// for additional texture scaling control
uniform float tscale_s = 1.0;
uniform float tscale_t = 1.0;

varying vec2 tc;
varying vec4 epos;
varying vec3 normal; // eye space

void main()
{
	tc = fg_TexCoord*vec2(tscale_s, tscale_t);
	normal = normalize(gl_NormalMatrix * fg_Normal);
	epos   = gl_ModelViewMatrix * fg_Vertex;
	gl_Position   = fg_ftransform();
	gl_FrontColor = fg_Color;
}
