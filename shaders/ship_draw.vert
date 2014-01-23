// for additional texture scaling control
uniform float tscale_s = 1.0;
uniform float tscale_t = 1.0;

varying vec4 epos;
varying vec3 normal; // eye space

void main()
{
	gl_TexCoord[0] = gl_MultiTexCoord0 + vec4(tscale_s, tscale_t, 0.0, 0.0);
	normal = normalize(gl_NormalMatrix * gl_Normal);
	epos   = gl_ModelViewMatrix * gl_Vertex;
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
}
