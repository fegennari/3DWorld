uniform sampler2D tex0;
uniform float min_alpha = 0.0;

varying vec3 normal;
varying vec2 tc;

void main()
{
	if (texture2D(tex0, tc).a <= min_alpha) discard; // transparent
	fg_FragColor = vec4(0.5*(normal + 1.0), 1.0); // Note: normal not normalized
}
