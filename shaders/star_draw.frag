uniform sampler2D tex0;
uniform vec4 colorA, colorB;

varying vec3 world_space_pos;
varying vec2 tc;

void main()
{
	vec4 texel   = texture2D(tex0, tc);
	gl_FragColor = texel * mix(colorA, colorB, gen_cloud_alpha(world_space_pos));
}
