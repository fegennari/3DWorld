uniform sampler2D tex0;
uniform vec4 colorA, colorB;

varying vec3 world_space_pos;
varying vec2 tc;

void main()
{
	fg_FragColor = texture2D(tex0, tc) * mix(colorA, colorB, gen_cloud_alpha(world_space_pos));
}
