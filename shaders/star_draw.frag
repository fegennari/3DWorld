uniform float radius = 1.0;
uniform sampler2D tex0;
uniform vec4 colorA, colorB;

in vec3 world_space_pos;
in vec2 tc;

void main()
{
	vec3 spos    = world_space_pos*radius/sqrt(length(world_space_pos));
	fg_FragColor = texture(tex0, tc) * mix(colorA, colorB, gen_cloud_alpha(spos));
}
