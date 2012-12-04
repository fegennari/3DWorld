uniform sampler2D tex0;
uniform vec4 colorA, colorB;
varying vec3 world_space_pos;

void main()
{
	vec4 texel   = texture2D(tex0, gl_TexCoord[0].st);
	gl_FragColor = texel * mix(colorA, colorB, gen_cloud_alpha(world_space_pos));
}
