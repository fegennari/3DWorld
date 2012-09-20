varying vec4 pos, color;

void main()
{
	gl_FragColor = apply_fog(color*vec4(1,1,1, gen_cloud_alpha(pos.xy)));
}
