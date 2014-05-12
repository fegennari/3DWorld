varying vec4 epos;

void main()
{
	fg_FragColor = apply_fog_epos(gl_Color, epos);
}
