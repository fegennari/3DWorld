uniform sampler2D tex0;
in vec2 tc;
in vec4 epos;

void main()
{
	if (gl_Color.a == 0.0) discard;
	fg_FragColor = apply_fog(gl_Color* texture(tex0, tc)); // Note: no custom fog
#ifndef NO_FOG
	fg_FragColor = apply_fog_epos(fg_FragColor, epos);
#endif
}
