uniform sampler2D tex0;
varying vec3 eye, vpos;

// add: volume texture, scene bounds, step delta

void main()
{
	vec4 texel = texture2D(tex0, gl_TexCoord[0].st);
	vec4 color = vec4(texel.rgb * gl_Color.rgb, texel.a * gl_Color.a);
	
	float density = 0.0;
	float scolor  = 0.0;
	vec3 pos = vpos;
	vec3 dir = eye - pos;
	
	// write smoke volume iteration using 3D texture, pos to eye
	
	float brightness = density*scolor;
	color.a = (1.0 - density)*color.a + density; // deal with transparent objects
	if (scolor < 1.0) color *= (1.0 - density)/(1.0 - brightness); // BLACK: attenuation
	float fog_coord = 15.0*brightness; // WHITE(GRAY): fog coord
	float fog = clamp((gl_Fog.end - fog_coord) * gl_Fog.scale, 0.0, 1.0);
	gl_FragColor = mix(gl_Fog.color, color, fog);
}
