uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta;
attribute float shadowed;
varying vec3 eye, vpos, spos;

#define ADD_LIGHT(en, N) if (en && (s&(1<<N)) == 0) color.rgb += add_light_comp(normal, N).rgb;

void main()
{
	if (use_texgen) {
		setup_texgen(0);
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}	
	gl_Position = ftransform();
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	vec4 color  = gl_Color;
	int s = int(round(shadowed));
	ADD_LIGHT(enable_light0, 0);
	ADD_LIGHT(enable_light1, 1);
	gl_FrontColor = clamp(color, 0.0, 1.0);
	
	spos = gl_Vertex.xyz + (0.25*step_delta)*normalize(gl_Normal); // move slightly away from the vertex
	
	if (!smoke_enabled) {
		eye  = vec3(0,0,0);
		vpos = vec3(0,0,0);
		return;
	}
	vec3 v2 = (inverse(gl_ModelViewMatrix) * vec4(0.0, 0.0, 0.0, 1.0)).xyz; // world space
	pt_pair res = clip_line(gl_Vertex.xyz, v2, smoke_bb);
	eye  = res.v1;
	vpos = res.v2;
} 
