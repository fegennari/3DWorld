uniform float alpha = 1.0;
uniform vec4 emission = vec4(0,0,0,1);
uniform float vertex_offset_scale = 0.0;// hack to make vertex_offset ignored when unused/unset
attribute vec3 vertex_offset;

varying vec2 tc;

void main()
{
	tc          = gl_MultiTexCoord0;
	vec4 epos   = gl_ModelViewMatrix * (gl_Vertex + vec4(vertex_offset_scale*vertex_offset, 0.0));
	gl_Position = gl_ProjectionMatrix * epos;
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
	vec4 color  = vec4(emission.rgb, alpha);
	if (enable_light0) color.rgb += add_light_comp0(normal).rgb;
	if (enable_light1) color.rgb += add_light_comp1(normal).rgb;
	gl_FrontColor = color;
}
