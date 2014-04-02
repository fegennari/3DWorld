uniform float alpha = 1.0;
uniform vec4 emission = vec4(0,0,0,1);
uniform float vertex_offset_scale = 0.0;// hack to make vertex_offset ignored when unused/unset
uniform vec3 xlate = vec3(0);
uniform vec3 scale = vec3(1);

attribute vec3 vertex_offset;

varying vec2 tc;

void main()
{
	tc          = fg_TexCoord;
	vec4 vertex = vec4((vertex_offset_scale*vertex_offset + xlate), 0.0) + (vec4(scale, 1.0) * fg_Vertex);
	vec4 epos   = gl_ModelViewMatrix * vertex;
	gl_Position = gl_ProjectionMatrix * epos;
	vec3 normal = normalize(gl_NormalMatrix * fg_Normal); // eye space
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
	vec4 color  = vec4(emission.rgb, alpha);
	if (enable_light0) color.rgb += add_light_comp0(normal).rgb;
	if (enable_light1) color.rgb += add_light_comp1(normal).rgb;
	gl_FrontColor = color;
}
