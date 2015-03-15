uniform float alpha = 1.0;
uniform vec4 emission = vec4(0,0,0,1);
uniform float vertex_offset_scale = 0.0;// hack to make vertex_offset ignored when unused/unset
uniform vec3 xlate = vec3(0);
uniform vec3 scale = vec3(1);

in vec3 vertex_offset;

out vec2 tc;

void main()
{
	tc          = fg_TexCoord;
	vec4 vertex = vec4((vertex_offset_scale*vertex_offset + xlate), 0.0) + (vec4(scale, 1.0) * fg_Vertex);
	vec4 epos   = fg_ModelViewMatrix * vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	vec3 normal = normalize(fg_NormalMatrix * fg_Normal); // eye space
	gl_FogFragCoord = length(epos.xyz); // set standard fog coord
	vec3 color  = emission.rgb;
	if (enable_light0) color += add_light_comp_pos0(normal, epos).rgb;
	if (enable_light1) color += add_light_comp_pos1(normal, epos).rgb;
	fg_Color_vf = vec4(color, alpha);
}
