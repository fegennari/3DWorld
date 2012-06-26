uniform float camera_facing_scale = 0.0;

void main()
{
	gl_TexCoord[0]  = gl_MultiTexCoord0;
	gl_Position     = ftransform();
	vec3 epos       = (gl_ModelViewMatrix * gl_Vertex).xyz;
	gl_FogFragCoord = length(epos); // set standard fog coord
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal);
	vec3 nz     = (gl_NormalMatrix * vec3(0.0, 0.0, gl_Normal.z));
	vec3 n_mod  = normalize(nz + vec3(0.0, 0.0, sqrt(1.0 - nz.x*nz.x - nz.y*nz.y - nz.z*nz.z))) * camera_facing_scale;
	n_mod      += normal * ((two_sided_lighting ? normal.z/abs(normal.z) : 1.0) * (1.0 - camera_facing_scale));
	vec3 color  = vec3(0,0,0);
	if (enable_light0) color += add_light_comp(n_mod, 0).rgb;
	if (enable_light1) color += add_light_comp(n_mod, 1).rgb;
	gl_FrontColor = vec4(color, 1.0);
} 
