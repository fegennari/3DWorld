uniform float radius_scale = 1.0;

varying float world_space_zval;

void main()
{
	set_tc0_from_vert_id();
	vec4 eye    = gl_ModelViewMatrixInverse[3]; // world space
	vec3 dir    = normalize(vec3(-(fg_Vertex.y - eye.y), (fg_Vertex.x - eye.x), 0.0)); // cross(z, fg_Vertex-eye)
	vec4 vertex = fg_Vertex + vec4(((2.0*tc.s - 1.0) * radius_scale * fg_Normal.x * dir), 0.0);
	world_space_zval = vertex.z;
	gl_Position = gl_ModelViewProjectionMatrix * vertex;
	gl_FogFragCoord = length((gl_ModelViewMatrix * vertex).xyz); // set standard fog coord
	//const float normal_z = 0.816; // high detail tree polygon normal
	const float normal_z = 0.95; // blends better in practice
	vec3 normal = normalize(normal_z*gl_NormalMatrix[2] + vec3(0.0, 0.0, sqrt(1.0 - normal_z*normal_z))) / ambient_scale;
	vec3 color  = vec3(0.0);
	if (enable_light0) color += add_light_comp0(normal).rgb;
	if (enable_light1) color += add_light_comp1(normal).rgb;
	gl_FrontColor = vec4(color, fg_Color.a);
} 
