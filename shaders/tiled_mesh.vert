uniform float water_plane_z = 0.0;
uniform float water_atten = 1.0;

// underwater attenuation code
void atten_color(inout vec4 color, in float dist) {
	color.rgb *= vec3(1,1,1) - min(vec3(0.98, 0.97, 0.95), vec3(1.5, 0.9, 0.5)*dist);
}

float integrate_water_dist(in vec3 targ_pos, in vec3 src_pos, in float water_z) {
	float t = clamp((water_z - targ_pos.z)/max(1.0E-6, abs(src_pos.z - targ_pos.z)), 0.0, 1.0);
	return t*distance(src_pos, targ_pos);
}

vec4 add_light_comp(in vec3 normal, in int i) {
	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(gl_LightSource[i].position.xyz);
		
	// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex.
	float NdotL = dot(normal, lightDir);
	
	// compute the ambient and diffuse lighting
	vec4 diffuse = gl_LightSource[i].diffuse;
	vec4 ambient = gl_LightSource[i].ambient;
	vec4 color   = (ambient + max(dot(normal, lightDir), 0.0)*diffuse);
	
	// apply underwater attenuation
	// Note: ok if vertex is above the water, dist will come out as 0
	vec4 eye   = gl_ModelViewMatrixInverse[3]; // world space
	vec4 light = gl_ModelViewMatrixInverse * gl_LightSource[i].position; // world space
	float dist = integrate_water_dist(gl_Vertex.xyz, eye.xyz, water_plane_z) + integrate_water_dist(gl_Vertex.xyz, light.xyz, water_plane_z);
	atten_color(color, dist*water_atten);
	return color;
}


// main code
void main()
{
	setup_texgen(0);
	setup_texgen(1);
	gl_Position = ftransform();
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space, not normalized
	vec4 color  = gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	gl_FrontColor = color;

	// calculate fog coord
#if 0
	gl_FogFragCoord = length((gl_ModelViewMatrix * gl_Vertex).xyz);
#else
	// clip the line to the water plane if the eye is above the water
	// could use different terms/fog scaling/color for inside/outside water?
	vec4 eye = gl_ModelViewMatrixInverse[3]; // world space
	float t = min(1.0, (eye.z - water_plane_z)/max(0.0, (eye.z - gl_Vertex.z)));
	vec4 clipped_vert = mix(eye, gl_Vertex, ((t < 0.0) ? 1.0 : t));
	gl_FogFragCoord = length((gl_ModelViewMatrix * clipped_vert).xyz);
#endif
} 
