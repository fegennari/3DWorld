uniform float water_plane_z = 0.0;
uniform float water_atten = 1.0;

// underwater attenuation code
vec4 atten_color(in vec4 color, in float dist) {
	color.r *= (1.0 - min(0.98, 1.5*dist));
	color.g *= (1.0 - min(0.97, 0.9*dist));
	color.b *= (1.0 - min(0.95, 0.5*dist));
	return color;
}

float integrate_water_dist(in vec3 targ_pos, in vec3 src_pos, in float water_z) {
	float t = clamp((water_z - targ_pos.z)/max(1.0E-6, abs(src_pos.z - targ_pos.z)), 0.0, 1.0);
	return t*length(src_pos - targ_pos);
}

vec4 add_light_comp(in vec3 normal, in int i) {
	// normalize the light's direction in eye space, directional light: position field is actually direction
	vec3 lightDir = normalize(gl_LightSource[i].position.xyz);
		
	// compute the cos of the angle between the normal and lights direction as a dot product, constant for every vertex.
	float NdotL = dot(normal, lightDir);
	
	// compute the ambient and diffuse lighting
	vec4 diffuse = gl_Color * gl_LightSource[i].diffuse;
	vec4 ambient = gl_Color * gl_LightSource[i].ambient;
	vec4 color = (ambient + max(dot(normal, lightDir), 0.0)*diffuse);
	
	// apply underwater attenuation
	if (gl_Vertex.z < water_plane_z) {
		//float dist = 2.5*(water_plane_z - gl_Vertex.z); // depth
		mat4 mm_inv = inverse(gl_ModelViewMatrix);
		vec4 eye = mm_inv * vec4(0.0, 0.0, 0.0, 1.0); // world space
		vec4 light = mm_inv * gl_LightSource[i].position; // world space
		float dist = integrate_water_dist(gl_Vertex.xyz, eye.xyz, water_plane_z) + integrate_water_dist(gl_Vertex.xyz, light.xyz, water_plane_z);
		color = atten_color(color, dist*water_atten);
	}
	return color;
}


// main code
void main()
{
	setup_texgen(0);
	setup_texgen(1);
	gl_Position = ftransform();
	vec3 normal = gl_NormalMatrix * gl_Normal; // eye space, not normalized
	vec4 color = gl_Color * gl_LightModel.ambient;
	if (enable_light0) color += add_light_comp(normal, 0);
	if (enable_light1) color += add_light_comp(normal, 1);
	gl_FrontColor = color;
	set_fog();
} 
