uniform vec3 planet_pos, sun_pos;
uniform float planet_radius, ring_ri, ring_ro, sun_radius;
uniform sampler2D noise_tex;
uniform sampler1D ring_tex;
varying vec4 epos;
varying vec3 normal, world_space_pos, vertex;

float pt_line_dist(in vec3 P, in vec3 L1, in vec3 L2) {
	return length(cross((L2 - L1), (L1 - P)))/length(L2 - L1);
}

vec4 add_light_rings(in vec3 n)
{
	vec3 light_dir = normalize(gl_LightSource[0].position.xyz - epos.xyz); // normalize the light's direction in eye space
	vec4 diffuse   = gl_Color * gl_LightSource[0].diffuse;
	vec4 ambient   = gl_Color * gl_LightSource[0].ambient;
	vec4 specular  = get_light_specular(n, light_dir, epos.xyz, 0);
	float atten    = calc_light_atten(epos, 0);
	
	if (sun_radius > 0.0) { // sun exists
		float dist_to_sun = length(sun_pos - world_space_pos);

		if (dist_to_sun > length(sun_pos - planet_pos)) { // behind the sun
			const float PI = 3.14159;
			float d = pt_line_dist(planet_pos, sun_pos, world_space_pos);
			float r = planet_radius;
			float R = sun_radius*length(planet_pos - world_space_pos)/dist_to_sun;
			float tot_area = PI*R*R;

			if (d < abs(R - r)) { // fully overlapped
				atten *= 1.0 - PI*min(r,R)*min(r,R)/tot_area;
			}
			else if (d < (r + R)) { // partially overlapped
				float shadowed_area = r*r*acos((d*d+r*r-R*R)/(2*d*r)) + R*R*acos((d*d+R*R-r*r)/(2*d*R)) - 0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
				atten *= 1.0 - clamp(shadowed_area/tot_area, 0.0, 1.0);
			}
		}
	}
	return (ambient + (abs(dot(n, light_dir))*diffuse + specular)) * atten;
}

void main()
{
	float rval = clamp((length(vertex) - ring_ri)/(ring_ro - ring_ri), 0.0, 1.0);
	vec4 texel = texture1D(ring_tex, rval);
	if (texel.a == 0.0) discard;

	vec2 tc = 16*gl_TexCoord[0].st;
	vec3 norm2 = normalize(normal + vec3(texture2D(noise_tex, tc).r-0.5, texture2D(noise_tex, tc+vec2(0.4,0.7)).r-0.5, texture2D(noise_tex, tc+vec2(0.3,0.8)).r-0.5));
	vec4 color = gl_FrontMaterial.emission;
	color.rgb += add_light_rings(norm2).rgb; // ambient, diffuse, and specular
	color.rgb += (gl_Color * gl_LightSource[1].ambient).rgb; // ambient only
	color.a   *= texture2D(noise_tex, 25*gl_TexCoord[0].st).r;
	gl_FragColor = apply_fog(color*texel);
}
