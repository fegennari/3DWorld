uniform float water_plane_z;
uniform vec3 cloud_offset = vec3(0,0,0);
uniform vec4 sun_color;
uniform vec3 sun_pos, camera_pos;

varying vec3 vertex;
varying vec4 color;

void main()
{
	float dp = max(0.0, dot(normalize(vertex - camera_pos), normalize(sun_pos - vertex)));
	vec4 color2 = 0.75*color + vec4(0.5*sun_color.rgb*pow(dp, 8.0), 0.0); // add sun glow
	vec3 pos = vertex + cloud_offset;
	float eye_z = gl_ModelViewMatrixInverse[3].z; // world space
	float t = min(1.0, (eye_z - water_plane_z)/max(0.0, (eye_z - pos.z)));
	float black_mix = ((underwater_atten && t > 0.0 && t < 1.0) ? 1.0 : 0.0);
	float alpha     = gen_cloud_alpha(pos.xy);
	gl_FragColor    = apply_fog(color2*mix(vec4(1,1,1, alpha), vec4(0,0,0,1), black_mix));
}
