uniform float water_plane_z;
varying vec4 pos, color;

void main()
{
	float eye_z = gl_ModelViewMatrixInverse[3].z; // world space
	float t = min(1.0, (eye_z - water_plane_z)/max(0.0, (eye_z - pos.z)));
	float black_mix = ((underwater_atten && t > 0.0 && t < 1.0) ? 1.0 : 0.0);
	float alpha     = gen_cloud_alpha(pos.xy);
	gl_FragColor    = apply_fog(color*mix(vec4(1,1,1, alpha), vec4(0,0,0,1), black_mix));
}
