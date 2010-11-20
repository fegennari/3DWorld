uniform float fog_scale = 0.0;
varying float fogFactor;

void set_fog()
{
	float z = length(gl_ModelViewMatrix * gl_Vertex);
	//fogFactor = exp2(-1.442695 * fog_scale * gl_Fog.density * z);
	fogFactor = mix(1.0, (gl_Fog.end - z)/(gl_Fog.end - gl_Fog.start), fog_scale);
	fogFactor = clamp(fogFactor, 0.0, 1.0);
}
