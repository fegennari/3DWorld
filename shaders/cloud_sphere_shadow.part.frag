uniform float atmosphere    = 1.0;
uniform float sphere_radius = 1.0;
uniform vec3 sphere_center  = vec3(0.0);
uniform float dx, dy;
uniform vec3 sun_pos, moon_pos;
uniform sampler2D cloud_tex;

float get_cloud_transparency(in vec3 pos, in vec3 lpos) {
	vec3 v1 = normalize(pos - lpos);
	vec3 v2 = pos - sphere_center; // dir from target sphere to object
	v2 -= v1*dot(v1, v2);
	vec3 vdir = v2 - v1*sqrt(max((sphere_radius*sphere_radius - dot(v2, v2)), 0.0));
	float val = texture(cloud_tex, vec2((vdir.x/sphere_radius + dx), (vdir.y/sphere_radius + dy))).a;
	return mix(1.0, (1.0 - clamp(1.7*val, 0.0, 0.99)), atmosphere);
}