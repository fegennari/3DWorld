uniform float fog_scale    = 0.0;
uniform float fog_fade_val = 8.0;

#ifdef USE_NONUNIFORM_FOG
// linear fog from 1.0 at/below fog_bot to 0.0 at/above fog_top
uniform float fog_bot = 0.0;
uniform float fog_top = 1.0;

float get_nonuniform_fog_scale(in float vert_z, in float eye_z, in float fog_bot, in float fog_top, in float density_exp) {
	float zmin   = min(vert_z, eye_z);
	float zmax   = max(vert_z, eye_z);
	float dz_inv = 1.0/(zmax - zmin);
	float tv     = min(fog_top, zmax);
	float bv     = max(fog_bot, zmin);
	float bot_w  = max(0.0, (min(fog_bot, zmax) - zmin)*dz_inv); // constant weight of 1
	float mid_w  = max(0.0, (tv - bv)*dz_inv)*pow((1.0 - (0.5*(tv + bv) - fog_bot)/(fog_top - fog_bot)), density_exp); // use density at midpoint
	return (bot_w + mid_w); // + top_w=0
}

float get_custom_fog_scale(in float vert_z) {
	float eye_z = gl_ModelViewMatrixInverse[3].z; // world space
	return get_nonuniform_fog_scale(vert_z, eye_z, fog_bot, fog_top, 1.0); // linear density
}

//fog *= 1.0 + 10.0*get_nonuniform_fog_scale(vert_z, eye_z, water_plane_z, water_plane_z+2.0, 1.0); // low fog in valleys

#else
float get_custom_fog_scale(in float vert_z) {return 1.0;}
#endif

float get_custom_fog_scale_epos(in vec4 epos) {
	return get_custom_fog_scale((gl_ModelViewMatrixInverse * epos).z);
}

// linear/quadratic fog
vec4 apply_fog_ffc(in vec4 color, in float ffc, in vec4 fog_color) {
	float fog = clamp((gl_Fog.end - ffc) * gl_Fog.scale, 0.0, 1.0);
#ifdef USE_QUADRATIC_FOG
	fog = 1.0 - (1.0-fog)*(1.0-fog); // quadratic term
#endif
	float fog_mix_val = mix(1.0, fog, fog_scale);
	//vec4 fin_color = vec4(mix(fog_color.rgb, color.rgb, fog_mix_val), color.a);
	vec4 fin_color = mix(fog_color, color, fog_mix_val); // fog_mix_val=0 => fog_color, fog_mix_val=1 => color
#ifdef FOG_FADE_TO_TRANSPARENT
	fin_color.a *= min(fog_fade_val*fog_mix_val, 1.0);
#endif
	return fin_color;
}

vec4 apply_fog(in vec4 color) {
	return apply_fog_ffc(color, gl_FogFragCoord, gl_Fog.color);
}

vec4 apply_fog_epos(in vec4 color, in vec4 epos) {
	return apply_fog_ffc(color, length(epos.xyz)*get_custom_fog_scale_epos(epos), gl_Fog.color);
}

vec4 apply_fog_scaled(in vec4 color, in float world_z) {
	return apply_fog_ffc(color, gl_FogFragCoord*get_custom_fog_scale(world_z), gl_Fog.color);
}

vec4 apply_fog_colored(in vec4  color,    // original color of the pixel
                       in float distance, // camera to point distance
                       in vec3  ray_dir,  // camera to point vector
                       in vec3  sun_dir)  // sun light direction
{
#ifdef USE_QUADRATIC_FOG // hack to make this only apply to tiled terrain mode, which uses custom quadratic fog
    float sun_amount = ((sun_dir.z > 0.0) ? max(dot(ray_dir, sun_dir), 0.0) : 0.0);
	vec4 sun_color   = vec4(max(gl_Fog.color.r, gl_Fog.color.b), max(gl_Fog.color.g, gl_Fog.color.b), 0.5*(gl_Fog.color.r + gl_Fog.color.g), gl_Fog.color.a);
    vec4 fog_color   = mix(gl_Fog.color, sun_color, pow(sun_amount, 8.0));
	return apply_fog_ffc(color, distance, fog_color);
#else
	return apply_fog_ffc(color, distance, gl_Fog.color);
#endif
}

vec4 apply_fog_colored(in vec4 color, in vec4 vertex) {
	vec4 epos       = gl_ModelViewMatrix * vertex;
	float fog_scale = gl_FogFragCoord*get_custom_fog_scale(vertex.z);
	return apply_fog_colored(color, fog_scale, -normalize(epos.xyz), normalize(epos.xyz - gl_LightSource[0].position.xyz));
}

vec4 apply_fog_colored_epos(in vec4 color, in vec4 epos) {
	float fog_coord = length(epos.xyz)*get_custom_fog_scale_epos(epos);
	return apply_fog_colored(color, fog_coord, -normalize(epos.xyz), normalize(epos.xyz - gl_LightSource[0].position.xyz));
}

