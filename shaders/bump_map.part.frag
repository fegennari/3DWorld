varying vec4 epos;
varying vec3 eye_norm;
varying vec2 tex_coord;

#ifdef USE_BUMP_MAP
uniform sampler2D bump_map;

#ifdef USE_TANGENT_VECTOR
varying vec4 tangent_v;

#else

uniform float bump_tb_scale = 1.0;

// http://www.thetenthplanet.de/archives/1180
mat3 cotangent_frame(in vec3 N, in vec3 p, in vec2 uv)
{
    // get edge vectors of the pixel triangle
    vec3 dp1  = dFdx(p);
    vec3 dp2  = dFdy(p);
    vec2 duv1 = dFdx(uv);
    vec2 duv2 = dFdy(uv);
 
    // solve the linear system
    vec3 dp2perp = cross(dp2, N);
    vec3 dp1perp = cross(N, dp1);
    vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
    vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;
 
    // construct a scale-invariant frame 
    float invmax = inversesqrt(max(dot(T,T), dot(B,B)));
    return mat3((bump_tb_scale * T * invmax), (bump_tb_scale * -B * invmax), N);
}
#endif

// Note: we assume the bump map tex coords are the same as the object diffuse tex coords
vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos) {
#ifdef USE_TANGENT_VECTOR
	mat3 TBN  = transpose(mat3(tangent_v.xyz*tangent_v.w, cross(eye_norm, tangent_v.xyz), eye_norm));
#else
	// assume N, the interpolated vertex normal and V, the view vector (vertex to eye / camera pos - vertex pos) from VS
    mat3 TBN  = transpose(cotangent_frame(eye_norm, epos.xyz, tex_coord));
#endif
	light_dir = TBN * light_dir;
	eye_pos   = TBN * eye_pos;
	return normalize(texture2D(bump_map, tex_coord).xyz * 2.0 - 1.0);
}
#endif
