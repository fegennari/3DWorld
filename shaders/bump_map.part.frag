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
mat3 cotangent_frame(in vec3 N, in vec3 p, in vec2 uv, in float bscale)
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
    float invmax = bump_tb_scale * inversesqrt(max(dot(T,T), dot(B,B)));
    return mat3((T * invmax), (-B * (bscale * invmax)), N);
}
#endif

mat3 get_tbn(in float bscale) {
#ifdef USE_TANGENT_VECTOR
	return transpose(mat3(tangent_v.xyz*tangent_v.w, bscale*cross(eye_norm, tangent_v.xyz), eye_norm));
#else
	// assume N, the interpolated vertex normal and V, the view vector (vertex to eye / camera pos - vertex pos) from VS
    return transpose(cotangent_frame(eye_norm, epos.xyz, tex_coord, bscale));
#endif
}

#ifdef ENABLE_PARALLAX_MAP
uniform sampler2D depth_map;
uniform float hole_depth = 1.0;

vec2 apply_parallax_map() {
    mat3 TBN    = get_tbn(-1.0); // FIXME: why is binormal inverted from bump map case?
	float depth = (texture2D(depth_map, tex_coord).w) * hole_depth; // Get depth from the alpha (w) of the relief map
    return tex_coord + depth * (TBN * -normalize(epos.xyz)).st; // transform view vector to tangent space and offset the uv
}
#endif

// Note: we assume the bump map tex coords are the same as the object diffuse tex coords
vec3 apply_bump_map(inout vec3 light_dir, inout vec3 eye_pos) {
	mat3 TBN  = get_tbn(1.0);
	light_dir = TBN * light_dir;
	eye_pos   = TBN * eye_pos;
	return normalize(texture2D(bump_map, tex_coord).xyz * 2.0 - 1.0);
}
#endif
