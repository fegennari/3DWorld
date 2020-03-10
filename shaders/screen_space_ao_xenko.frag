// https://old.reddit.com/r/opengl/comments/96api8/has_anyone_successfully_implemented_groundtruth_ao/
// https://github.com/xenko3d/xenko/blob/master/sources/engine/Xenko.Rendering/Rendering/Images/AmbientOcclusion/AmbientOcclusionRawAOShader.xksl

in vec2 tc;

// Copyright (c) Xenko contributors (https://xenko.com) and Silicon Studio Corp. (https://www.siliconstudio.co.jp)
// Distributed under the MIT license. See the LICENSE.md file in the project root for more information.

uniform float    ParamProjScale = 1;
uniform float    ParamIntensity = 1;
uniform float    ParamBias = 0.01f;
uniform float    ParamRadius = 1;
uniform float    ParamRadiusSquared = 1;

uniform int SamplesCount = 4;

//float reconstructCSZ(float depth) {return ZProjection.y / (depth - ZProjection.x);}

vec3 reconstructCSPosition(vec2 S, float z) {
    vec4 ProjInfo = vec4(znear*zfar, znear-zfar, zfar, 1.0); // .x = zN * zF, .y = zN - zF, .z = zF
    return vec3((S.xy * ProjInfo.xy + ProjInfo.zw) * z, z);
}

vec3 reconstructCSNormal(vec3 position) {
    return normalize(cross(dFdy(position), dFdx(position)));
}

float sampleAO(ivec2 screenPosition, vec3 viewPosition, vec3 viewNormal, float diskRadius, int i, float randomPatternRotationAngle) {
    //*****************************
    //  Sample Offset
    float alpha = 1 * (i + 0.5) * 0.675f / SamplesCount;
    float angle = 1 * 43.9822971503f * alpha + randomPatternRotationAngle;

    vec2 offset = vec2(cos(angle), sin(angle));
    float ssRadius = alpha * diskRadius;

    //*****************************
    //  Depth
    vec2 samplePos = tc + offset * ssRadius;
    //ivec2 samplePosInt = clamp((samplePos) / xy_step, 0.0, 1.0);     
    //float depth = Texture0.Load(ivec3(samplePosInt, 0));
    //float linearDepth = reconstructCSZ(depth);
    float linearDepth = get_linear_depth_zval(samplePos);

    //*****************************
    // View Position
    vec3 position = reconstructCSPosition(samplePos + vec2(0.5, 0.5), linearDepth);
    position.x *= -1;

    //*****************************
    // View Normal
    vec3 v = position - viewPosition;
    v.z *= -1;
            
    //*****************************
    // Ambient Occlusion
    float distSq = dot(v, v);
    float vn = dot(v, viewNormal);

    const float epsilon = 0.01;

    float f = max(ParamRadiusSquared - distSq, 0.0);

    return f * f * f * max((vn - ParamBias) / (epsilon + distSq), 0.0);
}

void main() {
    //*****************************
    // Reconstruct View space linear depth Z from the depth buffer
    //float depth = Texture0.SampleLevel(Sampler, tc, 0).x;
    //float linearDepth = reconstructCSZ(depth);
    float linearDepth = get_linear_depth_zval(tc);

    //*****************************
    // Reconstruct View space position XYZ
    vec2 screenPosition = tc * xy_step;
    vec3 viewPosition = reconstructCSPosition(screenPosition + vec2(0.5, 0.5), linearDepth);
    viewPosition.x *= -1;

    //*****************************
    // Reconstruct View space normal NxNyNz
    vec3 viewNormal = reconstructCSNormal(viewPosition.xyz);
    viewNormal.xy *= -1;

    //*****************************
    // Hash function used in the HPG12 AlchemyAO paper
    int linearDepthInt = int(linearDepth);
    //float randomPatternRotationAngle = (3 * screenPosition.x ^ screenPosition.y + screenPosition.x * screenPosition.y) * 10;
    float randomPatternRotationAngle = (15 * linearDepthInt + 3 * screenPosition.x + 2 * screenPosition.y + screenPosition.x * screenPosition.y) * 10;

    //*****************************
    // Choose a sample radius proportional to the projected area of the half-sphere
    //float diskRadius = -projScale * radius / linearDepth;
    float diskRadius = ParamProjScale / linearDepth;

    //*****************************
    // Compute the ambient occlusion
    float sum = 0.0;
    for (int i = 0; i < SamplesCount; i++) {
        sum += sampleAO(ivec2(screenPosition), viewPosition, viewNormal, diskRadius, i, randomPatternRotationAngle);
    }

    float temp = ParamRadiusSquared * ParamRadius;
    sum /= temp * temp;
	float A = max(0.0, 1.0 - sum * 5 * ParamIntensity / SamplesCount);
            
    float nearPlaneFade = clamp(linearDepth * 2.0 - 0.5, 0.0, 1.0);
    A = mix(1, A, nearPlaneFade);

    //*****************************
    // Bilateral box-filter over a quad for free, respecting depth edges
	// (the difference that this makes is subtle)
    if (abs(dFdx(linearDepth)) < 0.02) {
		A -= dFdx(A) * ((int(screenPosition.x) & 1) - 0.5);
	}
	if (abs(dFdy(linearDepth)) < 0.02) {
        A -= dFdy(A) * ((int(screenPosition.y) & 1) - 0.5);
	}

    //************************
    // A now contains the light intensity factor (0 to 1) which can be applied to the ambient light illuminating the pixel

    fg_FragColor = vec4(0.0, 0.0, 0.0, 1.0-A);
}
