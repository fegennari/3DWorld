in vec2 tc;

uniform sampler2D msdf;
uniform vec4 bg_color;
uniform vec4 fg_color;
uniform float px_range; // set to distance field's pixel range

float screen_px_range() {
    vec2 unit_range = vec2(px_range)/vec2(textureSize(msdf, 0));
    vec2 screen_tex_size = vec2(1.0)/fwidth(tc);
    return max(0.5*dot(unit_range, screen_tex_size), 1.0);
}

float median(float r, float g, float b) {
    return max(min(r, g), min(max(r, g), b));
}

void main() {
    vec3 msd = texture(msdf, tc).rgb;
    float sd = median(msd.r, msd.g, msd.b);
    float screen_px_distance = screen_px_range()*(sd - 0.5);
    float opacity = clamp(screen_px_distance + 0.5, 0.0, 1.0);
    fg_FragColor  = mix(bg_color, fg_color, opacity);
}

