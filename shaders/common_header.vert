layout(location = 0) in vec4 fg_Vertex;
layout(location = 1) in vec3 fg_Normal;
layout(location = 2) in vec4 fg_Color;
layout(location = 3) in vec2 fg_TexCoord;

vec4 fg_ftransform() {return fg_ModelViewProjectionMatrix * fg_Vertex;}

// required because some shared shader parts can be used with either vertex or fragment shaders
#define gl_Color fg_Color

out vec4 fg_Color_vf;
#define gl_FrontColor fg_Color_vf
