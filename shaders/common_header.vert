layout(location = 0) in vec4 fg_Vertex;
layout(location = 1) in vec3 fg_Normal;
layout(location = 2) in vec4 fg_Color;
layout(location = 3) in vec2 fg_TexCoord;

vec4 fg_ftransform() {return gl_ModelViewProjectionMatrix * fg_Vertex;}
//vec4 fg_ftransform() {return fg_ProjectionMatrix * fg_ModelViewMatrix * fg_Vertex;}

#define gl_Vertex fg_Vertex
#define gl_Normal fg_Normal
#define gl_Color  fg_Color
#define gl_MultiTexCoord0 fg_TexCoord
#define ftransform fg_ftransform
