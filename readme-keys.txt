3D World
Frank Gennari (cs184-hl)
Hiral Patel (cs184-bg)

Options:

-uel	enable record user mouse and keyboard events (save with key e)
<file>	play back recorded user event file named "file"



Mouse:

Air mode (God mode):

LEFT BUTTON		latitude(x)/longitude(y) rotate around camera origin
RIGHT BUTTON	zoom in/out (y), rotate up vector (x), fire (game mode only)
MIDDLE BUTTON   x/y translate of camera origin

Ground Mode (human mode):

LEFT BUTTON		change camera orientation vector



Key Bindings:

a	step left
b	toggle object enable/disable (rain, snow, hail, ball, etc., default = OFF)
c	toggle show smiley (default = ON)
d	step right
e	next weapon (default = ball)
f	print framerate to console along with other statistics
g	pause/resume playback of a user eventlist
h	toggle camera collision detection in ground mode (default = OFF)
i	switch between island, raised island, and land floating in space (must recreate mesh with c afterwards, default = island)
j	toggle camera real physics/collision (default = OFF)
k	switch between filled polygon and wireframe mesh in normal mode (default = filled), switch between color and grayscale in mapmode (default = COLOR)
l	toggle lightning OFF, ON, NO_SHADOWS (default = OFF)
m	toggle display fullscreen/windowed with a min time delay of 5 seconds between switch
n	toggle fog (default = OFF)
o	increase shadow detail
p	reset camera position and orientation
q	pervious weapon (default = ball)
r	redraw screen or frame advance (with animation disabled)
s	move camera backwards if on ground
t	freeze frame for objects, toggle stretched stars in universe mode (default = ON)
u	decrease shadow detail
v	toggle camera air/ground mode (default = air)
w	move camera forward if on ground
x	toggle object animation (default = ON)
y   reocord user event list for later playback in file "ueventlist"
z	zoom in 5X

A	decrease leaf red color component by 0.1
B	toggle show precipitation (default = ON)

D	take screenshot in .raw format
E	reset leaf color properties
F	take screenshot in .jpg (JPEG) format
G	toggle on-screen framerate display (default = OFF)
H	save mesh state
I	save mesh values
J	load mesh state
K	toggle overhead map mode (default = OFF)
L	increase terrain zoom
M	decrease precipitation rate by 1.5
N	increase precipitation rate by 1.5
O	decrease tree color coherence
P	increase tree color coherence
Q	decrease leaf color coherence by 0.5
R	toggle run mode 3 levels (default = OFF)
S	increase leaf red color component by 0.1
T	create/recreate tree
U	universe: toggle universe mode PLANET vs. UNIVERSE (default = UNIVERSE), landscape: toggle allow visibility to be calculated in another thread (default = ON)
V	toggle dynamic timestep (default = ON)
W	increase leaf color coherence by 0.5
X	increase leaf green color component by 0.1
Y	decrease terrain zoom
Z	decrease leaf green color component by 0.1, zoom
C	create new mesh
1	toggle mesh/rock display (default = mesh ON)
2	toggle tree display (default = OFF)
3	toggle water/ice effects/display (default = OFF)
4	toggle display peaks, ridges, valleys, etc. (slow, default = OFF)
5	toggle display watersheds (default = OFF)
6	switch between mesh lighting/shadow models (experimental, default is fastest model)
7	toggle show snow accumulation (slower, default = OFF)
8	toggle mouse draw mode (must enable snow with key 7, default = OFF)
9	toggle ocean waves (default = ON)
0	toggle stencil shadows off mesh/on water (default = ON)

-	decrease temperature by 10 degrees C (starts at 20 degrees C)
=	decrease temperature by 10 degrees C
[	increase sun rotation angle
]	decrease sun rotation angle
{	increase moon rotation angle
}	decrease moon rotation angle
;	(semicolon)	decrease simulation physics timestep/speed by 1.5
'	(quote)	increase simulation physics timestep/speed by 1.5
,	(comma) decrease mesh resolution by 2 (default = 200x200)
.	(period) increase mesh resolution by 2 (up to 200x200)
\	toggle show camera sphere (default = OFF)
<	decrease ball throw velocity (default is 0.0, drop)
>	increase ball throw velocity

LEFT_ARROW		decrease wind in x direction by 0.5
RIGHT_ARROW		increase wind in x direction by 0.5
UP_ARROW		decrease wind in y direction by 0.5
DOWN_ARROW		increase wind in y direction by 0.5

SPACE			fire (game mode only)
TAB				show scores (game mode only)

F1	switch mode: universe/planet, experimental land, infinite terrain, island/land (default = universe)
F2	toggle game mode/smileys (default = ON)
F3	switch world display between normal and "experimental" (default = normal)
F4  switch weapon firing mode
F5  toggle large/small trees (default = small)
F6  change gameplay moving/firing mode (3 modes)
F7  toggle auto day time advance: sun/moon position, precipitation, temperature, cloud cover (default = OFF)
F8	toggle spectator gameplay mode

<esc>	quit
