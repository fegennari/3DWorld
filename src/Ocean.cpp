// 3D World - OpenGL CS184 Computer Graphics Project
// by Hiral Patel
// 3/10/02

#include "3DWorld.h"
#include "mesh.h"
#include "inlines.h"


float const OCEAN_REPEAT  = 0.2;
float const SAND_REPEAT   = 0.5;
float const SAND_DEPTH    = 0.5;
float const MAX_VIS_DEPTH = 0.1;
float const MWAVE_HEIGHT  = 0.007;
float const MWAVE_PERIOD  = 0.05;
float const OCEAN_SKEW    = 0.99;

int const SIMPLE_OCEAN  = 0;
int const MOVING_OCEAN  = 1;
int const FORCE_ALPHA_1 = 0; // dynamic ocean
int const FIRST_OCEAN   = 1;
int const LAST_OCEAN    = 2;


struct WAVE {

   double direction, dirX, dirY;
   double lambda;    // wavelength
   double k;         // wave number: k=lambda/(2*Pi)
   double omega;     // angular frequency: omega=2*Pi*f
   double freq;      // frequency
   double periode;   // periode of the wave
   double amplitude; // amplitude of the wave
   double phase;     // phase of the wave
};


// extern global variables
extern int island, camera_mode, draw_model, display_mode, show_fog, frame_counter, DISABLE_WATER;
extern point ocean;
extern float zmin, light_factor, moon_rot, fticks, XY_SCENE_SIZE, OCEAN_DEPTH;


// ocean variables
double const alpha         = 9.81;
double const g             = 9.8;
int    const WaterX        = 60;
int    const WaterY        = 60;
float  const WindSpeed     = 1.1;
double const WindDirection = 200;
double const omega         = 10.0; // area of the water surface part
int    const maxwave       = 20;

// camera check
float const min_z_check(2.0);
float const max_z_check(10.0);
float const min_x1(-1.05*XY_SCENE_SIZE), min_x2(1.0*XY_SCENE_SIZE);
float const min_y1(-0.95*XY_SCENE_SIZE), min_y2(1.1*XY_SCENE_SIZE);
float const max_x1(-2.5*XY_SCENE_SIZE),  max_x2(2.5*XY_SCENE_SIZE);
float const max_y1(-2.5*XY_SCENE_SIZE),  max_y2(2.5*XY_SCENE_SIZE);
float const farx1(-50), farx2(50), fary1(-50), fary2(50);

int ocean_set(0), ocean_draw(0), ocean_do_repeat(1); // must be at least 1
float startx(-50), endx(50), starty(-50), endy(50);
float incrementx, incrementy, incxy;
float ocean_p1, ocean_p2, ocean_p3, ocean_p4;
int   can_see_ocean, can_see_far_ocean;


double const TWO_SQRT_PI = 2.0*sqrt(PI);
double const EPM_CONST   = alpha*g*g/pow(TWO_PI,4);
double const AMP_CONST   = g*PI*PI/omega;
double const lambda      = 5.0*PI/(pow((double)TWO_PI, 2.0)/g);
double const fPeak       = 0.13*g/WindSpeed; // calculate peak frequency


WAVE *w; // never freed
vector3d *WaterNormal;
vector3d *WaterVertNormal;
double *WaterHeight;


void draw_sand(colorRGBA &color, float cscale, int mode);
void draw_ocean2(point &camera, colorRGBA &color, float cscale);


#define fx(x) ((x<0)?(x+WaterX):((x>=WaterX)?(x-WaterX):x))
#define fy(y) ((y<0)?(y+WaterY):((y>=WaterY)?(y-WaterY):y))


void InitWater() { // *** not correct for fticks != 1.0 and X/Y/Z scene size != 1.0

	w               = new WAVE[maxwave];
	WaterNormal     = new vector3d[WaterX*WaterY];
	WaterVertNormal = new vector3d[WaterX*WaterY];
	WaterHeight     = new double[WaterX*WaterY];

	for (int nWaves = 0; nWaves < maxwave; ) {
		float const f(fPeak + (rand()&16383)/16384.0 * 0.5 - 0.5);

		if (f > 0.0) {
			w[nWaves].lambda    = lambda;
			w[nWaves].k         = TWO_PI/w[nWaves].lambda;
			w[nWaves].omega     = sqrt(g*w[nWaves].k);
			w[nWaves].freq      = w[nWaves].omega / TWO_PI;
			w[nWaves].periode   = 1.0/w[nWaves].freq;
			w[nWaves].direction = ((rand()&16383)/16384.0 + 0.5)*PI;
			w[nWaves].dirX		= cosf(w[nWaves].direction)*1.5;
			w[nWaves].dirY		= sinf(w[nWaves].direction/1.7)*0.1;
			float const phi0(((rand()&16383)/16384.0)*PI);
			float const a0(Amplitude(w[nWaves].freq, w[nWaves].direction-WindDirection, w[nWaves].k));
			w[nWaves].amplitude = a0*cosf(phi0);
			w[nWaves].phase     = a0*sinf(phi0);
			if (fabs(w[nWaves].amplitude) >= 0.0001) nWaves ++;
		}
	}
	for (int i = 0; i < WaterX*WaterY; i++) {
		WaterNormal[i] = plus_z;
	}
}

double gamma(double x) {

	--x;
	double faculty(1.0), floor_x(floor(x)), ulim((floor_x == x) ? (x + 1.0) : x);
	for(int i = 1; i < ulim; i++) faculty *= i;
	double const fact(floor_x*x + 1.0 - floor_x*floor_x);
	return faculty*fact;
}

double MitsuyasuDistribution(double f, double theta) {

	double const s((f < fPeak) ? 9.77*pow(f/fPeak, -2.5) : 6.97*pow(f/fPeak, 5.0)); // frequency peak fpeak is globally defined
	double const ang(cos(0.5*theta));
	return gamma(s+1)*pow(ang*ang,s)/(TWO_SQRT_PI*gamma(s+0.5));
}

inline double EPM(double f) {

	return EPM_CONST*exp(-1.25*pow(f/fPeak, -4)); // EPM=Energy according to the Pierson-Moskowitz formula
}

inline double EnergyDistribution(double f, double theta) {

   return EPM(f)*MitsuyasuDistribution(f, theta);
}

inline double Amplitude(double f, double theta, double k) {
              
   return sqrt(EnergyDistribution(f, theta)*AMP_CONST/(k*f)); // amplitude of the wave
}

void AnimateWater() {

	float const time(GET_TIME_MS()/2000.0);

	for (int q = 0; q < WaterX*WaterY; q++) {
		WaterHeight[q] = 0.0;
	}
	for (int k = 0; k < maxwave; k++) {
		float const const_sin(10.0*w[k].amplitude*sin(w[k].phase + w[k].omega*time));
		float const const_cos(10.0*w[k].amplitude*cos(w[k].phase + w[k].omega*time));

		for (int j = 0, q = 0; j < WaterY; j++) {
			float const cos_yj(cos(w[k].omega*w[k].dirY*j)), sin_yj(sin(w[k].omega*w[k].dirY*j));
			float const cos_x(cos(w[k].omega*w[k].dirX)), sin_x(sin(w[k].omega*w[k].dirX));
			float cos_next(1.0), sin_next(0.0);

			for (int i = 0; i < WaterX; i++, q++) {
				float const cos_nx(cos_next), sin_nx(sin_next);
				float const cos_phi(cos_yj*cos_next - sin_yj*sin_next);
				float const sin_phi(sin_yj*cos_next + cos_yj*sin_next);
				cos_next = cos_nx*cos_x - sin_nx*sin_x;
				sin_next = cos_nx*sin_x + sin_nx*cos_x;
				WaterHeight[q] += const_sin*cos_phi + const_cos*sin_phi;
			}
		}
	}
	calc_ocean_normals();
}


void calc_ocean_normals() {

	for (int j = 0; j < WaterY; j++) {
		int const jwx(j*WaterX), fyj(fy(j+1)*WaterX);

		for (int i = 0; i < WaterX; i++) {
			float const wmij(WaterHeight[i+jwx]);
			vector3d nv((WaterHeight[fx(i+1)+jwx] - wmij), (WaterHeight[i+fyj] - wmij), incxy);
			nv.normalize();
			WaterNormal[i+jwx] = nv;

			if (i > 0) {
				nv += WaterNormal[i-1+jwx];
				if (j > 0) nv += WaterNormal[i-1+(jwx - WaterX)];
			}
			if (j > 0)  nv    += WaterNormal[i+(jwx - WaterX)];
			if (i == 0) nv.z  += 2.0;
			if (j == 0) nv.z  += 2.0;
			nv.y = -nv.y;
			WaterVertNormal[i+jwx] = nv*0.25;
		}
	}
}


void set_ocean_z() {

	ocean.z = zmin + OCEAN_DEPTH;
}


void set_view_points() {

	point camera(get_camera_pos());
	can_see_ocean     = 1;
	can_see_far_ocean = 1;
	set_ocean_z();

	if (!camera_mode) { //in the air
		if (camera.z < min_z_check) {
			ocean_p1 = min_x1*2;
			ocean_p2 = min_x2*2;
			ocean_p3 = min_y1*2;
			ocean_p4 = min_y2*2;
		}
		else if (camera.z < max_z_check) {
			ocean_p1 = max_x1;
			ocean_p2 = max_x2;
			ocean_p3 = max_y1;
			ocean_p4 = max_y2;
		}
		else {
			ocean_p1 = farx1;
			ocean_p2 = farx2;
			ocean_p3 = fary1;
			ocean_p4 = fary2;
		}
		return;
	}
	//on the ground
	ocean_p1 = min_x1;
	ocean_p2 = min_x2;
	ocean_p3 = min_y1;
	ocean_p4 = min_y2;
	bool const visible[8] = {
		sphere_in_camera_view(point((min_x1 + 0.01), (min_y1 + 0.01), ocean.z), 0.1, 0),
		sphere_in_camera_view(point((min_x2 - 0.01), (min_y1 + 0.01), ocean.z), 0.1, 0),
		sphere_in_camera_view(point((min_x2 - 0.01), (min_y2 - 0.01), ocean.z), 0.1, 0),
		sphere_in_camera_view(point((min_x1 + 0.01), (min_y2 - 0.01), ocean.z), 0.1, 0),
		sphere_in_camera_view(point((-1 + 0.01),     (-1 + 0.01),     ocean.z), 0.1, 0),
		sphere_in_camera_view(point(( 1 - 0.01),     (-1 + 0.01),     ocean.z), 0.1, 0),
		sphere_in_camera_view(point(( 1 - 0.01),     ( 1 - 0.01),     ocean.z), 0.1, 0),
		sphere_in_camera_view(point((-1 + 0.01),     ( 1 - 0.01),     ocean.z), 0.1, 0)
	};
	if (visible[0] || visible[4]) {
		if (visible[5] || visible[7]) return;
		ocean_p1 = ocean_p3 = min_x1;
		ocean_p2 = ocean_p4 = min_y1;
	}
	else if (visible[1] || visible[5]) {
		if (visible[4] || visible[6]) return;
		ocean_p1 = ocean_p3 = min_x2;
		ocean_p2 = ocean_p4 = min_y1;
	}
	else if (visible[2] || visible[6]) {
		if (visible[5] || visible[7]) return;
		ocean_p1 = ocean_p3 = min_x2;
		ocean_p2 = ocean_p4 = min_y2;
	}
	else if (visible[3] || visible[7]) {
		if (visible[4] || visible[6]) return;
		ocean_p1 = ocean_p3 = min_x1;
		ocean_p2 = ocean_p4 = min_y2;
	}
	else {
		ocean_p1 = ocean_p2 = ocean_p3 = ocean_p4 = 0.0;
	}
}


void draw_sand(colorRGBA &color, float cscale, int mode) {

	float const oz3(zmin + OCEAN_DEPTH - SAND_DEPTH);
	(DISABLE_WATER ? WHITE : colorRGBA(0.25*cscale, 0.5*cscale, cscale, 1.0)).do_glColor();
	plus_z.do_glNormal();
	select_texture(SAND_TEX);
	setup_texgen(SAND_REPEAT, SAND_REPEAT, 0.0, 0.0);
	glBegin(GL_QUADS);
	if (mode == 1) draw_one_tquad(-0.25*ocean.x, -0.25*ocean.y, 0.25*ocean.x, 0.25*ocean.y, oz3, 0);
	
	for (unsigned i = 0; i < 4; ++i) {
		float const s0((-1.0 + 2.0*(i == 1 || i == 2))*(X_SCENE_SIZE - DX_VAL));
		float const s1((-1.0 + 2.0*(i < 2))*(Y_SCENE_SIZE - DY_VAL));
		glVertex3f( 4.0*s0, 4.0*s1, oz3);
		glVertex3f(-4.0*s1, 4.0*s0, oz3);
		glVertex3f(    -s1,     s0, zmin);
		glVertex3f(     s0,     s1, zmin);
	}
	{
		float val1(-X_SCENE_SIZE);
		float limit(-Y_SCENE_SIZE+DY_VAL);

		for (unsigned i = 1; i < (unsigned)MESH_X_SIZE; ++i) { // -y edge
			float const val2(val1 + DX_VAL);
			glVertex3f(val1, limit, zmin);
			glVertex3f(val2, limit, zmin);
			glVertex3f(val2, limit, mesh_height[1][i]);
			glVertex3f(val1, limit, mesh_height[1][i-1]);
			val1 = val2;
		}
	}
	{
		float val1(-Y_SCENE_SIZE);
		float limit(-X_SCENE_SIZE+DX_VAL);

		for (unsigned i = 1; i < (unsigned)MESH_Y_SIZE; ++i) { // -x edge
			float const val2(val1 + DY_VAL);
			glVertex3f(limit, val1, zmin);
			glVertex3f(limit, val2, zmin);
			glVertex3f(limit, val2, mesh_height[i][1]);
			glVertex3f(limit, val1, mesh_height[i-1][1]);
			val1 = val2;
		}
	}
	glEnd();
	disable_textures_texgen();

	if (!DISABLE_WATER) { // draw ocean barrier as well (upwards normal)
		glColor3f(0.3*color.R, color.G, color.B);
		draw_fast_cylinder(point(0.0, 0.0, oz3), point(0.0, 0.0, oz3+SAND_DEPTH-MWAVE_HEIGHT), 0.25*FAR_CLIP, 0.25*FAR_CLIP, N_CYL_SIDES, 0);
	}
}


float scale_color(colorRGBA &color) {

	float cscale;
	apply_red_sky(color);

	if (light_factor >= 0.6) { // sun
		cscale = light_factor;
	}
	else if (light_factor <= 0.4) { // moon
		cscale = 0.7*get_moon_light_factor();
	}
	else { // sun + moon
		cscale = 5.0*(light_factor*(light_factor - 0.4) + 0.7*get_moon_light_factor()*(0.6 - light_factor));
	}
	color *= cscale;
	return cscale;
}


void update_incs() {

	incrementx = ((float) (endx - startx))/WaterX;
	incrementy = ((float) (endy - starty))/WaterY;
	incxy      = incrementx*incrementx + incrementy*incrementy;
}


void set_ocean_alpha(colorRGBA &color, float &last_alpha, float zscale, int i, int j) {

	if (FORCE_ALPHA_1) {
		color.alpha = 1.0;
		return;
	}
	if (point_outside_mesh(i, j)) {
		color.alpha = 1.0;
	}
	else {
		color.alpha = max(0.4f, min(1.0f, (0.5f + zscale*(ocean.z - mesh_height[i][j]))));
	}
	if (fabs(last_alpha - color.alpha) > 0.02) {
		color.do_glColor();
		last_alpha = color.alpha;
	}
}


void draw_vertex(int index, float x, float y, float z, colorRGBA &color, float &last_alpha, float zscale) {

	WaterVertNormal[index].do_glNormal();
	set_ocean_alpha(color, last_alpha, zscale, int((y+Y_SCENE_SIZE)*DY_VAL_INV), int((x+X_SCENE_SIZE)*DX_VAL_INV));
	glVertex3f(x, y, (z + ocean.z + WaterHeight[index]));
}


void draw_ocean() {

	if (!island) {
		ocean_set = 0;
		return;
	}
	update_incs();
	float last_alpha(1.0);
	colorRGBA color(0.1, 0.25, 0.7, 1.0); // ocean color
	if (display_mode & 0x0100) color.G += 0.1;
	point camera(get_camera_pos()), pt(all_zeros);
	float const cscale(scale_color(color));

	if (DISABLE_WATER) {
		draw_sand(color, cscale, (camera.z <= ocean.z + 0.1));
		return;
	}
	if (camera.z < ocean.z) gen_bubble(camera);
	if (ocean_set == 0) InitWater();
	ocean_set = 1;
	pt.z     -= 0.3*(FAR_CLIP+X_SCENE_SIZE);
	ocean.assign((FAR_CLIP + X_SCENE_SIZE), (FAR_CLIP + Y_SCENE_SIZE), (zmin + OCEAN_DEPTH));
	if (distance_to_camera(pt) > (FAR_CLIP+X_SCENE_SIZE)/2.0) return;
	set_lighted_sides(2);
	set_fill_mode();

	if (display_mode & 0x0100) {
		color.do_glColor();
		draw_ocean2(camera, color, cscale);
		set_lighted_sides(1);
		plus_z.do_glNormal();
		return;
	}
	if (ocean_draw == 0) AnimateWater();
	set_view_points();
	if (!can_see_ocean && !can_see_far_ocean) return;
	float const zscale(0.5/min(MAX_VIS_DEPTH, ocean.z - zmin));
	glDisable(GL_LIGHTING);
	draw_sand(color, cscale, (camera.z <= ocean.z + 0.1));
	glEnable(GL_LIGHTING);
	enable_blend();
	if (show_fog) glDisable(GL_FOG);
	select_texture(WATER_TEX);
	plus_z.do_glNormal();
	glEnable(GL_COLOR_MATERIAL);
	color.do_glColor();
	setup_texgen(OCEAN_REPEAT, OCEAN_REPEAT, 0.0, 0.0);
	int min_startx(0), min_endx(WaterX), min_starty(0), min_endy(WaterY);
	
	for (int ocean_id = FIRST_OCEAN; ocean_id < LAST_OCEAN; ocean_id++) {
		if (ocean_id == 0 && can_see_ocean) {
			if (!camera_mode) {
				startx = camera.x + ocean_p1;
				endx   = camera.x + ocean_p2;
				starty = camera.y + ocean_p3;
				endy   = camera.y + ocean_p4;
			}
			else {
				startx = min_x1;
				endx   = min_x2;
				starty = min_y1;
				endy   = min_y2;
			}
		}
		else if (ocean_id == 1 && can_see_far_ocean) {
			startx = farx1;
			endx   = farx2;
			starty = fary1;
			endy   = fary2;
		}
		else {
			continue;
		}
		update_incs();

		//draw the detailed ocean
		glBegin(GL_QUADS);

		for (float j = min_starty; j < min_endy; j++) {
			float y(incrementy*j + starty);
			int const jwx(j*WaterX), fyj(fy(j+1)*WaterX);
			
			for (int i = min_startx; i < min_endx; i++) {
				float const x(incrementx*i + startx);
				draw_vertex(i+jwx,       x,            y,            0.0, color, last_alpha, zscale);
				draw_vertex(fx(i+1)+jwx, x+incrementx, y,            0.0, color, last_alpha, zscale);
				draw_vertex(fx(i+1)+fyj, x+incrementx, y+incrementy, 0.0, color, last_alpha, zscale);
				draw_vertex(i+fyj,       x,            y+incrementy, 0.0, color, last_alpha, zscale);
			}
		}
		glEnd();
	}
	disable_textures_texgen();
	glDisable(GL_COLOR_MATERIAL);
	disable_blend();
	set_lighted_sides(1);
	if (show_fog) glEnable(GL_FOG);
	++ocean_draw;
	ocean_draw = (ocean_draw % ocean_do_repeat);
}


void draw_ocean2(point &camera, colorRGBA &color, float cscale) {

	float last_alpha(1.0);
	static int lfc(0);
	static float phase(0.0);

	if (MOVING_OCEAN) {
		phase   += fticks*MWAVE_PERIOD*float(frame_counter - lfc);
		ocean.z += MWAVE_HEIGHT*sin(phase);
		lfc      = frame_counter;
	}
	glDisable(GL_LIGHTING);

	if (SIMPLE_OCEAN) {
		plus_z.do_glNormal();
		if (show_fog) glDisable(GL_FOG);
		select_texture(WATER_TEX);
		color.do_glColor();
		draw_tquad(ocean.x, ocean.y, ocean.z, 1, 0.0, 0.0, 100.0*OCEAN_REPEAT, 100.0*OCEAN_REPEAT);
		if (show_fog) glEnable(GL_FOG);
		if (camera.z <= ocean.z + 0.1) draw_sand(color, cscale, 1);
		glDisable(GL_TEXTURE_2D);
		glEnable(GL_LIGHTING);
		return;
	}
	draw_sand(color, cscale, (camera.z <= ocean.z + 0.1));
	if (show_fog) glDisable(GL_FOG);
	select_texture(WATER_TEX);
	setup_texgen(OCEAN_REPEAT, OCEAN_REPEAT, 0.0, 0.0);
	float const OCEAN_SKEW_X(OCEAN_SKEW*X_SCENE_SIZE), OCEAN_SKEW_Y(OCEAN_SKEW*Y_SCENE_SIZE);
	enable_blend();
	color.alpha = last_alpha;
	color.do_glColor();
	glPushMatrix();
	glTranslatef(0.0, 0.0, ocean.z);
	glBegin(GL_QUADS);
	glVertex2f(-ocean.x,      -ocean.y);
	glVertex2f( ocean.x,      -ocean.y);
	glVertex2f( ocean.x,      -OCEAN_SKEW_Y);
	glVertex2f(-ocean.x,      -Y_SCENE_SIZE);
	glVertex2f(-ocean.x,       OCEAN_SKEW_Y);
	glVertex2f( ocean.x,       Y_SCENE_SIZE);
	glVertex2f( ocean.x,       ocean.y);
	glVertex2f(-ocean.x,       ocean.y);
	glVertex2f(-ocean.x,      -Y_SCENE_SIZE);
	glVertex2f(-OCEAN_SKEW_X, -Y_SCENE_SIZE);
	glVertex2f(-X_SCENE_SIZE,  Y_SCENE_SIZE);
	glVertex2f(-ocean.x,       Y_SCENE_SIZE);
	glVertex2f( ocean.x,       Y_SCENE_SIZE);
	glVertex2f( OCEAN_SKEW_X,  Y_SCENE_SIZE);
	glVertex2f( X_SCENE_SIZE, -Y_SCENE_SIZE);
	glVertex2f( ocean.x,      -Y_SCENE_SIZE);
	glEnd();

	glBegin(GL_TRIANGLE_STRIP);
	float yval(-Y_SCENE_SIZE);
	float const zscale(0.5/min(MAX_VIS_DEPTH, ocean.z - zmin));

	for (int i = 1; i <= MESH_Y_SIZE; ++i) {
		float xval(-X_SCENE_SIZE);
		bool new_strip(1);

		for (int j = 1; j <= MESH_X_SIZE; ++j) {
			if (i == MESH_Y_SIZE || j == MESH_X_SIZE || z_min_matrix[i][j] <= ocean.z || (i > 0 && j > 0 && z_min_matrix[i-1][j-1] <= ocean.z)) {
				if (new_strip) {
					glEnd();
					glBegin(GL_TRIANGLE_STRIP);
					color.do_glColor();
					new_strip = 0;
					set_ocean_alpha(color, last_alpha, zscale, i, j);
					glVertex2f(xval, yval);
					set_ocean_alpha(color, last_alpha, zscale, i+1, j);
					glVertex2f(xval, (yval+DY_VAL));
				}
				set_ocean_alpha(color, last_alpha, zscale, i, j+1);
				glVertex2f((xval+DX_VAL), yval);
				set_ocean_alpha(color, last_alpha, zscale, i+1, j+1);
				glVertex2f((xval+DX_VAL), (yval+DY_VAL));
			}
			else {
				new_strip = 1;
			}
			xval += DX_VAL;
		}
		yval += DY_VAL;
	}
	glEnd();
	disable_blend();
	disable_textures_texgen();
	glEnable(GL_LIGHTING);
	glPopMatrix();
	if (show_fog) glEnable(GL_FOG);
}



