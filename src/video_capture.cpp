// 3D World - Video Capture using ffmpeg
// by Frank Gennari
// 10/26/15
#include "3DWorld.h"

using std::string;

extern int window_width, window_height;


class video_capture_t {

	FILE* ffmpeg;
	unsigned video_id, pbo, start_sz;

	static unsigned get_num_bytes() {return 4*window_width*window_height;}

public:
	video_capture_t() : ffmpeg(nullptr), video_id(0), pbo(0), start_sz(0) {}

	void start(std::string const &fn) {
		assert(ffmpeg == nullptr); // must end() before calling start() again
		assert(pbo == 0);
		// start ffmpeg telling it to expect raw RGBA, 60 FPS
		// -i - tells it to read frames from stdin
		std::ostringstream oss;
		oss << "ffmpeg.exe.lnk -r 60 -f rawvideo -pix_fmt rgba -s " << window_width << "x" << window_height << " -i - -threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip " << fn;
		// open pipe to ffmpeg's stdin in binary write mode
		start_sz = get_num_bytes();
		ffmpeg   = _popen(oss.str().c_str(), "wb");
		glGenBuffers(1, &pbo);
		glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo);
		glBufferData(GL_PIXEL_PACK_BUFFER, start_sz, NULL, GL_STREAM_READ);
	}
	void end(bool called_from_dtor=0) {
		if (ffmpeg == nullptr) return;
		glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
		_pclose(ffmpeg);
		ffmpeg = nullptr;
		if (called_from_dtor) return; // don't try to free the pbo
		glDeleteBuffers(1, &pbo);
		pbo = 0;
	}
	void toggle_start_stop() {
		if (ffmpeg != nullptr) {end(); return;} // start=>end
		std::ostringstream oss;
		oss << "video_out" << video_id++ << ".mp4";
		start(oss.str()); // end=>start
	}
	void end_frame() {
		if (ffmpeg == nullptr) return;
		//RESET_TIME;
		assert(start_sz == get_num_bytes());
		glReadBuffer(GL_FRONT);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glReadPixels(0, 0, window_width, window_height, GL_RGBA, GL_UNSIGNED_BYTE, nullptr); // use PBO
		void *ptr = glMapBufferRange(GL_PIXEL_PACK_BUFFER, 0, get_num_bytes(), GL_MAP_READ_BIT);
		//PRINT_TIME("Read");
		fwrite(ptr, get_num_bytes(), 1, ffmpeg);
		glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
		//PRINT_TIME("Write");
	}
	~video_capture_t() {end(1);}
};

video_capture_t video_capture;

// Note: not legal to resize the window between start() and end()
void start_video_capture(std::string const &fn) {video_capture.start(fn);}
void end_video_capture() {video_capture.end();}
void toggle_video_capture() {video_capture.toggle_start_stop();}
void video_capture_end_frame() {video_capture.end_frame();}

