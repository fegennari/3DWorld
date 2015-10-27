// 3D World - Video Capture using ffmpeg
// by Frank Gennari
// 10/26/15
#include "3DWorld.h"

using std::string;

extern int window_width, window_height;


class video_capture_t {

	FILE* ffmpeg;
	vector<char> buffer;
	unsigned video_id;

public:
	video_capture_t() : ffmpeg(nullptr), video_id(0) {}

	void start(std::string const &fn) {
		assert(ffmpeg == nullptr); // must end() before calling start() again
		// start ffmpeg telling it to expect raw RGBA, 30 FPS (60 is too high at fullscreen)
		// -i - tells it to read frames from stdin
		std::ostringstream oss;
		oss << "ffmpeg.exe.lnk -r 30 -f rawvideo -pix_fmt rgba -s " << window_width << "x" << window_height << " -i - -threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip " << fn;
		// open pipe to ffmpeg's stdin in binary write mode
		ffmpeg = _popen(oss.str().c_str(), "wb");
	}
	void end() {
		if (ffmpeg == nullptr) return;
		_pclose(ffmpeg);
		ffmpeg = nullptr;
		buffer.clear();
	}
	void toggle_start_stop() {
		if (ffmpeg != nullptr) {end(); return;} // start=>end
		std::ostringstream oss;
		oss << "video_out" << video_id++ << ".mp4";
		start(oss.str()); // end=>start
	}
	void end_frame() {
		if (ffmpeg == nullptr) return;
		buffer.resize(4*window_width*window_height);
		glReadBuffer(GL_FRONT);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glReadPixels(0, 0, window_width, window_height, GL_RGBA, GL_UNSIGNED_BYTE, &buffer.front()); // ffmpeg wants data in RGBA, not the normal RGB used by read_pixels()
		fwrite(&buffer.front(), sizeof(char), buffer.size(), ffmpeg);
	}
	~video_capture_t() {end();}
};

video_capture_t video_capture;

void start_video_capture(std::string const &fn) {video_capture.start(fn);}
void end_video_capture() {video_capture.end();}
void toggle_video_capture() {video_capture.toggle_start_stop();}
void video_capture_end_frame() {video_capture.end_frame();}

