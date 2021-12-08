// 3D World - Video Capture using ffmpeg
// by Frank Gennari
// 10/26/15
#include "3DWorld.h"
#include <list>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cstring> // for memcpy()
#include <chrono>

using namespace std;

unsigned const MAX_FRAMES_BUFFERED = 128;

extern int window_width, window_height;
extern unsigned video_framerate; // Note: should probably be either 30 or 60
extern unsigned num_video_threads; // defaults to 0 = max

void sleep_for_ms(unsigned milliseconds) {
	std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
}

void write_video();

// from http://vichargrave.com/multithreaded-work-queue-in-c/
template<typename T> class thread_safe_queue {
	list<T> m_queue;
	mutable std::mutex m_mutex;
	std::condition_variable m_condv;
public:
	void add(T const &item) {
		std::unique_lock<std::mutex> mlock(m_mutex);
		m_queue.push_back(item);
		mlock.unlock(); // unlock before notificiation to minimize mutex contention
		m_condv.notify_one(); // notify one waiting thread
	}
	T remove() {
		std::unique_lock<std::mutex> mlock(m_mutex);
		while (m_queue.empty()) {m_condv.wait(mlock);}
		T const item(m_queue.front());
		m_queue.pop_front();
		return item;
	}
	size_t size() const {
		std::unique_lock<std::mutex> mlock(m_mutex);
		size_t const size(m_queue.size());
		return size;
	}
	bool empty() const {
		std::unique_lock<std::mutex> mlock(m_mutex);
		bool const ret(m_queue.empty());
		return ret;
	}
};

class video_capture_t {

	class video_buffer {
		struct frame_t : public vector<unsigned char> {
			unsigned data_sz;
			frame_t(unsigned data_sz_) : vector<unsigned char>(data_sz_), data_sz(data_sz_) {}
		};
		typedef shared_ptr<frame_t> p_frame_t;
		thread_safe_queue<p_frame_t> frames, free_list; // queue of frames to compress

	public:
		void push_frame(void const *const data, unsigned data_sz) {
			assert(data_sz > 0);
			p_frame_t frame(free_list.empty() ? p_frame_t(new frame_t(data_sz)) : free_list.remove()); // take a frame from the free list if nonempty
			assert(frame->data_sz == data_sz); // free_list returned must be of the correct size
			memcpy(&frame->front(), data, data_sz);
			frames.add(frame);
		}
		void pop_and_send_frame(FILE *fp) {
			assert(fp != nullptr);
			p_frame_t const frame(frames.remove());
			fwrite(&frame->front(), frame->size(), 1, fp);
			if (free_list.size() < MAX_FRAMES_BUFFERED/4) {free_list.add(frame);} // recycle it
		}
		void write_frames(FILE *fp) {
			while (!frames.empty()) {pop_and_send_frame(fp);}
		}
		void empty_free_list() {
			// must be called on window resize, which is illegal to do when recording; okay/recommended to call between videos
			while (!free_list.empty()) {free_list.remove();}
		}
		size_t num_pending_frames() const {return frames.size();}
		bool empty() const {return frames.empty();}
	};

	unsigned video_id, pbo, data_sz;
	string filename;

	// multithreaded writing support
	bool is_recording, is_writing;
	video_buffer buffer;
	unique_ptr<std::thread> write_thread;

	void wait_for_write_complete() {
		if (!write_thread) return;
		is_recording = 0;
		if (is_writing) {cout << "Wating for " << buffer.num_pending_frames() << " video frames to be written" << endl;}
		write_thread->join();
		write_thread.reset();
		assert(!is_writing);
	}
	void queue_frame(void const *const data) {
		assert(data_sz == get_num_bytes());
		buffer.push_frame(data, data_sz);
		
		while (buffer.num_pending_frames() > MAX_FRAMES_BUFFERED) { // sleep for 10ms until buffer is partially emptied
			cout << "Waiting for video write buffer to empty" << endl;
			sleep_for_ms(10);
		}
	}
	static unsigned get_num_bytes() {return 4*window_width*window_height;}

public:
	video_capture_t() : video_id(0), pbo(0), data_sz(0), is_recording(0), is_writing(0) {}

	void start(string const &fn) {
		assert(!is_recording); // must end() before calling start() again
		wait_for_write_complete();
		assert(!is_writing);
		is_recording = 1;
		data_sz      = get_num_bytes();
		assert(pbo == 0);
		glGenBuffers(1, &pbo);
		glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo);
		glBufferData(GL_PIXEL_PACK_BUFFER, data_sz, NULL, GL_STREAM_READ);
		glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
		// start writing in a different thread
		filename = fn;
		assert(!write_thread);
		write_thread.reset(new std::thread(write_video));
	}
	void write_buffer() {
		assert(!filename.empty());
		// start ffmpeg telling it to expect raw RGBA, 60 FPS
		// -i - tells it to read frames from stdin
		// Note: 0 = max threads; the more threads the lower the frame rate, as video compression competes with 3DWorld for CPU cycles;
		// however, more threads is less likely to fill the buffer and block, producing heavy lag
		ostringstream oss;
		oss << " -r " << video_framerate << " -f rawvideo -pix_fmt rgba -s " << window_width << "x" << window_height
			<< " -i - -threads " << num_video_threads << " -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip " << filename;
		// open pipe to ffmpeg's stdin in binary write mode
#ifdef _WIN32
		string const cmd(string("ffmpeg.exe.lnk") + oss.str());
		FILE* ffmpeg = _popen(cmd.c_str(), "wb");
#else
		string const cmd(string("ffmpeg") + oss.str());
		FILE* ffmpeg = popen(cmd.c_str(), "w");
#endif
		if(ffmpeg == nullptr) {
		  cerr << "Error running ffmpeg command: " << cmd << endl;
		  is_writing = 0;
		  end();
		  return;
		}
		is_writing = 1;
		while (is_recording || !buffer.empty()) {buffer.write_frames(ffmpeg); sleep_for_ms(1);} // 1ms sleep
		buffer.empty_free_list(); // clear free list (to save memory, and in case the frame buffer size changes between recording sessions)
#ifdef _WIN32
		_pclose(ffmpeg);
#else
		pclose(ffmpeg);
#endif
		is_writing = 0;
	}
	void end() {
		is_recording = 0; // signal writer to finish
		glDeleteBuffers(1, &pbo);
		pbo = 0;
	}
	void toggle_start_stop() {
		if (is_recording) {end(); return;} // start=>end
		ostringstream oss;
		oss << "video_out" << video_id++ << ".mp4";
		start(oss.str()); // end=>start
	}
	void end_frame() {
		if (!is_recording) return;
		assert(pbo != 0);
		assert(data_sz == get_num_bytes()); // make sure the resolution hasn't changed since recording started
		//timer_t timer("Video Capture Frame"); // 13.7ms for 1920x1024, 10.9ms with free list
		glReadBuffer(GL_FRONT);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo);
		glReadPixels(0, 0, window_width, window_height, GL_RGBA, GL_UNSIGNED_BYTE, nullptr); // use PBO
		void *ptr = glMapBufferRange(GL_PIXEL_PACK_BUFFER, 0, data_sz, GL_MAP_READ_BIT); // this line takes most of the time
		queue_frame(ptr);
		glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
		glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
	}
	bool is_video_recording() const {return is_recording;}
	~video_capture_t() {wait_for_write_complete();} // wait for write to complete; don't try to free the pbo
};

video_capture_t video_capture;

// Note: must be a global function rather than member function for thread constructor
void write_video() {video_capture.write_buffer();}

// Note: not legal to resize the window between start() and end()
void start_video_capture(string const &fn) {video_capture.start(fn);}
void end_video_capture() {video_capture.end();}
void toggle_video_capture() {video_capture.toggle_start_stop();}
void video_capture_end_frame() {video_capture.end_frame();}
bool is_video_recording() {return video_capture.is_video_recording();}

