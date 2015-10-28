// 3D World - Video Capture using ffmpeg
// by Frank Gennari
// 10/26/15
#include "3DWorld.h"
#include "openal_wrap.h" // for alut_sleep()
#include <pthread.h>
#include <list>

using namespace std;

bool const USE_WRITE_THREAD = 1;

extern int window_width, window_height;

void *write_video(void *data);

// from http://vichargrave.com/multithreaded-work-queue-in-c/
template<typename T> class thread_safe_queue {
	list<T> m_queue;
	mutable pthread_mutex_t m_mutex;
	pthread_cond_t m_condv;

public:
	thread_safe_queue() {
		pthread_mutex_init(&m_mutex, nullptr);
		pthread_cond_init (&m_condv, nullptr);
	}
	~thread_safe_queue() {
		pthread_mutex_destroy(&m_mutex);
		pthread_cond_destroy (&m_condv);
	}
	void add(T const &item) {
		pthread_mutex_lock(&m_mutex);
		m_queue.push_back(item);
		pthread_cond_signal(&m_condv);
		pthread_mutex_unlock(&m_mutex);
	}
	T remove() {
		pthread_mutex_lock(&m_mutex);
		while (m_queue.size() == 0) {pthread_cond_wait(&m_condv, &m_mutex);}
		T const item(m_queue.front());
		m_queue.pop_front();
		pthread_mutex_unlock(&m_mutex);
		return item;
	}
	size_t size() const {
        pthread_mutex_lock(&m_mutex);
        size_t const size(m_queue.size());
        pthread_mutex_unlock(&m_mutex);
        return size;
    }
	bool empty() const {
        pthread_mutex_lock(&m_mutex);
        bool const ret(m_queue.empty());
        pthread_mutex_unlock(&m_mutex);
        return ret;
    }
};

class video_capture_t {

	class video_buffer {
		typedef vector<unsigned char> frame_t;
		typedef shared_ptr<frame_t> p_frame_t;
		thread_safe_queue<p_frame_t> frames; // queue of frames to compress

	public:
		void push_frame(void const *const data, unsigned data_sz) {
			assert(data_sz > 0);
			p_frame_t frame(new frame_t(data_sz));
			memcpy(&frame->front(), data, data_sz);
			frames.add(frame);
		}
		void pop_and_send_frame(FILE *fp) {
			assert(fp != nullptr);
			p_frame_t const frame(frames.remove());
			fwrite(&frame->front(), frame->size(), 1, fp);
		}
		void write_frames(FILE *fp) {
			while (!frames.empty()) {pop_and_send_frame(fp);}
		}
		size_t num_pending_frames() const {return frames.size();}
		bool empty() const {return frames.empty();}
	};

	unsigned video_id, pbo, start_sz;
	FILE* ffmpeg;

	// multithreaded writing support
	bool end_data, write_active;
	video_buffer buffer;
	pthread_t write_thread;

	void wait_for_write_complete() {
		end_data = 1;
		cout << "Wating for " << buffer.num_pending_frames() << " frames to be written" << endl;
		pthread_join(write_thread, nullptr);
		assert(!write_active);
		end_data = 0;
	}
	void queue_frame(void const *const data) {buffer.push_frame(data, get_num_bytes());}

	static unsigned get_num_bytes() {return 4*window_width*window_height;}

public:
	video_capture_t() : video_id(0), pbo(0), start_sz(0), ffmpeg(nullptr), end_data(0), write_active(0) {}

	void start(string const &fn) {
		assert(ffmpeg == nullptr); // must end() before calling start() again
		assert(pbo == 0);
		assert(!write_active);
		// start ffmpeg telling it to expect raw RGBA, 60 FPS
		// -i - tells it to read frames from stdin
		ostringstream oss;
		oss << "ffmpeg.exe.lnk -r 60 -f rawvideo -pix_fmt rgba -s " << window_width << "x" << window_height << " -i - -threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip " << fn;
		// open pipe to ffmpeg's stdin in binary write mode
		start_sz = get_num_bytes();
		ffmpeg   = _popen(oss.str().c_str(), "wb");
		glGenBuffers(1, &pbo);
		glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo);
		glBufferData(GL_PIXEL_PACK_BUFFER, start_sz, NULL, GL_STREAM_READ);
		glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

		if (USE_WRITE_THREAD) { // start in a different thread
			pthread_attr_t attr;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			int const rc(pthread_create(&write_thread, &attr, write_video, nullptr)); 
			if (rc) {cout << "Error: Return code from pthread_create() is " << rc << endl; assert(0);}
			pthread_attr_destroy(&attr);
		}
	}
	void write_buffer() {
		write_active = 1;
		while (!end_data) {buffer.write_frames(ffmpeg);}
		write_active = 0;
	}
	void end(bool called_from_dtor=0) {
		if (ffmpeg == nullptr) return;
		// Note: if a write thread is active, this will block until writing is complete
		// FIXME: if !called_from_dtor, we may want to let writing continue in the background, but force it to finish if we plan to start recording again
		if (USE_WRITE_THREAD) {wait_for_write_complete();} // only needs to wait if the write thread is active
		_pclose(ffmpeg);
		ffmpeg = nullptr;
		if (called_from_dtor) return; // don't try to free the pbo
		glDeleteBuffers(1, &pbo);
		pbo = 0;
	}
	void toggle_start_stop() {
		if (ffmpeg != nullptr) {end(); return;} // start=>end
		ostringstream oss;
		oss << "video_out" << video_id++ << ".mp4";
		start(oss.str()); // end=>start
	}
	void end_frame() {
		if (ffmpeg == nullptr) return;
		assert(pbo != 0);
		assert(start_sz == get_num_bytes());
		//RESET_TIME;
		glReadBuffer(GL_FRONT);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo);
		glReadPixels(0, 0, window_width, window_height, GL_RGBA, GL_UNSIGNED_BYTE, nullptr); // use PBO
		void *ptr = glMapBufferRange(GL_PIXEL_PACK_BUFFER, 0, get_num_bytes(), GL_MAP_READ_BIT);
		//PRINT_TIME("Read");
		if (write_active) {queue_frame(ptr);} else {fwrite(ptr, get_num_bytes(), 1, ffmpeg);} // queue it or write it directly
		//PRINT_TIME("Write");
		glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
		glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
	}
	~video_capture_t() {end(1);}
};

video_capture_t video_capture;

// Note: must be a global function rather than member function for pthread_create(); argument and return value are unused
void *write_video(void *data) {video_capture.write_buffer(); return nullptr;}

// Note: not legal to resize the window between start() and end()
void start_video_capture(string const &fn) {video_capture.start(fn);}
void end_video_capture() {video_capture.end();}
void toggle_video_capture() {video_capture.toggle_start_stop();}
void video_capture_end_frame() {video_capture.end_frame();}

