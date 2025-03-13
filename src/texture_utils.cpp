// 3D World - Texture Compression and Mipmap Utility Functions
// by Frank Gennari
// 4/3/22
#include "3DWorld.h"
#include "function_registry.h"
//#include "profiler.h"

#define STB_DXT_IMPLEMENTATION
#include "stb_dxt.h"


void dxt_texture_compress(uint8_t const *const data, vector<uint8_t> &comp_data, int width, int height, int ncolors) {
	//highres_timer_t timer("stb_dxt Texture Compress", 1, 1); // enabled, no loading screen
	assert(width > 0 && height > 0);
	assert(ncolors == 3 || ncolors == 4);
	assert(data != nullptr);
	// RGB=DXT1/BC1, RGBA=DXT5/BC3
	bool const has_alpha(ncolors == 4);
	unsigned const block_sz(has_alpha ? 16 : 8), x_blocks((width + 3)/4), y_blocks((height + 3)/4); // take ceil()
	comp_data.resize(x_blocks*y_blocks*block_sz);

#pragma omp parallel for schedule(static)
	for (int y = 0; y < height; y += 4) {
		uint8_t block[4*4*4] = {};

		for (int x = 0; x < width; x += 4) {
			for (int yy = 0; yy < 4; ++yy) {
				for (int xx = 0; xx < 4; ++xx) {
					unsigned const bix(4*(4*yy + xx));
					// clamp to valid input texture range in case width and height are not a multiple of 4, which duplicates rows and columns
					unsigned const dix(ncolors*(width*min(y+yy, height-1) + min(x+xx, width-1)));
					for (int c = 0; c < ncolors; ++c) {block[bix + c] = data[dix + c];}
					if (!has_alpha) {block[bix + 3] = 255;} // set alpha=255
				}
			}
			unsigned const comp_offset(((y/4)*x_blocks + (x/4))*block_sz);
			assert(comp_offset < comp_data.size());
			stb_compress_dxt_block(&comp_data[comp_offset], block, has_alpha, /*STB_DXT_NORMAL*/STB_DXT_HIGHQUAL);
		} // for x
	} // for y
}

void texture_t::compress_and_send_texture() {
	//highres_timer_t timer("compress_and_send_texture", 1, 1); // enabled, no loading screen; 680ms
	vector<uint8_t> comp_data; // reuse across calls doesn't seem to help much
	dxt_texture_compress(data, comp_data, width, height, ncolors);
	GL_CHECK(glCompressedTexImage2D(GL_TEXTURE_2D, 0, calc_internal_format(), width, height, 0, comp_data.size(), comp_data.data());)
}

void texture_t::create_compressed_mipmaps() {
	//highres_timer_t timer("create_compressed_mipmaps", 1, 1); // enabled, no loading screen; 1350ms total for city + cars + people
	assert(is_allocated());
	vector<uint8_t> idatav, odata, comp_data; // reuse across calls doesn't seem to help much

	for (unsigned w = width, h = height, level = 1; w > 1 || h > 1; w >>= 1, h >>= 1, ++level) {
		unsigned const w1(max(w, 1U)), h1(max(h, 1U)), w2(max(w>>1, 1U)), h2(max(h>>1, 1U));
		unsigned const xinc((w2 < w1) ? ncolors : 0), yinc((h2 < h1) ? ncolors*w1 : 0);
		uint8_t const *const idata((level == 1) ? data : idatav.data());
		odata.resize(ncolors*w2*h2);

#pragma omp parallel for schedule(static)
		for (int y = 0; y < (int)h2; ++y) { // simple 2x2 box filter
			for (int x = 0; x < (int)w2; ++x) {
				unsigned const ix1(ncolors*(y*w2+x)), ix2(ncolors*((y<<1)*w1+(x<<1)));

				for (int n = 0; n < ncolors; ++n) {
					odata[ix1+n] = uint8_t(((unsigned)idata[ix2+n] + idata[ix2+xinc+n] + idata[ix2+yinc+n] + idata[ix2+yinc+xinc+n]) >> 2);
				}
			}
		}
		dxt_texture_compress(odata.data(), comp_data, w2, h2, ncolors);
		GL_CHECK(glCompressedTexImage2D(GL_TEXTURE_2D, level, calc_internal_format(), w2, h2, 0, comp_data.size(), comp_data.data());)
		idatav.swap(odata);
	} // for level
}

void texture_t::create_custom_mipmaps() { // 350ms total for city + cars + people

	//highres_timer_t timer("Create Custom Mipmaps");
	assert(is_allocated());
	vector<uint8_t> idatav, odata, comp_data;
	color_wrapper cw(color); // for use_mipmaps == 4 with RGBA

	for (unsigned w = width, h = height, level = 1; w > 1 || h > 1; w >>= 1, h >>= 1, ++level) {
		unsigned const w1(max(w, 1U)), h1(max(h, 1U)), w2(max(w>>1, 1U)), h2(max(h>>1, 1U));
		unsigned const xinc((w2 < w1) ? ncolors : 0), yinc((h2 < h1) ? ncolors*w1 : 0);
		uint8_t const *const idata((level == 1) ? data : idatav.data());
		odata.resize(ncolors*w2*h2);

#pragma omp parallel for schedule(static)
		for (int y = 0; y < (int)h2; ++y) {
			for (unsigned x = 0; x < w2; ++x) {
				unsigned const ix1(ncolors*(y*w2+x)), ix2(ncolors*((y<<1)*w1+(x<<1)));

				if (ncolors == 1) {
					odata[ix1] = uint8_t(((unsigned)idata[ix2] + idata[ix2+xinc] + idata[ix2+yinc] + idata[ix2+yinc+xinc]) >> 2);
				}
				else if (ncolors == 3) {
					UNROLL_3X(odata[ix1+i_] = uint8_t(((unsigned)idata[ix2+i_] + idata[ix2+xinc+i_] + idata[ix2+yinc+i_] + idata[ix2+yinc+xinc+i_]) >> 2);)
				}
				else { // custom alpha mipmaps
					assert(ncolors == 4);
					unsigned const a1(idata[ix2+3]), a2(idata[ix2+xinc+3]), a3(idata[ix2+yinc+3]), a4(idata[ix2+yinc+xinc+3]);
					unsigned const a_sum(a1 + a2 + a3 + a4);

					if (a_sum == 0) { // fully transparent
						if (use_mipmaps == 4) {UNROLL_3X(odata[ix1+i_] = cw.c[i_];)} // use average texture color
						else { // color is average of all 4 values
							UNROLL_3X(odata[ix1+i_] = uint8_t(((unsigned)idata[ix2+i_] + idata[ix2+xinc+i_] + idata[ix2+yinc+i_] + idata[ix2+yinc+xinc+i_]) / 4);)
						}
						odata[ix1+3] = 0;
					}
					else { // pre-multiplied and normalized colors
						if (use_mipmaps == 4) {
							unsigned const a_cw(1020 - a_sum); // use average texture color for transparent pixels
							UNROLL_3X(odata[ix1+i_] = uint8_t((a1*idata[ix2+i_] + a2*idata[ix2+xinc+i_] + a3*idata[ix2+yinc+i_] + a4*idata[ix2+yinc+xinc+i_] + a_cw*cw.c[i_]) / 1020);)
						}
						else {
							UNROLL_3X(odata[ix1+i_] = uint8_t((a1*idata[ix2+i_] + a2*idata[ix2+xinc+i_] + a3*idata[ix2+yinc+i_] + a4*idata[ix2+yinc+xinc+i_]) / a_sum);)
						}
						odata[ix1+3] = min(255U, min(max(max(a1, a2), max(a3, a4)), unsigned(mipmap_alpha_weight*a_sum)));
					}
				}
			} // for x
		} // for y
		if (is_texture_compressed()) { // uses stb
			dxt_texture_compress(odata.data(), comp_data, w2, h2, ncolors);
			GL_CHECK(glCompressedTexImage2D(GL_TEXTURE_2D, level, calc_internal_format(), w2, h2, 0, comp_data.size(), comp_data.data());)
		}
		else {
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // needed for mipmap levels where width*ncolors is not aligned
			glTexImage2D(GL_TEXTURE_2D, level, calc_internal_format(), w2, h2, 0, calc_format(), get_data_format(), odata.data());
			glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
		}
		idatav.swap(odata);
	} // for level
}


