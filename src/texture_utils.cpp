// 3D World - Texture Compression and Mipmap Utility Functions
// by Frank Gennari
// 4/3/22
#include "3DWorld.h"
#include "function_registry.h"
#include "profiler.h"

#define STB_DXT_IMPLEMENTATION
#include "stb_dxt.h"
bool const USE_STB_DXT = 1;

bool gen_mipmaps(unsigned dim=2);


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

void create_one_mipmap(uint8_t const *const idata, vector<uint8_t> &odata, unsigned w1, unsigned h1, unsigned w2, unsigned h2,
	unsigned ncolors, unsigned use_mipmaps, colorRGBA const &color, float mipmap_alpha_weight)
{
	unsigned const xinc((w2 < w1) ? ncolors : 0), yinc((h2 < h1) ? ncolors*w1 : 0);
	odata.resize(ncolors*w2*h2);

	if (use_mipmaps == 3 || use_mipmaps == 4) { // custom mipmap path
		color_wrapper cw(color); // for use_mipmaps == 4 with RGBA

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
	}
	else { // simple 2x2 box filter path
#pragma omp parallel for schedule(static)
		for (int y = 0; y < (int)h2; ++y) { // simple 2x2 box filter
			for (int x = 0; x < (int)w2; ++x) {
				unsigned const ix1(ncolors*(y*w2+x)), ix2(ncolors*((y<<1)*w1+(x<<1)));

				for (int n = 0; n < ncolors; ++n) {
					odata[ix1+n] = uint8_t(((unsigned)idata[ix2+n] + idata[ix2+xinc+n] + idata[ix2+yinc+n] + idata[ix2+yinc+xinc+n]) >> 2);
				}
			}
		} // for y
	}
}

void texture_t::compress_and_send_texture_with_mipmaps() {
	assert(is_allocated());
	assert(width > 0 && height > 0);
	bool const compressed(is_texture_compressed()), use_custom_compress(USE_STB_DXT && compressed && (ncolors == 3 || ncolors == 4));
	bool const use_normal_mipmaps(use_mipmaps == 1 || use_mipmaps == 2), use_custom_mipmaps(use_mipmaps == 3 || use_mipmaps == 4);
	vector<uint8_t> comp_data, idatav, odata; // reused across calls; doesn't seem to help much

	if (use_custom_compress) { // compressed RGB or RGBA + mipmaps
		//highres_timer_t timer("compress_and_send_texture", 1, 1); // enabled, no loading screen; 1850ms / 1000ms with PBO (but total load time is not actually faster)
		if (use_normal_mipmaps) {
			//highres_timer_t timer("create_compressed_mipmaps", 1, 1); // enabled, no loading screen; 1350ms total for city + cars + people

			for (unsigned w = width, h = height, level = 1; w > 1 || h > 1; w >>= 1, h >>= 1, ++level) {
				unsigned const w1(max(w, 1U)), h1(max(h, 1U)), w2(max(w>>1, 1U)), h2(max(h>>1, 1U));
				uint8_t const *const idata((level == 1) ? data : idatav.data());
				create_one_mipmap(idata, odata, w1, h1, w2, h2, ncolors, use_mipmaps, color, mipmap_alpha_weight);
				dxt_texture_compress(odata.data(), comp_data, w2, h2, ncolors);
				GL_CHECK(glCompressedTexImage2D(GL_TEXTURE_2D, level, calc_internal_format(), w2, h2, 0, comp_data.size(), comp_data.data()););
				idatav.swap(odata);
			} // for level
		}
		// sending the main texture last seems to be slightly faster because it will block on the first mipmap if sent first; if sent last it may overlap with compress
		dxt_texture_compress(data, comp_data, width, height, ncolors); // 640ms
		GL_CHECK(glCompressedTexImage2D(GL_TEXTURE_2D, 0, calc_internal_format(), width, height, 0, comp_data.size(), comp_data.data());); // 36ms
	}
	else { // font atlas and noise gen texture
		glTexImage2D(GL_TEXTURE_2D, 0, calc_internal_format(), width, height, 0, calc_format(), get_data_format(), data); // 44ms
		if (use_normal_mipmaps) {gen_mipmaps();} // Note: compressed mipmaps are created inside compress_and_send_texture()
	}
	if (use_custom_mipmaps) { // custom mipmap creation for compressed and non-compressed textures
		//highres_timer_t timer("Create Custom Mipmaps"); // 350ms total for city + cars + people

		for (unsigned w = width, h = height, level = 1; w > 1 || h > 1; w >>= 1, h >>= 1, ++level) {
			unsigned const w1(max(w, 1U)), h1(max(h, 1U)), w2(max(w>>1, 1U)), h2(max(h>>1, 1U));
			uint8_t const *const idata((level == 1) ? data : idatav.data());
			create_one_mipmap(idata, odata, w1, h1, w2, h2, ncolors, use_mipmaps, color, mipmap_alpha_weight);

			if (compressed) { // uses stb
				dxt_texture_compress(odata.data(), comp_data, w2, h2, ncolors);
				GL_CHECK(glCompressedTexImage2D(GL_TEXTURE_2D, level, calc_internal_format(), w2, h2, 0, comp_data.size(), comp_data.data());) // 140ms
			}
			else { // 10ms
				glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // needed for mipmap levels where width*ncolors is not aligned
				glTexImage2D(GL_TEXTURE_2D, level, calc_internal_format(), w2, h2, 0, calc_format(), get_data_format(), odata.data());
				glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
			}
			idatav.swap(odata);
		} // for level
	}
}

