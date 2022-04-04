// 3D World - Texture Compression and Mipmap Utility Functions
// by Frank Gennari
// 4/3/22
#include "3DWorld.h"
#include "function_registry.h"


void texture_t::create_custom_mipmaps() {

	assert(is_allocated());
	GLenum const format(calc_format());
	unsigned const tsize(num_bytes());
	vector<unsigned char> idata, odata;
	idata.resize(tsize);
	memcpy(&idata.front(), data, tsize);
	color_wrapper cw; cw.set_c4(color);

	for (unsigned w = width, h = height, level = 1; w > 1 || h > 1; w >>= 1, h >>= 1, ++level) {
		unsigned const w1(max(w,    1U)), h1(max(h,    1U));
		unsigned const w2(max(w>>1, 1U)), h2(max(h>>1, 1U));
		unsigned const xinc((w2 < w1) ? ncolors : 0), yinc((h2 < h1) ? ncolors*w1 : 0);
		odata.resize(ncolors*w2*h2);

		for (unsigned y = 0; y < h2; ++y) {
			for (unsigned x = 0; x < w2; ++x) {
				unsigned const ix1(ncolors*(y*w2+x)), ix2(ncolors*((y<<1)*w1+(x<<1)));

				if (ncolors == 1) {
					odata[ix1] = (unsigned char)(((unsigned)idata[ix2] + idata[ix2+xinc] + idata[ix2+yinc] + idata[ix2+yinc+xinc]) >> 2);
				}
				else if (ncolors == 3) {
					UNROLL_3X(odata[ix1+i_] = (unsigned char)(((unsigned)idata[ix2+i_] + idata[ix2+xinc+i_] + idata[ix2+yinc+i_] + idata[ix2+yinc+xinc+i_]) >> 2);)
				}
				else { // custom alpha mipmaps
					assert(ncolors == 4);
					unsigned const a1(idata[ix2+3]), a2(idata[ix2+xinc+3]), a3(idata[ix2+yinc+3]), a4(idata[ix2+yinc+xinc+3]);
					unsigned const a_sum(a1 + a2 + a3 + a4);

					if (a_sum == 0) { // fully transparent
						if (use_mipmaps == 4) {UNROLL_3X(odata[ix1+i_] = cw.c[i_];)} // use average texture color
						else { // color is average of all 4 values
							UNROLL_3X(odata[ix1+i_] = (unsigned char)(((unsigned)idata[ix2+i_] + idata[ix2+xinc+i_] + idata[ix2+yinc+i_] + idata[ix2+yinc+xinc+i_]) / 4);)
						}
						odata[ix1+3] = 0;
					}
					else { // pre-multiplied and normalized colors
						if (use_mipmaps == 4) {
							unsigned const a_cw(1020 - a_sum); // use average texture color for transparent pixels
							UNROLL_3X(odata[ix1+i_] = (unsigned char)((a1*idata[ix2+i_] + a2*idata[ix2+xinc+i_] + a3*idata[ix2+yinc+i_] + a4*idata[ix2+yinc+xinc+i_] + a_cw*cw.c[i_]) / 1020);)
						}
						else {
							UNROLL_3X(odata[ix1+i_] = (unsigned char)((a1*idata[ix2+i_] + a2*idata[ix2+xinc+i_] + a3*idata[ix2+yinc+i_] + a4*idata[ix2+yinc+xinc+i_]) / a_sum);)
						}
						odata[ix1+3] = min(255U, min(max(max(a1, a2), max(a3, a4)), unsigned(mipmap_alpha_weight*a_sum)));
					}
				}
			} // for x
		} // for y
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // needed for mipmap levels where width*ncolors is not aligned
		glTexImage2D(GL_TEXTURE_2D, level, calc_internal_format(), w2, h2, 0, format, get_data_format(), &odata.front());
		glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
		idata.swap(odata);
	} // for w
#if 0
	for (unsigned level = 0; level < mm_offsets.size(); ++level) {
		unsigned const tsz(width >> level);
		assert(tsz > 1);
		unsigned char const *const src((level == 0) ? data : (mm_data + mm_offsets[level-1]));
		int const ret(gluScaleImage(format, tsz,   tsz,   get_data_format(), src,
			tsz/2, tsz/2, get_data_format(), (mm_data + mm_offsets[level])));
		if (ret) cout << "GLU error during mipmap image scale: " << gluErrorString(ret) << "." << endl;
	}
#endif
}


