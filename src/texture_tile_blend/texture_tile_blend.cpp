// 3D World - Texture Tiling and Blending from Deliot2019
// edited by Frank Gennari
// https://eheitzresearch.wordpress.com/738-2/
// 2/16/19
#include "tlingandblending.h"
#include "../function_registry.h"
#include "../gl_ext_arb.h"
#include "../shaders.h"

extern unsigned context_clear_count;

GLuint CreateGLTextureFromTextureDataStruct(const TextureDataFloat& im, GLenum wrapMode, bool generateMips) {

	if (im.data.empty()) return 0;
	GLuint texture;
	glGenTextures(1, &texture);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, generateMips ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMode);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMode);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, im.width, im.height, 0, GL_RGB, GL_FLOAT, im.data.data());
	if (generateMips) {glGenerateMipmap(GL_TEXTURE_2D);}
	return texture;
}

void tile_blend_tex_data_t::bind_shader(shader_t &s) const {

	assert(textures_valid());
	unsigned const lut_tu_id = 21; // use a value that won't conflict with other
	bind_texture_tu(tid_tinput, 0);
	bind_texture_tu(tid_lut,    lut_tu_id);
	s.add_uniform_int("Tinput", 0);
	s.add_uniform_int("invT",   lut_tu_id);
	s.add_uniform_vector3d("_colorSpaceVector1", colorSpaceVector1);
	s.add_uniform_vector3d("_colorSpaceVector2", colorSpaceVector2);
	s.add_uniform_vector3d("_colorSpaceVector3", colorSpaceVector3);
	s.add_uniform_vector3d("_colorSpaceOrigin",  colorSpaceOrigin );
}

void tile_blend_tex_data_t::create_textures(texture_t const &texture) {

	timer_t timer("Create Tile Blend Textures");
	assert(texture.ncolors == 3); // only supports RGB for now
	unsigned const num_bytes(texture.num_bytes());
	TextureDataFloat input(texture.width, texture.height, texture.ncolors);
	for (unsigned i = 0; i < num_bytes; ++i) {input.data[i] = texture.get_data()[i] / 255.0;} // convert unsigned char to bytes
	TextureDataFloat Tinput;
	TextureDataFloat lut;
	Precomputations(input, Tinput, lut, colorSpaceVector1, colorSpaceVector2, colorSpaceVector3, colorSpaceOrigin);
	tid_tinput = CreateGLTextureFromTextureDataStruct(Tinput, GL_REPEAT, true);
	tid_lut    = CreateGLTextureFromTextureDataStruct(lut, GL_CLAMP_TO_EDGE, false);
}

void tile_blend_tex_data_t::ensure_textures(unsigned tid) {
	if (context_count != context_clear_count) {
		clear_context();
		context_count = context_clear_count;
	}
	if (!textures_valid()) {create_textures(get_texture_by_id(tid));}
}

void tile_blend_tex_data_t::clear_context() {
	free_texture(tid_tinput);
	free_texture(tid_lut);
}

