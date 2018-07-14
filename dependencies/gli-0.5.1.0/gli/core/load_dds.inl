///////////////////////////////////////////////////////////////////////////////////
/// OpenGL Image (gli.g-truc.net)
///
/// Copyright (c) 2008 - 2013 G-Truc Creation (www.g-truc.net)
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
/// THE SOFTWARE.
///
/// @ref core
/// @file gli/core/load_dds.inl
/// @date 2010-09-26 / 2013-01-28
/// @author Christophe Riccio
///////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <cassert>

namespace gli{
namespace detail
{
	// DDS Documentation
	/*
		http://msdn.microsoft.com/en-us/library/bb943991(VS.85).aspx#File_Layout1
		http://msdn.microsoft.com/en-us/library/bb943992.aspx
	*/

	enum ddsCubemapflag
	{
		DDSCAPS2_CUBEMAP				= 0x00000200,
		DDSCAPS2_CUBEMAP_POSITIVEX		= 0x00000400,
		DDSCAPS2_CUBEMAP_NEGATIVEX		= 0x00000800,
		DDSCAPS2_CUBEMAP_POSITIVEY		= 0x00001000,
		DDSCAPS2_CUBEMAP_NEGATIVEY		= 0x00002000,
		DDSCAPS2_CUBEMAP_POSITIVEZ		= 0x00004000,
		DDSCAPS2_CUBEMAP_NEGATIVEZ		= 0x00008000,
		DDSCAPS2_VOLUME					= 0x00200000
	};

	glm::uint32 const DDSCAPS2_CUBEMAP_ALLFACES = (
		DDSCAPS2_CUBEMAP_POSITIVEX | DDSCAPS2_CUBEMAP_NEGATIVEX |
		DDSCAPS2_CUBEMAP_POSITIVEY | DDSCAPS2_CUBEMAP_NEGATIVEY |
		DDSCAPS2_CUBEMAP_POSITIVEZ | DDSCAPS2_CUBEMAP_NEGATIVEZ);

	enum ddsFlag
	{
		DDSD_CAPS			= 0x00000001,
		DDSD_HEIGHT			= 0x00000002,
		DDSD_WIDTH			= 0x00000004,
		DDSD_PITCH			= 0x00000008,
		DDSD_PIXELFORMAT	= 0x00001000,
		DDSD_MIPMAPCOUNT	= 0x00020000,
		DDSD_LINEARSIZE		= 0x00080000,
		DDSD_DEPTH			= 0x00800000
	};

	enum ddsSurfaceflag
	{
		DDSCAPS_COMPLEX				= 0x00000008,
		DDSCAPS_MIPMAP				= 0x00400000,
		DDSCAPS_TEXTURE				= 0x00001000
	};

	struct ddsPixelFormat
	{
		glm::uint32 size; // 32
		glm::uint32 flags;
		glm::uint32 fourCC;
		glm::uint32 bpp;
		glm::uint32 redMask;
		glm::uint32 greenMask;
		glm::uint32 blueMask;
		glm::uint32 alphaMask;
	};

	struct ddsHeader
	{
		glm::uint32 size;
		glm::uint32 flags;
		glm::uint32 height;
		glm::uint32 width;
		glm::uint32 pitch;
		glm::uint32 depth;
		glm::uint32 mipMapLevels;
		glm::uint32 reserved1[11];
		ddsPixelFormat format;
		glm::uint32 surfaceFlags;
		glm::uint32 cubemapFlags;
		glm::uint32 reserved2[3];
	};

	enum D3D10_RESOURCE_DIMENSION 
	{
		D3D10_RESOURCE_DIMENSION_UNKNOWN     = 0,
		D3D10_RESOURCE_DIMENSION_BUFFER      = 1,
		D3D10_RESOURCE_DIMENSION_TEXTURE1D   = 2,
		D3D10_RESOURCE_DIMENSION_TEXTURE2D   = 3,
		D3D10_RESOURCE_DIMENSION_TEXTURE3D   = 4 
	};

	enum D3D10_RESOURCE_MISC_FLAG 
	{
		D3D10_RESOURCE_MISC_GENERATE_MIPS       = 0x1L,
		D3D10_RESOURCE_MISC_SHARED              = 0x2L,
		D3D10_RESOURCE_MISC_TEXTURECUBE         = 0x4L,
		D3D10_RESOURCE_MISC_SHARED_KEYEDMUTEX   = 0x10L,
		D3D10_RESOURCE_MISC_GDI_COMPATIBLE      = 0x20L 
	};

	struct ddsHeader10
	{
		ddsHeader10() :
			Format(DXGI_FORMAT_UNKNOWN),
			resourceDimension(D3D10_RESOURCE_DIMENSION_UNKNOWN),
			miscFlag(0),
			arraySize(1),
			reserved(0)
		{}

		dxgiFormat					Format;
		D3D10_RESOURCE_DIMENSION	resourceDimension;
		glm::uint32					miscFlag; // D3D10_RESOURCE_MISC_GENERATE_MIPS
		glm::uint32					arraySize;
		glm::uint32					reserved;
	};

	inline gli::format format_fourcc2gli_cast(glm::uint32 const & Flags, glm::uint32 const & FourCC)
	{
		switch(FourCC)
		{
		case D3DFMT_DXT1:
			return Flags & DDPF_ALPHAPIXELS ? RGBA_DXT1 : RGB_DXT1;
		case D3DFMT_DXT2:
		case D3DFMT_DXT3:
			return RGBA_DXT3;
		case D3DFMT_DXT4:
		case D3DFMT_DXT5:
			return RGBA_DXT5;
		case D3DFMT_ATI1:
		case D3DFMT_AT1N:
			return R_ATI1N_UNORM;
		case D3DFMT_ATI2:
		case D3DFMT_AT2N:
			return RG_ATI2N_UNORM;
		case D3DFMT_R16F:
			return R16F;
		case D3DFMT_G16R16F:
			return RG16F;
		case D3DFMT_A16B16G16R16F:
			return RGBA16F;
		case D3DFMT_R32F:
			return R32F;
		case D3DFMT_G32R32F:
			return RG32F;
		case D3DFMT_A32B32G32R32F:
			return RGBA32F;
		case D3DFMT_R8G8B8:
			return RGB8U;
		case D3DFMT_A8R8G8B8:
		case D3DFMT_X8R8G8B8:
		case D3DFMT_A8B8G8R8:
		case D3DFMT_X8B8G8R8:
			return RGBA8U;
		case D3DFMT_R5G6B5:
			return R5G6B5;
		case D3DFMT_A4R4G4B4:
		case D3DFMT_X4R4G4B4:
			return RGBA4;
		case D3DFMT_G16R16:
			return RG16U;
		case D3DFMT_A16B16G16R16:
			return RGBA16U;
		case D3DFMT_A2R10G10B10:
		case D3DFMT_A2B10G10R10:
			return RGB10A2;
		default:
			assert(0);
			return FORMAT_NULL;
		}
	}

	inline gli::format format_dds2gli_cast(dxgiFormat const & Format)
	{
		static gli::format const Cast[] = 
		{
			gli::FORMAT_NULL,	//DXGI_FORMAT_UNKNOWN                      = 0,
			gli::RGBA32U,		//DXGI_FORMAT_R32G32B32A32_TYPELESS        = 1,
			gli::RGBA32F,		//DXGI_FORMAT_R32G32B32A32_FLOAT           = 2,
			gli::RGBA32U,		//DXGI_FORMAT_R32G32B32A32_UINT            = 3,
			gli::RGBA32I,		//DXGI_FORMAT_R32G32B32A32_SINT            = 4,
			gli::RGB32U,			//DXGI_FORMAT_R32G32B32_TYPELESS           = 5,
			gli::RGB32F,			//DXGI_FORMAT_R32G32B32_FLOAT              = 6,
			gli::RGB32U,			//DXGI_FORMAT_R32G32B32_UINT               = 7,
			gli::RGB32I,			//DXGI_FORMAT_R32G32B32_SINT               = 8,
			gli::RGBA16U,		//DXGI_FORMAT_R16G16B16A16_TYPELESS        = 9,
			gli::RGBA16F,		//DXGI_FORMAT_R16G16B16A16_FLOAT           = 10,
			gli::RGBA16U,		//DXGI_FORMAT_R16G16B16A16_UNORM           = 11,
			gli::RGBA16I,		//DXGI_FORMAT_R16G16B16A16_UINT            = 12,
			gli::RGBA16I,		//DXGI_FORMAT_R16G16B16A16_SNORM           = 13,
			gli::RGBA16I,		//DXGI_FORMAT_R16G16B16A16_SINT            = 14,
			gli::RG32U,			//DXGI_FORMAT_R32G32_TYPELESS              = 15,
			gli::RG32F,			//DXGI_FORMAT_R32G32_FLOAT                 = 16,
			gli::RG32U,			//DXGI_FORMAT_R32G32_UINT                  = 17,
			gli::RG32I,			//DXGI_FORMAT_R32G32_SINT                  = 18,
			gli::FORMAT_NULL,	//DXGI_FORMAT_R32G8X24_TYPELESS            = 19,
			gli::D32FS8X24,		//DXGI_FORMAT_D32_FLOAT_S8X24_UINT         = 20,
			gli::FORMAT_NULL,	//DXGI_FORMAT_R32_FLOAT_X8X24_TYPELESS     = 21,
			gli::FORMAT_NULL,	//DXGI_FORMAT_X32_TYPELESS_G8X24_UINT      = 22,
			gli::RGB10A2,		//DXGI_FORMAT_R10G10B10A2_TYPELESS         = 23,
			gli::RGB10A2,		//DXGI_FORMAT_R10G10B10A2_UNORM            = 24,
			gli::RGB10A2,		//DXGI_FORMAT_R10G10B10A2_UINT             = 25,
			gli::RG11B10F,		//DXGI_FORMAT_R11G11B10_FLOAT              = 26,
			gli::RGBA8U,			//DXGI_FORMAT_R8G8B8A8_TYPELESS            = 27,
			gli::RGBA8U,			//DXGI_FORMAT_R8G8B8A8_UNORM               = 28,
			gli::RGBA8U,			//DXGI_FORMAT_R8G8B8A8_UNORM_SRGB          = 29,
			gli::RGBA8U,			//DXGI_FORMAT_R8G8B8A8_UINT                = 30,
			gli::RGBA8I,			//DXGI_FORMAT_R8G8B8A8_SNORM               = 31,
			gli::RGBA8I,			//DXGI_FORMAT_R8G8B8A8_SINT                = 32,
			gli::RG16U,			//DXGI_FORMAT_R16G16_TYPELESS              = 33,
			gli::RG16F,			//DXGI_FORMAT_R16G16_FLOAT                 = 34,
			gli::RG16U,			//DXGI_FORMAT_R16G16_UNORM                 = 35,
			gli::RG16U,			//DXGI_FORMAT_R16G16_UINT                  = 36,
			gli::RG16I,			//DXGI_FORMAT_R16G16_SNORM                 = 37,
			gli::RG16I,			//DXGI_FORMAT_R16G16_SINT                  = 38,
			gli::R32F,			//DXGI_FORMAT_R32_TYPELESS                 = 39,
			gli::D32F,			//DXGI_FORMAT_D32_FLOAT                    = 40,
			gli::R32F,			//DXGI_FORMAT_R32_FLOAT                    = 41,
			gli::R32U,			//DXGI_FORMAT_R32_UINT                     = 42,
			gli::R32I,			//DXGI_FORMAT_R32_SINT                     = 43,
			gli::FORMAT_NULL,	//DXGI_FORMAT_R24G8_TYPELESS               = 44,
			gli::FORMAT_NULL,	//DXGI_FORMAT_D24_UNORM_S8_UINT            = 45,
			gli::FORMAT_NULL,	//DXGI_FORMAT_R24_UNORM_X8_TYPELESS        = 46,
			gli::FORMAT_NULL,	//DXGI_FORMAT_X24_TYPELESS_G8_UINT         = 47,
			gli::RG8U,			//DXGI_FORMAT_R8G8_TYPELESS                = 48,
			gli::RG8U,			//DXGI_FORMAT_R8G8_UNORM                   = 49,
			gli::RG8U,			//DXGI_FORMAT_R8G8_UINT                    = 50,
			gli::RG8I,			//DXGI_FORMAT_R8G8_SNORM                   = 51,
			gli::RG8I,			//DXGI_FORMAT_R8G8_SINT                    = 52,
			gli::R16U,			//DXGI_FORMAT_R16_TYPELESS                 = 53,
			gli::R16F,			//DXGI_FORMAT_R16_FLOAT                    = 54,
			gli::D16,			//DXGI_FORMAT_D16_UNORM                    = 55,
			gli::R16U,			//DXGI_FORMAT_R16_UNORM                    = 56,
			gli::R16U,			//DXGI_FORMAT_R16_UINT                     = 57,
			gli::R16I,			//DXGI_FORMAT_R16_SNORM                    = 58,
			gli::R16I,			//DXGI_FORMAT_R16_SINT                     = 59,
			gli::R8U,			//DXGI_FORMAT_R8_TYPELESS                  = 60,
			gli::R8U,			//DXGI_FORMAT_R8_UNORM                     = 61,
			gli::R8U,			//DXGI_FORMAT_R8_UINT                      = 62,
			gli::R8I,			//DXGI_FORMAT_R8_SNORM                     = 63,
			gli::R8I,			//DXGI_FORMAT_R8_SINT                      = 64,
			gli::R8U,			//DXGI_FORMAT_A8_UNORM                     = 65,
			gli::FORMAT_NULL,	//DXGI_FORMAT_R1_UNORM                     = 66,
			gli::RGB9E5,			//DXGI_FORMAT_R9G9B9E5_SHAREDEXP           = 67,
			gli::FORMAT_NULL,		//DXGI_FORMAT_R8G8_B8G8_UNORM              = 68,
			gli::FORMAT_NULL,		//DXGI_FORMAT_G8R8_G8B8_UNORM              = 69,
			gli::RGBA_DXT1,			//DXGI_FORMAT_BC1_TYPELESS                 = 70,
			gli::RGBA_DXT1,			//DXGI_FORMAT_BC1_UNORM                    = 71,
			gli::RGBA_DXT1,			//DXGI_FORMAT_BC1_UNORM_SRGB               = 72,
			gli::RGBA_DXT3,			//DXGI_FORMAT_BC2_TYPELESS                 = 73,
			gli::RGBA_DXT3,			//DXGI_FORMAT_BC2_UNORM                    = 74,
			gli::RGBA_DXT3,			//DXGI_FORMAT_BC2_UNORM_SRGB               = 75,
			gli::RGBA_DXT5,			//DXGI_FORMAT_BC3_TYPELESS                 = 76,
			gli::RGBA_DXT5,			//DXGI_FORMAT_BC3_UNORM                    = 77,
			gli::RGBA_DXT5,			//DXGI_FORMAT_BC3_UNORM_SRGB               = 78,
			gli::R_ATI1N_UNORM,		//DXGI_FORMAT_BC4_TYPELESS                 = 79,
			gli::R_ATI1N_UNORM,		//DXGI_FORMAT_BC4_UNORM                    = 80,
			gli::R_ATI1N_SNORM,		//DXGI_FORMAT_BC4_SNORM                    = 81,
			gli::RG_ATI2N_UNORM,	//DXGI_FORMAT_BC5_TYPELESS                 = 82,
			gli::RG_ATI2N_UNORM,	//DXGI_FORMAT_BC5_UNORM                    = 83,
			gli::RG_ATI2N_SNORM,	//DXGI_FORMAT_BC5_SNORM                    = 84,
			gli::FORMAT_NULL,		//DXGI_FORMAT_B5G6R5_UNORM                 = 85,
			gli::FORMAT_NULL,		//DXGI_FORMAT_B5G5R5A1_UNORM               = 86,
			gli::RGBA8_UNORM,			//DXGI_FORMAT_B8G8R8A8_UNORM               = 87,
			gli::RGB8_UNORM,				//DXGI_FORMAT_B8G8R8X8_UNORM               = 88,
			gli::FORMAT_NULL,		//DXGI_FORMAT_R10G10B10_XR_BIAS_A2_UNORM   = 89,
			gli::RGBA8_UNORM,			//DXGI_FORMAT_B8G8R8A8_TYPELESS            = 90,
			gli::RGBA8_UNORM,			//DXGI_FORMAT_B8G8R8A8_UNORM_SRGB          = 91,
			gli::RGB8_UNORM,				//DXGI_FORMAT_B8G8R8X8_TYPELESS            = 92,
			gli::SRGB8,						//DXGI_FORMAT_B8G8R8X8_UNORM_SRGB          = 93,
			gli::RGB_BP_UNSIGNED_FLOAT,		//DXGI_FORMAT_BC6H_TYPELESS                = 94,
			gli::RGB_BP_UNSIGNED_FLOAT,		//DXGI_FORMAT_BC6H_UF16                    = 95,
			gli::RGB_BP_SIGNED_FLOAT,		//DXGI_FORMAT_BC6H_SF16                    = 96,
			gli::RGB_BP_UNORM,				//DXGI_FORMAT_BC7_TYPELESS                 = 97,
			gli::RGB_BP_UNORM,				//DXGI_FORMAT_BC7_UNORM                    = 98,
			gli::RGB_BP_UNORM,				//DXGI_FORMAT_BC7_UNORM_SRGB               = 99,
			gli::R32U						//DXGI_FORMAT_FORCE_UINT                   = 0xffffffffUL 
		};

		return Cast[Format];
	}
}//namespace detail

inline storage load_dds
(
	char const * Filename
)
{
	std::ifstream FileIn(Filename, std::ios::in | std::ios::binary);
	assert(!FileIn.fail());

	if(FileIn.fail())
		return storage();

	detail::ddsHeader HeaderDesc;
	detail::ddsHeader10 HeaderDesc10;
	char Magic[4]; 

	//* Read magic number and check if valid .dds file 
	FileIn.read((char*)&Magic, sizeof(Magic));

	assert(strncmp(Magic, "DDS ", 4) == 0);

	// Get the surface descriptor 
	FileIn.read((char*)&HeaderDesc, sizeof(HeaderDesc));
	if(HeaderDesc.format.flags & detail::DDPF_FOURCC && HeaderDesc.format.fourCC == detail::D3DFMT_DX10)
		FileIn.read((char*)&HeaderDesc10, sizeof(HeaderDesc10));

	gli::format Format(gli::FORMAT_NULL);
	if(HeaderDesc.format.fourCC == detail::D3DFMT_DX10)
		Format = detail::format_dds2gli_cast(HeaderDesc10.Format);
	else if(HeaderDesc.format.flags & detail::DDPF_FOURCC)
		Format = detail::format_fourcc2gli_cast(HeaderDesc.format.flags, HeaderDesc.format.fourCC);
	else if(HeaderDesc.format.flags & detail::DDPF_RGB)
	{
		switch(HeaderDesc.format.bpp)
		{
		case 8:
			Format = R8_UNORM;
			break;
		case 16:
			Format = RG8_UNORM;
			break;
		case 24:
			Format = RGB8_UNORM;
			break;
		case 32:
			Format = RGBA8_UNORM;
			break;
		}
	}
	else
		assert(0);

	std::streamoff Curr = FileIn.tellg();
	FileIn.seekg(0, std::ios_base::end);
	std::streamoff End = FileIn.tellg();
	FileIn.seekg(Curr, std::ios_base::beg);

	storage::size_type const MipMapCount = (HeaderDesc.flags & detail::DDSD_MIPMAPCOUNT) ? 
		HeaderDesc.mipMapLevels : 1;

	storage::size_type FaceCount(1);
	if(HeaderDesc.cubemapFlags & detail::DDSCAPS2_CUBEMAP)
		FaceCount = int(glm::bitCount(HeaderDesc.cubemapFlags & detail::DDSCAPS2_CUBEMAP_ALLFACES));

	storage::size_type DepthCount = 1;
	if(HeaderDesc.cubemapFlags & detail::DDSCAPS2_VOLUME)
			DepthCount = HeaderDesc.depth;

	storage Storage(
		HeaderDesc10.arraySize, 
		FaceCount,
		MipMapCount,
		Format,
		storage::dimensions_type(HeaderDesc.width, HeaderDesc.height, DepthCount));

	FileIn.read((char*)Storage.data(), std::size_t(End - Curr));

	return Storage;
}

}//namespace gli
