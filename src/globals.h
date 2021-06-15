// 3D World - Global Defines
// by Frank Gennari
// 12/12/10
#pragma once

// fixes for MS Visual Studio 2005 performance problems
#define _HAS_ITERATOR_DEBUGGING 0
// MSVS 2017 and above
#define _ITERATOR_DEBUG_LEVEL 0
//#define _SECURE_SCL 0

// required for 64-bit builds to avoid size_t to unsigned conversion warnings
#ifdef _WIN64
#pragma warning (disable : 4267)
#pragma warning (disable : 4244)
#endif

