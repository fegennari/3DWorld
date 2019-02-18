/*****************************************************************************/
/*****************************************************************************/
/********************************* Includes **********************************/
/*****************************************************************************/
/*****************************************************************************/
#include "../3DWorld.h"
#include "jacobi.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cfloat> // for FLT_MAX
using namespace std;

/*****************************************************************************/
/*****************************************************************************/
/******************************** Parameters *********************************/
/*****************************************************************************/
/*****************************************************************************/

#define USE_DXT_COMPRESSION false // Use DXT1 (true) or GL_RGB8 (false) (Section 1.6)
#define GAUSSIAN_AVERAGE 0.5f // Expectation of the Gaussian distribution
#define GAUSSIAN_STD 0.16666f // Std of the Gaussian distribution
#define LUT_WIDTH 128 // Size of the look-up table

/*****************************************************************************/
/*****************************************************************************/
/*************************** Types and Structures ****************************/
/*****************************************************************************/
/*****************************************************************************/

typedef vector3d vec3;
typedef vector2d vec2;

struct TextureDataFloat
{
	TextureDataFloat() : data(), width(0), height(0), channels(0) {}
	TextureDataFloat(const int w, const int h, const int c) :
		data(w * h * c), width(w), height(h), channels(c)
	{
	}
	float GetPixel(int w, int h, int c)
	{
		return data[h * width * channels + w * channels + c];
	}
	vec3 GetColorAt(int w, int h)
	{
		return vec3(
			data[h * width * channels + w * channels + 0],
			data[h * width * channels + w * channels + 1],
			data[h * width * channels + w * channels + 2]);
	}
	void SetPixel(int w, int h, int c, float value)
	{
		data[h * width * channels + w * channels + c] = value;
	}
	void SetColorAt(int w, int h, vec3 const &value)
	{
		data[h * width * channels + w * channels + 0] = value.x;
		data[h * width * channels + w * channels + 1] = value.y;
		data[h * width * channels + w * channels + 2] = value.z;
	}
	vector<float> data;
	int width;
	int height;
	int channels;
};

struct PixelSortStruct
{
	int x;
	int y;
	float value;
	bool operator < (const PixelSortStruct& other) const
	{
		return (value < other.value);
	}
};


/*****************************************************************************/
/*****************************************************************************/
/**************** Section 1.3.1 Target Gaussian distribution *****************/
/*****************************************************************************/
/*****************************************************************************/

float Erf(float x)
{
	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = abs(x);

	// A&S formula 7.1.26
	float t = 1.0f / (1.0f + 0.3275911f * x);
	float y = 1.0f - (((((1.061405429f * t + -1.453152027f) * t) + 1.421413741f)
		* t + -0.284496736f) * t + 0.254829592f) * t * exp(-x * x);

	return sign * y;
}

float ErfInv(float x)
{
	float w, p;
	w = -log((1.0f - x) * (1.0f + x));
	if (w < 5.000000f)
	{
		w = w - 2.500000f;
		p = 2.81022636e-08f;
		p = 3.43273939e-07f + p * w;
		p = -3.5233877e-06f + p * w;
		p = -4.39150654e-06f + p * w;
		p = 0.00021858087f + p * w;
		p = -0.00125372503f + p * w;
		p = -0.00417768164f + p * w;
		p = 0.246640727f + p * w;
		p = 1.50140941f + p * w;
	}
	else
	{
		w = sqrt(w) - 3.000000f;
		p = -0.000200214257f;
		p = 0.000100950558f + p * w;
		p = 0.00134934322f + p * w;
		p = -0.00367342844f + p * w;
		p = 0.00573950773f + p * w;
		p = -0.0076224613f + p * w;
		p = 0.00943887047f + p * w;
		p = 1.00167406f + p * w;
		p = 2.83297682f + p * w;
	}
	return p * x;
}

float CDF(float x, float mu, float sigma)
{
	float U = 0.5f * (1 + Erf((x-mu)/(sigma*sqrtf(2.0f))));
	return U;
}

float invCDF(float U, float mu, float sigma)
{
	float x = sigma*sqrtf(2.0f) * ErfInv(2.0f*U-1.0f) + mu;
	return x;
}

/*****************************************************************************/
/*****************************************************************************/
/**** Section 1.3.2 Applying the histogram transformation T on the input *****/
/*****************************************************************************/
/*****************************************************************************/

void ComputeTinput(TextureDataFloat& input, TextureDataFloat& T_input, int channel)
{	
	// Sort pixels of example image
	vector<PixelSortStruct> sortedInputValues;
	sortedInputValues.resize(input.width * input.height);
	for (int y = 0; y < input.height; y++)
	for (int x = 0; x < input.width; x++)
	{
		sortedInputValues[y * input.width + x].x = x;
		sortedInputValues[y * input.width + x].y = y;
		sortedInputValues[y * input.width + x].value = input.GetPixel(x, y, channel);
	}
	sort(sortedInputValues.begin(), sortedInputValues.end());

	// Assign Gaussian value to each pixel
	for (unsigned int i = 0; i < sortedInputValues.size() ; i++)
	{
		// Pixel coordinates
		int x = sortedInputValues[i].x;
		int y = sortedInputValues[i].y;
		// Input quantile (given by its order in the sorting)
		float U = (i + 0.5f) / (sortedInputValues.size());
		// Gaussian quantile
		float G = invCDF(U, GAUSSIAN_AVERAGE, GAUSSIAN_STD);
		// Store
		T_input.SetPixel(x, y, channel, G);
	}
}

/*****************************************************************************/
/*****************************************************************************/
/*  Section 1.3.3 Precomputing the inverse histogram transformation T^{-1}   */
/*****************************************************************************/
/*****************************************************************************/

void ComputeinvT(TextureDataFloat& input, TextureDataFloat& Tinv, int channel)
{
	// Sort pixels of example image
	vector<float> sortedInputValues;
	sortedInputValues.resize(input.width * input.height);
	for (int y = 0; y < input.height; y++)
	for (int x = 0; x < input.width; x++)
	{
		sortedInputValues[y * input.width + x] = input.GetPixel(x, y, channel);
	}
	sort(sortedInputValues.begin(), sortedInputValues.end());

	// Generate Tinv look-up table 
	for (int i = 0; i < Tinv.width; i++)
	{
		// Gaussian value in [0, 1]
		float G = (i + 0.5f) / (Tinv.width);
		// Quantile value 
		float U = CDF(G, GAUSSIAN_AVERAGE, GAUSSIAN_STD);
		// Find quantile in sorted pixel values
		int index = (int)floor(U * sortedInputValues.size());
		// Get input value 
		float I = sortedInputValues[index];
		// Store in LUT
		Tinv.SetPixel(i, 0, channel, I);
	}
}

/*****************************************************************************/
/*****************************************************************************/
/******** Section 1.4 Improvement: using a decorrelated color space **********/
/*****************************************************************************/
/*****************************************************************************/

// Compute the eigen vectors of the histogram of the input
void ComputeEigenVectors(TextureDataFloat& input, vec3 eigenVectors[3])
{
	// First and second order moments
	float R=0, G=0, B=0, RR=0, GG=0, BB=0, RG=0, RB=0, GB=0;
	for (int y = 0; y < input.height; y++)
	for (int x = 0; x < input.width; x++)
	{
		vec3 col = input.GetColorAt(x, y);
		R += col.x;
		G += col.y;
		B += col.z;
		RR += col.x * col.x;
		GG += col.y * col.y;
		BB += col.z * col.z;
		RG += col.x * col.y;
		RB += col.x * col.z;
		GB += col.y * col.z;
	}	
	R /= (float)(input.width * input.height);
	G /= (float)(input.width * input.height);
	B /= (float)(input.width * input.height);
	RR /= (float)(input.width * input.height);
	GG /= (float)(input.width * input.height);
	BB /= (float)(input.width * input.height);
	RG /= (float)(input.width * input.height);
	RB /= (float)(input.width * input.height);
	GB /= (float)(input.width * input.height);
	
	// Covariance matrix
	double covarMat[3][3];
	covarMat[0][0] = RR - R*R;
	covarMat[0][1] = RG - R*G;
	covarMat[0][2] = RB - R*B;
	covarMat[1][0] = RG - R*G;
	covarMat[1][1] = GG - G*G;
	covarMat[1][2] = GB - G*B;
	covarMat[2][0] = RB - R*B;
	covarMat[2][1] = GB - G*B;
	covarMat[2][2] = BB - B*B;

	// Find eigen values and vectors using Jacobi algorithm
	double eigenVectorsTemp[3][3];
	double eigenValuesTemp[3];
	ComputeEigenValuesAndVectors(covarMat, eigenVectorsTemp, eigenValuesTemp);

	// Set return values
	eigenVectors[0] = vec3((float)eigenVectorsTemp[0][0], (float)eigenVectorsTemp[1][0], (float)eigenVectorsTemp[2][0]);
	eigenVectors[1] = vec3((float)eigenVectorsTemp[0][1], (float)eigenVectorsTemp[1][1], (float)eigenVectorsTemp[2][1]);
	eigenVectors[2] = vec3((float)eigenVectorsTemp[0][2], (float)eigenVectorsTemp[1][2], (float)eigenVectorsTemp[2][2]);
}

// Main function of Section 1.4
void DecorrelateColorSpace(
 TextureDataFloat& input,			  // input: example image
 TextureDataFloat& input_decorrelated,// output: decorrelated input 
 vec3& colorSpaceVector1,			  // output: color space vector1 
 vec3& colorSpaceVector2,			  // output: color space vector2
 vec3& colorSpaceVector3,			  // output: color space vector3
 vec3& colorSpaceOrigin)			  // output: color space origin
{
	// Compute the eigenvectors of the histogram
	vec3 eigenvectors[3];
	ComputeEigenVectors(input, eigenvectors);

	// Rotate to eigenvector space and 
	for (int y = 0; y < input.height; y++)
	for (int x = 0; x < input.width; x++)
	for(int channel = 0 ; channel < 3 ; ++channel)
	{
		// Get current color
		vec3 color = input.GetColorAt(x, y);
		// Project on eigenvector 
		float new_channel_value = color.dot(eigenvectors[channel]);
		// Store
		input_decorrelated.SetPixel(x, y, channel, new_channel_value);
	}

	// Compute ranges of the new color space
	vec2 colorSpaceRanges[3] = {vec2(FLT_MAX,FLT_MIN), vec2(FLT_MAX,FLT_MIN), vec2(FLT_MAX,FLT_MIN)};
	for (int y = 0; y < input.height; y++)
	for (int x = 0; x < input.width; x++)
	for(int channel = 0 ; channel < 3 ; ++channel)
	{
		colorSpaceRanges[channel].x = min(colorSpaceRanges[channel].x, input_decorrelated.GetPixel(x,y,channel));
		colorSpaceRanges[channel].y = max(colorSpaceRanges[channel].y, input_decorrelated.GetPixel(x,y,channel));
	}

	// Remap range to [0, 1]
	for (int y = 0; y < input.height; y++)
	for (int x = 0; x < input.width; x++)
	for(int channel = 0 ; channel < 3 ; ++channel)
	{
		// Get current value
		float value = input_decorrelated.GetPixel(x,y,channel);
		// Remap in [0, 1]
		float remapped_value = (value - colorSpaceRanges[channel].x) / (colorSpaceRanges[channel].y - colorSpaceRanges[channel].x);
		// Store
		input_decorrelated.SetPixel(x, y, channel, remapped_value);
	}

	// Compute color space origin and vectors scaled for the normalized range
	colorSpaceOrigin.x = colorSpaceRanges[0].x * eigenvectors[0].x + colorSpaceRanges[1].x * eigenvectors[1].x + colorSpaceRanges[2].x * eigenvectors[2].x;
	colorSpaceOrigin.y = colorSpaceRanges[0].x * eigenvectors[0].y + colorSpaceRanges[1].x * eigenvectors[1].y + colorSpaceRanges[2].x * eigenvectors[2].y;
	colorSpaceOrigin.z = colorSpaceRanges[0].x * eigenvectors[0].z + colorSpaceRanges[1].x * eigenvectors[1].z + colorSpaceRanges[2].x * eigenvectors[2].z;
	colorSpaceVector1.x = eigenvectors[0].x * (colorSpaceRanges[0].y - colorSpaceRanges[0].x);
	colorSpaceVector1.y = eigenvectors[0].y * (colorSpaceRanges[0].y - colorSpaceRanges[0].x);
	colorSpaceVector1.z = eigenvectors[0].z * (colorSpaceRanges[0].y - colorSpaceRanges[0].x);
	colorSpaceVector2.x = eigenvectors[1].x * (colorSpaceRanges[1].y - colorSpaceRanges[1].x);
	colorSpaceVector2.y = eigenvectors[1].y * (colorSpaceRanges[1].y - colorSpaceRanges[1].x);
	colorSpaceVector2.z = eigenvectors[1].z * (colorSpaceRanges[1].y - colorSpaceRanges[1].x);
	colorSpaceVector3.x = eigenvectors[2].x * (colorSpaceRanges[2].y - colorSpaceRanges[2].x);
	colorSpaceVector3.y = eigenvectors[2].y * (colorSpaceRanges[2].y - colorSpaceRanges[2].x);
	colorSpaceVector3.z = eigenvectors[2].z * (colorSpaceRanges[2].y - colorSpaceRanges[2].x);
}

/*****************************************************************************/
/*****************************************************************************/
/* ===== Section 1.5 Improvement: prefiltering the look-up table =========== */
/*****************************************************************************/
/*****************************************************************************/

// Compute average subpixel variance at a given LOD
float ComputeLODAverageSubpixelVariance(TextureDataFloat& image, int LOD, int channel)
{
	// Window width associated with
	int windowWidth = 1 << LOD;
	
	// Compute average variance in all the windows
	float average_window_variance = 0.0;

	// Loop over al the windows
	for(int window_y = 0 ; window_y < image.height ; window_y += windowWidth)
	for(int window_x = 0 ; window_x < image.width  ; window_x += windowWidth)
	{
		// Estimate variance of current window
		float v = 0.0f; 
		float v2 = 0.0f; 
		for(int y = 0 ; y < windowWidth ; y++)
		for(int x = 0 ; x < windowWidth ; x++)
		{
			float value = image.GetPixel(window_x + x, window_y + y, channel);
			v  += value;
			v2 += value*value;
		}
		v  /= (float)(windowWidth*windowWidth);
		v2 /= (float)(windowWidth*windowWidth);
		float window_variance = max(0.0f, v2 - v*v);
		
		// Update average
		average_window_variance += window_variance / (image.width*image.height/windowWidth/windowWidth);
	}
	
	return average_window_variance;
}

// Filter LUT by sampling a Gaussian N(mu, std²)
float FilterLUTValueAtx(TextureDataFloat& LUT, float x, float std, int channel)
{
	// Number of samples for filtering (heuristic: twice the LUT resolution)
	const int numberOfSamples = 2 * LUT_WIDTH;

	// Filter
	float filtered_value = 0.0f;
	for (int sample = 0; sample < numberOfSamples; sample++)
	{
		// Quantile used to sample the Gaussian
		float U = (sample + 0.5f) / numberOfSamples;
		// Sample the Gaussian 
		float sample_x = invCDF(U, x, std);
		// Find sample texel in LUT (the LUT covers the domain [0, 1])
		int sample_texel = max(0, min(LUT_WIDTH - 1, (int)floor(sample_x * LUT_WIDTH)));
		// Fetch LUT at level 0
		float sample_value = LUT.GetPixel(sample_texel, 0, channel);
		// Accumulate
		filtered_value += sample_value;
	}
	// Normalize and return
	filtered_value /= (float) numberOfSamples;
	return filtered_value;
}

// Main function of section 1.5
void PrefilterLUT(TextureDataFloat& image_T_Input, TextureDataFloat& LUT_Tinv, int channel)
{
	// Prefilter 
	for(int LOD = 1 ; LOD < LUT_Tinv.height ; LOD++)
	{
		// Compute subpixel variance at LOD 
		float window_variance = ComputeLODAverageSubpixelVariance(image_T_Input, LOD, channel);
		float window_std = sqrtf(window_variance);

		// Prefilter LUT with Gaussian kernel of this variance
		for (int i = 0; i < LUT_Tinv.width; i++)
		{
			// Texel position in [0, 1]
			float x_texel = (i+0.5f) / LUT_Tinv.width;
			// Filter look-up table around this position with Gaussian kernel
			float filteredValue = FilterLUTValueAtx(LUT_Tinv, x_texel, window_std, channel);
			// Store filtered value
			LUT_Tinv.SetPixel(i, LOD, channel, filteredValue);
		}
	}
}

/*********************************************************************/
/*********************************************************************/
/*************************** Main Function ***************************/
/*********************************************************************/
/*********************************************************************/

void Precomputations(
	TextureDataFloat& input,   // input: example image
	TextureDataFloat& Tinput,  // output: T(input) image
	TextureDataFloat& Tinv,	   // output: T^{-1} look-up table
	vec3& colorSpaceVector1,   // output: color space vector1 
	vec3& colorSpaceVector2,   // output: color space vector2
	vec3& colorSpaceVector3,   // output: color space vector3
	vec3& colorSpaceOrigin)    // output: color space origin
{	
	// Section 1.4 Improvement: using a decorrelated color space
	TextureDataFloat input_decorrelated = TextureDataFloat(input.width, input.height, 3);
	DecorrelateColorSpace(input, input_decorrelated, colorSpaceVector1, colorSpaceVector2, colorSpaceVector3, colorSpaceOrigin);

	// Section 1.3.2 Applying the histogram transformation T on the input
	Tinput = TextureDataFloat(input.width, input.height, 3);
#pragma omp parallel for schedule(static,1)
	for(int channel = 0 ; channel < 3 ; channel++) {
		ComputeTinput(input_decorrelated, Tinput, channel);
	}

	// Section 1.3.3 Precomputing the inverse histogram transformation T^{-1}
	Tinv = TextureDataFloat(LUT_WIDTH, 1, 3);
#pragma omp parallel for schedule(static,1)
	for(int channel = 0 ; channel < 3 ; channel++) {
		ComputeinvT(input_decorrelated, Tinv, channel);
	}

	// Section 1.5 Improvement: prefiltering the look-up table
	// Compute number of prefiltered levels and resize LUT
	Tinv.height = (int)(log((float)Tinput.width)/log(2.0f));
	Tinv.data.resize(3 * Tinv.width * Tinv.height);

#pragma omp parallel for schedule(static,1)
	for(int channel = 0 ; channel < 3 ; channel++) {		
		PrefilterLUT(Tinput, Tinv, channel);
	}
}
