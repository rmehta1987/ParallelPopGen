/*
 * shared.cuh
 *
 *      Author: David Lawrie
 *      for cuda and rand functions used by both go_fish and by sfs
 */

#ifndef SHARED_CUH_
#define SHARED_CUH_

//includes below in sfs & go_fish
#include <cuda_runtime.h>
#include "../_outside_libraries/helper_math.h"
#include <limits.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "../_outside_libraries/Random123/philox.h"
#include "../_outside_libraries/Random123/features/compilerfeatures.h"

/* ----- cuda error checking & device setting ----- */
#define __DEBUG__ false
#define cudaCheckErrors(expr1,expr2,expr3) { cudaError_t e = expr1; int g = expr2; int p = expr3; if (e != cudaSuccess) { fprintf(stderr,"error %d %s\tfile %s\tline %d\tgeneration %d\t population %d\n", e, cudaGetErrorString(e),__FILE__,__LINE__, g,p); exit(1); } }
#define cudaCheckErrorsAsync(expr1,expr2,expr3) { cudaCheckErrors(expr1,expr2,expr3); if(__DEBUG__){ cudaCheckErrors(cudaDeviceSynchronize(),expr2,expr3); } }

__forceinline__ cudaDeviceProp set_cuda_device(int & cuda_device){
	int cudaDeviceCount;
	cudaCheckErrorsAsync(cudaGetDeviceCount(&cudaDeviceCount),-1,-1);
	if(cuda_device >= 0 && cuda_device < cudaDeviceCount){ cudaCheckErrors(cudaSetDevice(cuda_device),-1,-1); } //unless user specifies, driver auto-magically selects free GPU to run on
	int myDevice;
	cudaCheckErrorsAsync(cudaGetDevice(&myDevice),-1,-1);
	cudaDeviceProp devProp;
	cudaCheckErrors(cudaGetDeviceProperties(&devProp, myDevice),-1,-1);
	cuda_device = myDevice;
	return devProp;
}

/* ----- end cuda error checking ----- */

/* ----- random number generation ----- */

namespace RNG{
#define RNG_MEAN_BOUNDARY_NORM 6
#define RNG_N_BOUNDARY_POIS_BINOM 100 //binomial calculation starts to become numerically unstable for large values of N, not sure where that starts but is between 200 and 200,000

// uint_float_01: Input is a W-bit integer (unsigned).  It is multiplied
// by Float(2^-W) and added to Float(2^(-W-1)).  A good compiler should
// optimize it down to an int-to-float conversion followed by a multiply
// and an add, which might be fused, depending on the architecture.
//
// If the input is a uniformly distributed integer, then the
// result is a uniformly distributed floating point number in [0, 1].
// The result is never exactly 0.0.
// The smallest value returned is 2^-W.
// Let M be the number of mantissa bits in Float.
// If W>M  then the largest value retured is 1.0.
// If W<=M then the largest value returned is the largest Float less than 1.0.
__host__ __device__ __forceinline__ float uint_float_01(unsigned int in){
	//(mostly) stolen from Philox code "uniform.hpp"
	R123_CONSTEXPR float factor = float(1.)/(UINT_MAX + float(1.));
	R123_CONSTEXPR float halffactor = float(0.5)*factor;
    return in*factor + halffactor;
}


__host__ __device__ __forceinline__  uint4 Philox(int2 seed, int k, int step, int population, int round){
	typedef r123::Philox4x32_R<10> P; //can change the 10 rounds of bijection down to 8 (lowest safe limit) to get possible extra speed!
	P rng;

	P::key_type key = {{seed.x, seed.y}}; //random int to set key space + seed
	P::ctr_type count = {{k, step, population, round}};

	union {
		P::ctr_type c;
		uint4 i;
	}u;

	u.c = rng(count, key);

	return u.i;
}

__host__ __device__ __forceinline__ void binom_iter(float j, float x, float n, float & emu, float & cdf){
	emu *= ((n+1.f-j)*x)/(j*(1-x));
	cdf += emu;
}

__host__ __device__ __forceinline__ int binomcdfinv(float r, float mean, float x, float n){
	float emu = powf(1-x,n);
	if(emu == 1) { emu = expf(-1 * mean);  }
	float cdf = emu;
	if(cdf >= r){ return 0; }

	binom_iter(1.f, x, n, emu, cdf); if(cdf >= r){ return 1; }
	binom_iter(2.f, x, n, emu, cdf); if(cdf >= r){ return 2; }
	binom_iter(3.f, x, n, emu, cdf); if(cdf >= r){ return 3; }
	binom_iter(4.f, x, n, emu, cdf); if(cdf >= r){ return 4; }
	binom_iter(5.f, x, n, emu, cdf); if(cdf >= r){ return 5; }
	binom_iter(6.f, x, n, emu, cdf); if(cdf >= r){ return 6; }
	binom_iter(7.f, x, n, emu, cdf); if(cdf >= r){ return 7; }
	binom_iter(8.f, x, n, emu, cdf); if(cdf >= r){ return 8; }
	binom_iter(9.f, x, n, emu, cdf); if(cdf >= r){ return 9; }
	binom_iter(10.f, x, n, emu, cdf); if(cdf >= r){ return 10; }
	binom_iter(11.f, x, n, emu, cdf); if(cdf >= r || mean <= 1){ return 11; }
	binom_iter(12.f, x, n, emu, cdf); if(cdf >= r){ return 12; }
	binom_iter(13.f, x, n, emu, cdf); if(cdf >= r){ return 13; }
	binom_iter(14.f, x, n, emu, cdf); if(cdf >= r || mean <= 2){ return 14; }
	binom_iter(15.f, x, n, emu, cdf); if(cdf >= r){ return 15; }
	binom_iter(16.f, x, n, emu, cdf); if(cdf >= r){ return 16; }
	binom_iter(17.f, x, n, emu, cdf); if(cdf >= r || mean <= 3){ return 17; }
	binom_iter(18.f, x, n, emu, cdf); if(cdf >= r){ return 18; }
	binom_iter(19.f, x, n, emu, cdf); if(cdf >= r){ return 19; }
	binom_iter(20.f, x, n, emu, cdf); if(cdf >= r || mean <= 4){ return 20; }
	binom_iter(21.f, x, n, emu, cdf); if(cdf >= r){ return 21; }
	binom_iter(22.f, x, n, emu, cdf); if(cdf >= r || mean <= 5){ return 22; }
	binom_iter(23.f, x, n, emu, cdf); if(cdf >= r){ return 23; }
	binom_iter(24.f, x, n, emu, cdf); if(cdf >= r || mean <= 6){ return 24; }
	binom_iter(25.f, x, n, emu, cdf); if(cdf >= r){ return 25; }
	binom_iter(26.f, x, n, emu, cdf); if(cdf >= r || mean <= 7){ return 26; }
	binom_iter(27.f, x, n, emu, cdf); if(cdf >= r){ return 27; }
	binom_iter(28.f, x, n, emu, cdf); if(cdf >= r || mean <= 8){ return 28; }
	binom_iter(29.f, x, n, emu, cdf); if(cdf >= r){ return 29; }
	binom_iter(30.f, x, n, emu, cdf); if(cdf >= r || mean <= 9){ return 30; }
	binom_iter(31.f, x, n, emu, cdf); if(cdf >= r){ return 31; }
	binom_iter(32.f, x, n, emu, cdf); if(cdf >= r || mean <= 10){ return 32; }
	binom_iter(33.f, x, n, emu, cdf); if(cdf >= r){ return 33; }
	binom_iter(34.f, x, n, emu, cdf); if(cdf >= r || mean <= 11){ return 34; }
	binom_iter(35.f, x, n, emu, cdf); if(cdf >= r){ return 35; }
	binom_iter(36.f, x, n, emu, cdf); if(cdf >= r || mean <= 12){ return 36; }
	binom_iter(37.f, x, n, emu, cdf); if(cdf >= r){ return 37; }
	binom_iter(38.f, x, n, emu, cdf); if(cdf >= r || mean <= 13){ return 38; }
	binom_iter(39.f, x, n, emu, cdf); if(cdf >= r){ return 39; }
	binom_iter(40.f, x, n, emu, cdf); if(cdf >= r || mean <= 14){ return 40; }
	binom_iter(41.f, x, n, emu, cdf); if(cdf >= r || mean <= 15){ return 41; }
	binom_iter(42.f, x, n, emu, cdf); if(cdf >= r){ return 42; }
	binom_iter(43.f, x, n, emu, cdf); if(cdf >= r || mean <= 16){ return 43; }
	binom_iter(44.f, x, n, emu, cdf); if(cdf >= r){ return 44; }
	binom_iter(45.f, x, n, emu, cdf); if(cdf >= r || mean <= 17){ return 45; }
	binom_iter(46.f, x, n, emu, cdf); if(cdf >= r || mean <= 18){ return 46; }
	binom_iter(47.f, x, n, emu, cdf); if(cdf >= r){ return 47; }
	binom_iter(48.f, x, n, emu, cdf); if(cdf >= r || mean <= 19){ return 48; }
	binom_iter(49.f, x, n, emu, cdf); if(cdf >= r){ return 49; }
	binom_iter(50.f, x, n, emu, cdf); if(cdf >= r || mean <= 20){ return 50; }
	binom_iter(51.f, x, n, emu, cdf); if(cdf >= r || mean <= 21){ return 51; }
	binom_iter(52.f, x, n, emu, cdf); if(cdf >= r){ return 52; }
	binom_iter(53.f, x, n, emu, cdf); if(cdf >= r || mean <= 22){ return 53; }
	binom_iter(54.f, x, n, emu, cdf); if(cdf >= r){ return 54; }
	binom_iter(55.f, x, n, emu, cdf); if(cdf >= r || mean <= 23){ return 55; }
	binom_iter(56.f, x, n, emu, cdf); if(cdf >= r || mean <= 24){ return 56; }
	binom_iter(57.f, x, n, emu, cdf); if(cdf >= r){ return 57; }
	binom_iter(58.f, x, n, emu, cdf); if(cdf >= r || mean <= 25){ return 58; }
	binom_iter(59.f, x, n, emu, cdf); if(cdf >= r || mean <= 26){ return 59; }
	binom_iter(60.f, x, n, emu, cdf); if(cdf >= r){ return 60; }
	binom_iter(61.f, x, n, emu, cdf); if(cdf >= r || mean <= 27){ return 61; }
	binom_iter(62.f, x, n, emu, cdf); if(cdf >= r || mean <= 28){ return 62; }
	binom_iter(63.f, x, n, emu, cdf); if(cdf >= r){ return 63; }
	binom_iter(64.f, x, n, emu, cdf); if(cdf >= r || mean <= 29){ return 64; }
	binom_iter(65.f, x, n, emu, cdf); if(cdf >= r || mean <= 30){ return 65; }
	binom_iter(66.f, x, n, emu, cdf); if(cdf >= r){ return 66; }
	binom_iter(67.f, x, n, emu, cdf); if(cdf >= r || mean <= 31){ return 67; }
	binom_iter(68.f, x, n, emu, cdf); if(cdf >= r || mean <= 32){ return 68; }
	binom_iter(69.f, x, n, emu, cdf); if(cdf >= r){ return 69; }

	return 70; //17 for mean <= 3, 24 limit for mean <= 6, 32 limit for mean <= 10, 36 limit for mean <= 12, 41 limit for mean <= 15, 58 limit for mean <= 25, 70 limit for mean <= 33; max float between 0 and 1 is 0.99999999
}

__host__ __device__ __forceinline__ void pois_iter(float j, float mean, float & emu, float & cdf){
	emu *= mean/j;
	cdf += emu;
}

__host__ __device__ __forceinline__ int poiscdfinv(float r, float mean){
	float emu = expf(-1 * mean);
	float cdf = emu;
	if(cdf >= r){ return 0; }

	pois_iter(1.f, mean, emu, cdf); if(cdf >= r){ return 1; }
	pois_iter(2.f, mean, emu, cdf); if(cdf >= r){ return 2; }
	pois_iter(3.f, mean, emu, cdf); if(cdf >= r){ return 3; }
	pois_iter(4.f, mean, emu, cdf); if(cdf >= r){ return 4; }
	pois_iter(5.f, mean, emu, cdf); if(cdf >= r){ return 5; }
	pois_iter(6.f, mean, emu, cdf); if(cdf >= r){ return 6; }
	pois_iter(7.f, mean, emu, cdf); if(cdf >= r){ return 7; }
	pois_iter(8.f, mean, emu, cdf); if(cdf >= r){ return 8; }
	pois_iter(9.f, mean, emu, cdf); if(cdf >= r){ return 9; }
	pois_iter(10.f, mean, emu, cdf); if(cdf >= r){ return 10; }
	pois_iter(11.f, mean, emu, cdf); if(cdf >= r || mean <= 1){ return 11; }
	pois_iter(12.f, mean, emu, cdf); if(cdf >= r){ return 12; }
	pois_iter(13.f, mean, emu, cdf); if(cdf >= r){ return 13; }
	pois_iter(14.f, mean, emu, cdf); if(cdf >= r || mean <= 2){ return 14; }
	pois_iter(15.f, mean, emu, cdf); if(cdf >= r){ return 15; }
	pois_iter(16.f, mean, emu, cdf); if(cdf >= r){ return 16; }
	pois_iter(17.f, mean, emu, cdf); if(cdf >= r || mean <= 3){ return 17; }
	pois_iter(18.f, mean, emu, cdf); if(cdf >= r){ return 18; }
	pois_iter(19.f, mean, emu, cdf); if(cdf >= r){ return 19; }
	pois_iter(20.f, mean, emu, cdf); if(cdf >= r || mean <= 4){ return 20; }
	pois_iter(21.f, mean, emu, cdf); if(cdf >= r){ return 21; }
	pois_iter(22.f, mean, emu, cdf); if(cdf >= r || mean <= 5){ return 22; }
	pois_iter(23.f, mean, emu, cdf); if(cdf >= r){ return 23; }
	pois_iter(24.f, mean, emu, cdf); if(cdf >= r || mean <= 6){ return 24; }
	pois_iter(25.f, mean, emu, cdf); if(cdf >= r){ return 25; }
	pois_iter(26.f, mean, emu, cdf); if(cdf >= r || mean <= 7){ return 26; }
	pois_iter(27.f, mean, emu, cdf); if(cdf >= r){ return 27; }
	pois_iter(28.f, mean, emu, cdf); if(cdf >= r || mean <= 8){ return 28; }
	pois_iter(29.f, mean, emu, cdf); if(cdf >= r){ return 29; }
	pois_iter(30.f, mean, emu, cdf); if(cdf >= r || mean <= 9){ return 30; }
	pois_iter(31.f, mean, emu, cdf); if(cdf >= r){ return 31; }
	pois_iter(32.f, mean, emu, cdf); if(cdf >= r || mean <= 10){ return 32; }
	pois_iter(33.f, mean, emu, cdf); if(cdf >= r){ return 33; }
	pois_iter(34.f, mean, emu, cdf); if(cdf >= r || mean <= 11){ return 34; }
	pois_iter(35.f, mean, emu, cdf); if(cdf >= r){ return 35; }
	pois_iter(36.f, mean, emu, cdf); if(cdf >= r || mean <= 12){ return 36; }
	pois_iter(37.f, mean, emu, cdf); if(cdf >= r){ return 37; }
	pois_iter(38.f, mean, emu, cdf); if(cdf >= r || mean <= 13){ return 38; }
	pois_iter(39.f, mean, emu, cdf); if(cdf >= r){ return 39; }
	pois_iter(40.f, mean, emu, cdf); if(cdf >= r || mean <= 14){ return 40; }
	pois_iter(41.f, mean, emu, cdf); if(cdf >= r || mean <= 15){ return 41; }
	pois_iter(42.f, mean, emu, cdf); if(cdf >= r){ return 42; }
	pois_iter(43.f, mean, emu, cdf); if(cdf >= r || mean <= 16){ return 43; }
	pois_iter(44.f, mean, emu, cdf); if(cdf >= r){ return 44; }
	pois_iter(45.f, mean, emu, cdf); if(cdf >= r || mean <= 17){ return 45; }
	pois_iter(46.f, mean, emu, cdf); if(cdf >= r || mean <= 18){ return 46; }
	pois_iter(47.f, mean, emu, cdf); if(cdf >= r){ return 47; }
	pois_iter(48.f, mean, emu, cdf); if(cdf >= r || mean <= 19){ return 48; }
	pois_iter(49.f, mean, emu, cdf); if(cdf >= r){ return 49; }
	pois_iter(50.f, mean, emu, cdf); if(cdf >= r || mean <= 20){ return 50; }
	pois_iter(51.f, mean, emu, cdf); if(cdf >= r || mean <= 21){ return 51; }
	pois_iter(52.f, mean, emu, cdf); if(cdf >= r){ return 52; }
	pois_iter(53.f, mean, emu, cdf); if(cdf >= r || mean <= 22){ return 53; }
	pois_iter(54.f, mean, emu, cdf); if(cdf >= r){ return 54; }
	pois_iter(55.f, mean, emu, cdf); if(cdf >= r || mean <= 23){ return 55; }
	pois_iter(56.f, mean, emu, cdf); if(cdf >= r || mean <= 24){ return 56; }
	pois_iter(57.f, mean, emu, cdf); if(cdf >= r){ return 57; }
	pois_iter(58.f, mean, emu, cdf); if(cdf >= r || mean <= 25){ return 58; }
	pois_iter(59.f, mean, emu, cdf); if(cdf >= r || mean <= 26){ return 59; }
	pois_iter(60.f, mean, emu, cdf); if(cdf >= r){ return 60; }
	pois_iter(61.f, mean, emu, cdf); if(cdf >= r || mean <= 27){ return 61; }
	pois_iter(62.f, mean, emu, cdf); if(cdf >= r || mean <= 28){ return 62; }
	pois_iter(63.f, mean, emu, cdf); if(cdf >= r){ return 63; }
	pois_iter(64.f, mean, emu, cdf); if(cdf >= r || mean <= 29){ return 64; }
	pois_iter(65.f, mean, emu, cdf); if(cdf >= r || mean <= 30){ return 65; }
	pois_iter(66.f, mean, emu, cdf); if(cdf >= r){ return 66; }
	pois_iter(67.f, mean, emu, cdf); if(cdf >= r || mean <= 31){ return 67; }
	pois_iter(68.f, mean, emu, cdf); if(cdf >= r || mean <= 32){ return 68; }
	pois_iter(69.f, mean, emu, cdf); if(cdf >= r){ return 69; }

	return 70; //17 for mean <= 3, 24 limit for mean <= 6, 32 limit for mean <= 10, 36 limit for mean <= 12, 41 limit for mean <= 15, 58 limit for mean <= 25, 70 limit for mean <= 33; max float between 0 and 1 is 0.99999999
}
//define RNG_MEAN_BOUNDARY_NORM 6 -- just to see values
//RNG_N_BOUNDARY_POIS_BINOM 100 - -- just to see values
__host__ __device__ __forceinline__ int ApproxRandPois1(float mean, float var, float p, float N, int2 seed, int id, int generation, int population){
	uint4 i = Philox(seed, id, generation, population, 0);
	//printf("The mean in approxrandpois %f\n", mean);
	if(mean <= RNG_MEAN_BOUNDARY_NORM)
	{ 
		return poiscdfinv(uint_float_01(i.x), mean); 
	}
	else if(mean >= N-RNG_MEAN_BOUNDARY_NORM)  //flip side of poisson, when 1-p is small
	{ 
		return N - poiscdfinv(uint_float_01(i.x), N-mean); 
	} 
	return round(normcdfinv(uint_float_01(i.x))*sqrtf(var)+mean);
}

__host__ __device__ __forceinline__ int ApproxRandBinom1(float mean, float var, float p, float N, int2 seed, int id, int generation, int population){
	uint4 i = Philox(seed, id, generation, population, 0);

	//printf("The mean in ApproxRandBinom1 %f\n", mean);
	
	if(mean <= RNG_MEAN_BOUNDARY_NORM)
	{
		if(N < RNG_N_BOUNDARY_POIS_BINOM)
		{ 
			return binomcdfinv(uint_float_01(i.x), mean, mean/N, N); 
		} 
		else
		{ 
			return poiscdfinv(uint_float_01(i.x), mean); 
		}
	}
	else if(mean >= N-RNG_MEAN_BOUNDARY_NORM)
	{ //flip side of binomial, when 1-p is small
		if(N < RNG_N_BOUNDARY_POIS_BINOM)
		{ 
			return N - binomcdfinv(uint_float_01(i.x), N-mean, (N-mean)/N, N); 
		} 
		else
		{ 
			return N - poiscdfinv(uint_float_01(i.x), N-mean); 
		}
	}
	return round(normcdfinv(uint_float_01(i.x))*sqrtf(var)+mean);
}

//faster on if don't inline on both GPUs!
__device__ int ApproxRandBinomHelper(unsigned int i, float mean, float var, float N);

__device__ __forceinline__ int4 ApproxRandBinom4(float4 mean, float4 var, float4 p, float N, int2 seed, int id, int generation, int population)
{
	uint4 i = Philox(seed, id, generation, population, 0);
	return make_int4(ApproxRandBinomHelper(i.x, mean.x, var.x, N), ApproxRandBinomHelper(i.y, mean.y, var.y, N), ApproxRandBinomHelper(i.z, mean.z, var.z, N), ApproxRandBinomHelper(i.w, mean.w, var.w, N));
}
/* ----- end random number generation ----- */


/**** THIS USES THE POISSONINV CDF FROM MIKE GILES SEE CODE AND REFERNCES BELOW*/

/*

    This software was written by Mike Giles, copyright University of Oxford,
    and is provided under the terms of the GNU GPLv3 license:
    http://www.gnu.org/licenses/gpl.html

    Commercial users who would like to use the software under a more
    permissive license, such as BSD, should contact the author:
    mike.giles@maths.ox.ac.uk

*/

// include CUDA header which defines Inf and NaN constants



//
// This single precision function computes the inverse
// of the Poisson CDF, using at most 10 registers
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// max |error| no more than 1 for lam < ~1.6e+07
//
// Compilation with use_fast_math. Extra error for 5.5 < lam <= 10 due to
// conversion of epxf(0.5f*lam) into single precision intrinsics
// __expf(0.5f*lam)in the direct summation.
// ave |error| < 6.2e-07*sqrt(max(1e-02,lam)   for 0   <= lam <= 5.5
//             < 2.0e-06                       for 5.5 <  lam <= 10
//             < 1.8e-07*sqrt(max(69   ,lam)   for 10  <  lam <  ~1.6e+07
//
// Compilation without use_fast_math
// ave |error| < 6.2e-07*sqrt(max(1e-02,lam)   for 0    <= lam <= 7.5
//             < 1.8e-07*sqrt(max(89   ,lam)   for 7.5  <  lam <  ~1.6e+07
//
// For lam > ~1.6e+07, the errors will be about 1 ulp.
//

__device__ static double rlog1(double);

__device__ inline float poissinvf(float u, float lam) {

    float s, t, x=0.0f;

    // handle exceptions

    if (lam < 0.0f || isinf(lam)) return CUDART_NAN_F;
    if (u  < 0.0f || u > 1.0f) return CUDART_NAN_F;
    if (u == 0.0f || lam == 0.0f) return 0.0f;
    if (u == 1.0f) return CUDART_INF_F;

    // large lam
    // NOTE: need threshold for large Lam > 4.0.
    // For fixed lam the threshold Lam > 10.0 is faster
    // but for mixed lam it is better to choose Lam > 4.0 due
    // to warp divergence

    if (lam > 4.0f) {
        s = normcdfinvf(u)*rsqrtf(lam);

        // use polynomial approximations in central region

        if ((s > -0.6833501f) && (s < 1.777993f)) {
            float rm;

            //  polynomial approximation to f^{-1}(s) - 1

            rm =  2.82298751e-07f;
            rm = -2.58136133e-06f + rm*s;
            rm =  1.02118025e-05f + rm*s;
            rm = -2.37996199e-05f + rm*s;
            rm =  4.05347462e-05f + rm*s;
            rm = -6.63730967e-05f + rm*s;
            rm =  0.000124762566f + rm*s;
            rm = -0.000256970731f + rm*s;
            rm =  0.000558953132f + rm*s;
            rm =  -0.00133129194f + rm*s;
            rm =   0.00370367937f + rm*s;
            rm =   -0.0138888706f + rm*s;
            rm =     0.166666667f + rm*s;
            rm =             s + s*(rm*s);

            //  polynomial approximation to correction c0(r)

            t  =   1.86386867e-05f;
            t  =  -0.000207319499f + t*rm;
            t  =     0.0009689451f + t*rm;
            t  =   -0.00247340054f + t*rm;
            t  =    0.00379952985f + t*rm;
            t  =   -0.00386717047f + t*rm;
            t  =    0.00346960934f + t*rm;
            t  =   -0.00414125511f + t*rm;
            t  =    0.00586752093f + t*rm;
            t  =   -0.00838583787f + t*rm;
            t  =     0.0132793933f + t*rm;
            t  =     -0.027775536f + t*rm;
            t  =      0.333333333f + t*rm;

            //  O(1/lam) correction

            x  =  -0.000134529865f;
            x  =     0.0013711034f + x*rm;
            x  =   -0.00583140335f + x*rm;
            x  =      0.013452415f + x*rm;
            x  =    -0.0185780057f + x*rm;
            x  =     0.0169826221f + x*rm;
            x  =    -0.0135307848f + x*rm;
            x  =     0.0135629517f + x*rm;
            x  =    -0.0155146497f + x*rm;
            x  =     0.0174051522f + x*rm;
            x  =    -0.0198016963f + x*rm;
            x  = __fdividef(x,lam);

            //    sum from smallest to largest to minimise rounding error;
            //    use of __fadd_rd to round down final sum is important for
            //    very large values of lambda to ensure correct rounding

            x = truncf( __fadd_rd(lam, (x+t)+lam*rm) );
        }

        // otherwise use Newton iteration

        else if (s > -sqrtf(2.0f)) {
            float r, r2, s2;

            r = 1.0f + s;
            if (r < 0.1f) r = 0.1f;

            do {
                t  = __logf(r);
                r2 = r;
                s2 = sqrtf(2.0f*((1.0f-r) + r*t));
                if (r < 1.0f) s2 = -s2;
                r = r2 - (s2-s)*s2/t;
                if (r < 0.1f*r2) r = 0.1f*r2;
            } while (fabsf(r-r2) > 1e-5f);

            t = __logf(r);
            x = lam*r + __logf(sqrtf(2.0f*r*((1.0f-r)+r*t))/fabsf(r-1.0f)) / t;
            x = truncf( x - 0.0218f/(x+0.065f*lam) );
        }
    }

    // bottom-up summation

    if (x < 10.0f) {
        float del;

        x    = 0.0f;
        t    = expf(0.5f*lam);
        del  = 0.0f;
        if (u > 0.5f)
            del  = t*(1e-6f*t);
        s    = 1.0f - t*(u*t) + del;

        while (s < 0.0f) {
            x  += 1.0f;
            t   = x/lam;
            del = t*del;
            s   = t*s + 1.0f;
        }

        // top-down summation if needed

        if (s < 2.0f*del) {
            del = 1e6f*del;
            t   = 1e7f*del;
            del = (1.0f-u)*del;

            while (del < t) {
                x   += 1.0f;
                del *= x/lam;
            }

            s = del;
            t = 1.0f;
            while (s > 0.0f) {
                t *= x/lam;
                s -= t;
                x -= 1.0f;
            }
        }
    }

    return x;
}


//
// This double precision function computes the inverse
// of the Poisson CDF, using about 30 registers
//
// u   = CDF value in range (0,1)
// lam = Poisson rate
//
// max |error| no more than 1 for lam < ~9e+15
// ave |error| < 1.2e-15*sqrt(max(1e-02,lam))    for 0 <= lam <= 6
//             < 2.0e-16*sqrt(max(216  ,lam))    for 6 <  lam < ~9e+15
//
// For lam >= ~9e+15, the errors will be about 1 ulp.
//

//
// naming convention: double precision variables are capitalised
//

__device__ inline double poissinv(double U, double Lam) {

    double X=0.0, Xi, S, T, Rm;

    // handle exceptions

    if ((Lam < 0.0) || isinf(Lam)) return CUDART_NAN;
    if (U <  0.0 || U > 1.0) return CUDART_NAN;
    if (U == 0.0 || Lam == 0.0) return 0.0;
    if (U == 1.0) return CUDART_INF;

    // large lam
    // NOTE: need threshold for large Lam > 4.0.
    // For fixed lam the threshold Lam > 15.0 is faster
    // but for mixed lam it is better to choose Lam > 4.0 due
    // to warp divergence

    if (Lam > 4.0) {
        float s, t, del, rm;

        S   = normcdfinv(U)*rsqrt(Lam);
        s   = (float) S;

        // del = 2.095e-06f is sufficient for U > 1e-195, but for very small U
        // a larger del = 3.3e-06f is correct. Minimal effect on samples/sec
        // from changing del to the smaller of the two.

        del = 3.3e-6f;

        // use polynomial approximations in central region

        if ((s > -0.6833501f) && (s < 1.777993f)) {
            float x;

            //  polynomial approximation to f^{-1}(s) - 1

            rm =    2.82298751e-07f;
            rm =   -2.58136133e-06f + rm*s;
            rm =    1.02118025e-05f + rm*s;
            rm =   -2.37996199e-05f + rm*s;
            rm =    4.05347462e-05f + rm*s;
            rm =   -6.63730967e-05f + rm*s;
            rm =    0.000124762566f + rm*s;
            rm =   -0.000256970731f + rm*s;
            rm =    0.000558953132f + rm*s;
            rm =    -0.00133129194f + rm*s;
            rm =     0.00370367937f + rm*s;
            rm =     -0.0138888706f + rm*s;
            Rm = 0.1666666666666667 + rm*s;

            S +=                  S*(Rm*S);
            Rm = S;
            rm = (float) Rm;

            //  polynomial approximation to correction c0(r)

            t  =   1.86386867e-05f;
            t  =  -0.000207319499f + t*rm;
            t  =     0.0009689451f + t*rm;
            t  =   -0.00247340054f + t*rm;
            t  =    0.00379952985f + t*rm;
            t  =   -0.00386717047f + t*rm;
            t  =    0.00346960934f + t*rm;
            t  =   -0.00414125511f + t*rm;
            t  =    0.00586752093f + t*rm;
            t  =   -0.00838583787f + t*rm;
            t  =     0.0132793933f + t*rm;
            t  =     -0.027775536f + t*rm;
            t  =      0.333333333f + t*rm;

            //  O(1/lam) correction

            x  =  -0.000134549989f;
            x  =    0.00137126732f + x*rm;
            x  =    -0.0058318548f + x*rm;
            x  =     0.0134526835f + x*rm;
            x  =    -0.0185770659f + x*rm;
            x  =     0.0169808982f + x*rm;
            x  =    -0.0135302102f + x*rm;
            x  =     0.0135636388f + x*rm;
            x  =      -0.01551508f + x*rm;
            x  =     0.0174051003f + x*rm;
            x  =    -0.0198016502f + x*rm;

            x  = __fdividef(x,(float) Lam);

            //    sum from smallest to largest to minimise rounding error;
            S = ((x + del) + t) + Lam*Rm;
        }

        // otherwise use Newton iteration

        else if (s > -sqrtf(2.0f)) {
            float r, r2, s2;

            r = 1.0f + s;
            if (r < 0.1f) r = 0.1f;

            do {
                t  = __logf(r);
                r2 = r;
                s2 = sqrtf(2.0f*((1.0f-r) + r*t));
                if (r < 1.0f) s2 = -s2;
                r = r2 - (s2-s)*s2/t;
                if (r < 0.1f*r2) r = 0.1f*r2;
            } while (fabsf(r-r2) > 1e-5f);

            t   = __logf(r);
            rm  = r - 1.0f;
            t   = __logf(sqrtf(2.0f*r*(-rm+r*t))/fabsf(rm)) / t;

            // O(1/lam) ad-hoc correction
            t  += -0.0218/(Lam*(0.065 + r));
            del =  0.01/(Lam*r);

            S = (del + t) + Lam*rm;
        }

        // if x>10, round down to nearest integer, and check accuracy

        // NOTE: use of __dadd_rd to round down final sum is important
        // for large values of lambda to ensure fast execution and accuracy
        // If Lam >> 1 then del can be smaller than ULP and the standard
        // IEEE_754 rounding mode (round to nearest, ties to even)
        // causes issues when the asymptotic approximation returns an
        // answer very close to an integer value. Notably it results in a
        // significant slowdown as the CDF must be checked increasingly
        // often. Round down mode for the final step of the summation in
        // the asymptotic approximation (somewhat) mitigates this.
        X = __dadd_rd(Lam, S);
        X = trunc(X);

        if (X > 10.0) {
            // If NO_CDF is defined we skip this step, this increases speed
            // but can significantly impact accuracy
#ifdef NO_CDF
            return X;
#endif

            // Accurate enough, no need to check the CDF
            if (S - 2.0f*del >= X - Lam) return X;

            // correction procedure based on Temme approximation (double precision)

            if (X > 0.5*Lam && X < 2.0*Lam) {
                double Eta, B0, B1;

                Xi = 1.0 / X;
                // Use rlog1(z) = z - log(1+z) function, a fast rational
                // minimax approximation for -0.6 < z < 1.6. This is a lot
                // more accurate when x is near 0 (here X near Lambda).
                Eta = sqrt(2.0*rlog1((Lam-X)*Xi));
                if (X > Lam) Eta = -Eta;

                B1 =  8.0995211567045583e-16;              S = B1;
                B0 = -1.9752288294349411e-15;              S = B0 + S*Eta;
                B1 = -5.1391118342426808e-16 + 25.0*B1*Xi; S = B1 + S*Eta;
                B0 =  2.8534893807047458e-14 + 24.0*B0*Xi; S = B0 + S*Eta;
                B1 = -1.3923887224181616e-13 + 23.0*B1*Xi; S = B1 + S*Eta;
                B0 =  3.3717632624009806e-13 + 22.0*B0*Xi; S = B0 + S*Eta;
                B1 =  1.1004392031956284e-13 + 21.0*B1*Xi; S = B1 + S*Eta;
                B0 = -5.0276692801141763e-12 + 20.0*B0*Xi; S = B0 + S*Eta;
                B1 =  2.4361948020667402e-11 + 19.0*B1*Xi; S = B1 + S*Eta;
                B0 = -5.8307721325504166e-11 + 18.0*B0*Xi; S = B0 + S*Eta;
                B1 = -2.5514193994946487e-11 + 17.0*B1*Xi; S = B1 + S*Eta;
                B0 =  9.1476995822367933e-10 + 16.0*B0*Xi; S = B0 + S*Eta;
                B1 = -4.3820360184533521e-09 + 15.0*B1*Xi; S = B1 + S*Eta;
                B0 =  1.0261809784240299e-08 + 14.0*B0*Xi; S = B0 + S*Eta;
                B1 =  6.7078535434015332e-09 + 13.0*B1*Xi; S = B1 + S*Eta;
                B0 = -1.7665952736826086e-07 + 12.0*B0*Xi; S = B0 + S*Eta;
                B1 =  8.2967113409530833e-07 + 11.0*B1*Xi; S = B1 + S*Eta;
                B0 = -1.8540622107151585e-06 + 10.0*B0*Xi; S = B0 + S*Eta;
                B1 = -2.1854485106799979e-06 +  9.0*B1*Xi; S = B1 + S*Eta;
                B0 =  3.9192631785224383e-05 +  8.0*B0*Xi; S = B0 + S*Eta;
                B1 = -0.00017875514403292177 +  7.0*B1*Xi; S = B1 + S*Eta;
                B0 =  0.00035273368606701921 +  6.0*B0*Xi; S = B0 + S*Eta;
                B1 =   0.0011574074074074078 +  5.0*B1*Xi; S = B1 + S*Eta;
                B0 =   -0.014814814814814815 +  4.0*B0*Xi; S = B0 + S*Eta;
                B1 =    0.083333333333333329 +  3.0*B1*Xi; S = B1 + S*Eta;
                B0 =    -0.33333333333333331 +  2.0*B0*Xi; S = B0 + S*Eta;
                S  = S / (1.0 + B1*Xi);

                S = S*exp(-0.5*X*Eta*Eta)*rsqrt(2.0*3.141592653589793*X);
                if (X < Lam) {
                    S += 0.5*erfc(Eta*sqrt(0.5*X));
                    if (S > U) X -= 1.0;
                }
                else {
                    S -= 0.5*erfc(-Eta*sqrt(0.5*X));
                    if (S > U-1.0) X -= 1.0;
                }
            }

            // sum downwards or upwards

            else {
                Xi = 1.0 / X;
                S = - (691.0/360360.0);
                S =   (1.0/1188.0) + S*Xi*Xi;
                S = - (1.0/1680.0) + S*Xi*Xi;
                S =   (1.0/1260.0) + S*Xi*Xi;
                S = - (1.0/360.0)  + S*Xi*Xi;
                S =   (1.0/12.0)   + S*Xi*Xi;
                S =                  S*Xi;
                S = (X - Lam) - X*log(X/Lam) - S;

                if (X < Lam) {
                    T  = exp(-0.5*S);
                    S  = 1.0 - T*(U*T) * sqrt(2.0*3.141592653589793*Xi) * Lam;
                    T  = 1.0;
                    Xi = X;
                    for (int i=1; i<50; i++) {
                        Xi -= 1.0;
                        T  *= Xi/Lam;
                        S  += T;
                    }
                    if (S > 0.0) X -= 1.0;
                }

                else {
                    T  = exp(-0.5*S);
                    S  = 1.0 - T*((1.0-U)*T) * sqrt(2.0*3.141592653589793*X);
                    Xi = X;
                    for (int i=0; i<50; i++) {
                        Xi += 1.0;
                        S   = S*Xi/Lam + 1.0;
                    }
                    if (S < 0.0) X -= 1.0;
                }
            }
            return X;
        }
    }

    // bottom-up summation

    double Del;

    X   = 0.0;
    T   = exp(0.5*Lam);
    Del = 0.0;
    if (U > 0.5) Del = T*(1e-13*T);
    S   = 1.0 - T*(U*T) + Del;

    while (S < 0.0) {
        X  += 1.0;
        T   = X/Lam;
        Del = T*Del;
        S   = T*S + 1.0;
    }

    // top-down summation if needed

    if (S < 2.0*Del) {
        Del = 1e13*Del;
        T   = 1e17*Del;
        Del = (1.0-U)*Del;

        while (Del < T) {
            X   += 1.0;
            Del *= X/Lam;
        }

        S = Del;
        T = 1.0;
        while (S > 0.0) {
            T *= X/Lam;
            S -= T;
            X -= 1.0;
        }
    }
    return X;
}

// rlog1(x) = x - log(1+x) for -0.618<x<1.618 using rational minimax approx.
__device__ static inline double rlog1(double X) {

#define P0  0.2000000000000000
#define P1 -0.3636535967811319
#define P2  0.2015244511825799
#define P3 -0.03274937605228191
#define P4  0.00004542775258423288
#define Q1 -2.532553698191349
#define Q2  2.261033627633934
#define Q3 -0.8263420776370621
#define Q4  0.1008870710149048

    double R, T, W;

    // rlog1(x) = f(r) = 2*r^2*(1/(1-r) + r*f1(r))
    // where f1(r) = (log((1+r)/(1-r)) - 2*r)/(2*r^3)
    // which has series expansion sum_{n=0} r^(2n)/(2n + 3)
    // Calculate f1(r) = 1/3 + r^2*P(r^2)/Q(r^2) where P,Q follow
    // from rational minimax approximation for 0 <= r^2 <= 0.2
    R = X/(X + 2.0);
    T = R*R;
    W = (1.0/3.0) + T*((((P4*T + P3)*T + P2)*T + P1)*T + P0 ) /
                      ((((Q4*T + Q3)*T + Q2)*T + Q1)*T + 1.0);
    return T*((X + 2.0) - 2.0*R*W);
} /* rlog1 */
#undef P0
#undef P1
#undef P2
#undef P3
#undef P4
#undef Q1
#undef Q2
#undef Q3
#undef Q4


} /* ----- end namespace RNG ----- */

#endif /* SHARED_CUH_ */
