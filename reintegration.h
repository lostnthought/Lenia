#ifndef REINTEGRATION_H
#define REINTEGRATION_H

#include <xmmintrin.h>
#include <intrin.h>

static inline void hsum(__m128 in, float * out)
{
	__m128 shuf = _mm_shuffle_ps(in, in, _MM_SHUFFLE(2, 3, 0, 1));
	__m128 sums = _mm_add_ps(in, shuf);
	shuf = _mm_shuffle_ps(sums, sums, _MM_SHUFFLE(0, 1, 2, 3));
	sums = _mm_add_ps(sums, shuf);
	_mm_store_ss(out, sums);
}

static inline __m128 fastexp(__m128 x)
{
	const __m128 a0 = _mm_set1_ps(1.0f);
	const __m128 a1 = _mm_set1_ps(1.0f);
	const __m128 a2 = _mm_set1_ps(1.0f / 2);
	const __m128 a3 = _mm_set1_ps(1.0f / 2 / 3);
	const __m128 a4 = _mm_set1_ps(1.0f / 2 / 3 / 4);
	const __m128 a5 = _mm_set1_ps(1.0f / 2 / 3 / 4 / 5);
	const __m128 a6 = _mm_set1_ps(1.0f / 2 / 3 / 4 / 5 / 6);
	const __m128 a7 = _mm_set1_ps(1.0f / 2 / 3 / 4 / 5 / 6 / 7);

	__m128 ret = _mm_add_ps(_mm_mul_ps(a7, x), a6);
	ret = _mm_add_ps(_mm_mul_ps(ret, x), a5);
	ret = _mm_add_ps(_mm_mul_ps(ret, x), a4);
	ret = _mm_add_ps(_mm_mul_ps(ret, x), a3);
	ret = _mm_add_ps(_mm_mul_ps(ret, x), a2);
	ret = _mm_add_ps(_mm_mul_ps(ret, x), a1);
	ret = _mm_add_ps(_mm_mul_ps(ret, x), a0);

	return ret; 

}

// no normalization due to some weirdness with the math currently
static inline void compute_vel(int pcidx, __m128 DAMP)
{
	//static const __m128 LENMAXRG = _mm_set1_ps(0.0f);


	//__m128 len = _mm_sqrt_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(FX, FX), _mm_mul_ps(FY, FY)), _mm_mul_ps(FZ, FZ)));

	
	//__m128 FXM = _mm_cmpgt_ps(FX, LENMAXRG);
	//__m128 DIVX = _mm_and_ps(FXM, len);
	//__m128 FYM = _mm_cmpgt_ps(FY, LENMAXRG);
	//__m128 DIVY = _mm_and_ps(FYM, len);
	//__m128 FZM = _mm_cmpgt_ps(FZ, LENMAXRG);
	//__m128 DIVZ = _mm_and_ps(FZM, len);
	//FX = _mm_div_ps(FX, DIVX);
	//FY = _mm_div_ps(FY, DIVY);
	//FZ = _mm_div_ps(FZ, DIVZ);

	__m128 VX = _mm_load_ps(&V_x[pcidx]);
	__m128 VY = _mm_load_ps(&V_y[pcidx]);
	__m128 VZ = _mm_load_ps(&V_z[pcidx]);
	
	VX = _mm_mul_ps(_mm_add_ps(VX, _mm_load_ps(&F_x[pcidx])), DAMP);
	VY = _mm_mul_ps(_mm_add_ps(VY, _mm_load_ps(&F_y[pcidx])), DAMP);
	VZ = _mm_mul_ps(_mm_add_ps(VZ, _mm_load_ps(&F_z[pcidx])), DAMP);

	_mm_store_ps(&V_x[pcidx], VX);
	_mm_store_ps(&V_y[pcidx], VY);
	_mm_store_ps(&V_z[pcidx], VZ);

	

}	

void precompute_sobel_registers();
void precompute_blur_kernel(Kernel * kernel);
void compute_growth_rt(int section, int end, int c);
void compute_soft_clip_rt(int section, int end, int c);
void compute_growth_tr(int section, int end, int c);
void compute_soft_clip_tr(int section, int end, int c);


void h_sum_blur(TD * t);
// gradient pre-blur
void h_blur(int section, int end);
void sum_blur(int section, int end);
void h_blur_border_wall(int section, int end);
void sum_blur_border_wall(int section, int end);
void h_blur_border_none(int section, int end);
void sum_blur_border_none(int section, int end);

// gradient calc
void h_sobel(int section, int end);
void sum_sobel(int section, int end);
void h_central_diff(int section, int end);
void sum_central_diff(int section, int end);

void compute_velocity_mu_wall(TD * t);
void compute_velocity_mu(TD * t);
void compute_f_mu_wall(TD * t);
void compute_f_mu(TD * t);

// reintegration
void srt(float * ARR, int tid, int section, int end, float clip, float sn);
void rt_border_torus(float * ARR, int tid, int section, int end, float clip, float sn);
void rt_border_other(float * ARR, int tid, int section, int end, float clip, float sn);
void adv_vel(float * ARR, int tid, int section, int end, float clip, float sn);
void adv_f(float * ARR, int tid, int section, int end, float clip, float sn);


#endif