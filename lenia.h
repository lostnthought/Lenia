#ifndef LENIA_H
#define LENIA_H
#include <intrin.h>
#include "render.h"
void init_conway();
void fft_step(TD * t);
void mt_fft_step(TD * t);
void init_kernel(Kernel * kernel);
void init_fftw_plans();
void generate_random();
void init_pos_array();
void init_offset_table(Button * b);
void sobel_convolve(TD * t);
void sobel_convolve_h(TD * t);
void transpose_to_channel_major(TD * t);
void transpose_to_row_major(TD * t);
void compute_f(TD * t);
void compute_mu(TD * t);
void update_velocity(TD * t);
void compute_dpmu();
void reintegration_step(TD * t);
void transpose_sobel_borders(TD * t);
void transpose_borders(TD * t);
void transpose_mesh_borders(TD * t);
void init_gaussian_kernel();
void accumulate_reintegration();
extern float BELL_LUT[BELL_LUT_SIZE];
void init_bell_lut();

// z, y, x technically, we need to load in for simd
static const float sx[3][3][3] =
{
	{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}},
	{{-2, 0, 2}, {-4, 0, 4}, {-2, 0, 2}},
	{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}},
};

static const float sx2[3][3][3] =
{
	{{-1, -2, -1}, {-2, -4, -2}, {-1, -2, -1}},
	{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
	{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}},
};

static const float sy[3][3][3] =
{
	{{-1, -2, -1}, {-2, -4, -2}, {-1, -2, -1}},
	{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
	{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}},
};

static const float sy2[3][3][3] =
{
	{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}},
	{{-2, 0, 2}, {-4, 0, 4}, {-2, 0, 2}},
	{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}}
};

static const float sz[3][3][3] =
{
	{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}},
	{{-2, -4, -2}, {0, 0, 0}, {2, 4, 2}},
	{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}},
};

static const float sz2[3][3][3] =
{
	{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}},
	{{-2, -4, -2}, {0, 0, 0}, {2, 4, 2}},
	{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}},
};

// does not bounds check
static inline int calc_pos(int x, int y, int z)
{
	return (int)(x + y * GRID_WIDTH + z * GRID_WIDTH * GRID_HEIGHT);
}

static inline float bell_lookup(float u, float m, float s)
{

	return BELL_LUT[((int)(m * 255) * 256 * 256) + ((int)(s * 255) * 256) + (int)(u * 255)];

}
static inline float wrap(float x, float b)
{
	x -= floorf(x / b) * b;
	return x;
}

static inline __m128 simd_floor(__m128 x)
{
	__m128 dir = _mm_cvtepi32_ps(_mm_cvttps_epi32(x));
	return _mm_sub_ps(dir, _mm_and_ps(_mm_cmplt_ps(x, dir), _mm_set1_ps(1)));
}


static inline bool central_diff_border_torus(int x, int y, int z, int dx, int dy, int dz, ivec3 * m, ivec3 * i)
{
	//m->x = (x - 1 + GRID_WIDTH) % GRID_WIDTH; 
	//i->x = (x + 1 + GRID_WIDTH) % GRID_WIDTH; 
	//m->y = (y - 1 + GRID_HEIGHT) % GRID_HEIGHT; 
	//i->y = (y + 1 + GRID_HEIGHT) % GRID_HEIGHT; 
	//m->z = (z - 1 + GRID_WIDTH) % GRID_WIDTH;
	//i->z = (z + 1 + GRID_WIDTH) % GRID_WIDTH;
}

static inline ivec3 central_diff_border_wall(int x, int y, int z, int dx, int dy, int dz, ivec3 * m, ivec3 * i)
{
	//i->x = x + 1; i->y = y + 1; iz = z + 1;
	//if (ix >= GRID_WIDTH) ix = x;
	//if (iy >= GRID_HEIGHT) iy = y;
	//if (iz >= GRID_WIDTH) iz = z;
	//mx = x - 1; my = y - 1; mz = z - 1;
	//if (mx < 0) mx = x;
	//if (my < 0) my = y;
	//if (mz < 0) mz = z;
}

static inline float alpha_padded_fn(float * a, int x, int y, int z, int c)
{
	return clamp(powf(a[CPIDX(PIDX(x, y, z), c)] / THETA_A, N), 0, 1);
}
static inline float alpha_unpadded_fn(float * a, int x, int y, int z, int c)
{
	return clamp(powf(a[CIDX(IDX(x, y, z), c)] / THETA_A, N), 0, 1);
}
//
//static inline float alpha_grid(float a)
//{
//	return clamp(powf(a / THETA_A, N), 0, 1);
//}

#endif