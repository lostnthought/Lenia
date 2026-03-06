#include "render.h"
#include "lenia.h"
#include "reintegration.h"
#include <xmmintrin.h>

__declspec(align(16)) fftwf_complex H_PREBLUR[MAX_CHANNELS][GRID_SIZE];
__declspec(align(16)) fftwf_complex H_BLUR[MAX_CHANNELS][GRID_SIZE];
__declspec(align(16)) fftwf_complex SUM_PREBLUR[GRID_SIZE];
__declspec(align(16)) fftwf_complex SUM_BLUR[GRID_SIZE];
__declspec(align(16)) fftwf_complex TEMP_H_BLUR[MAX_CHANNELS][GRID_SIZE];
__declspec(align(16)) fftwf_complex TEMP_BLUR[GRID_SIZE];

void precompute_blur_kernel(Kernel * kernel)
{
	memset(kernel->kg, 0, sizeof(fftwf_complex) * GRID_SIZE);
	memset(kernel->real_k, 0, sizeof(float) * GRID_SIZE);
	memset(kernel->grid, 0, sizeof(fftwf_complex) * GRID_SIZE);
	int r = kernel->radius;
	float k_cells = 0.0f;
	int d = r * 2 + 1;
	int bound = GRID_WIDTH;
	int offset_w = GRID_WIDTH / 2;
	int offset_h = GRID_HEIGHT / 2;

	int beta_size = kernel->beta_size;
	const float s = 1.5f;
	/*for (float dx = -1; dx <= 1; dx++)
	for (float dy = -1; dy <= 1; dy++)
	for (float dz = -1; dz <= 1; dz++)
	{
		float d2 = dx * dx + dy * dy + dz * dz;
		gaussian_kernel[(int)(dx + 1)][(int)(dy + 1)][(int)(dz + 1)] = exp(-d2 / (s * s * s));
		int o = 1;
	}*/
	
	int abs_x = 0, abs_y = 0, abs_z = 0;
	for (int r_x = -2; r_x <= 2; r_x++)
	{
	abs_y = 0;
	for (int r_y = -2; r_y <= 2; r_y++)
	{
	abs_z = 0;
	for (int r_z = -2; r_z <= 2; r_z++)
	{
		float x = fabsf((float)r_x);
		float y = fabsf((float)r_y);
		float z = fabsf((float)r_z);
		
		float dist = (x * x + y * y + z * z);
		float g = expf(-dist / (s * s * s));


		int real_x = emod(r_x, GRID_WIDTH);
		int real_y = emod(r_y, GRID_HEIGHT);
		int real_z = emod(r_z, GRID_WIDTH);


		kernel->grid[(real_x) + (real_y) * GRID_WIDTH + (real_z) * GRID_WIDTH * GRID_HEIGHT][0] = g;
		kernel->real_k[r_x + offset_w + (r_y + offset_h) * GRID_WIDTH + (r_z + offset_w) * GRID_WIDTH * GRID_HEIGHT] = g; 
		kernel->grid[(real_x) + (real_y) * GRID_WIDTH + (real_z) * GRID_WIDTH * GRID_HEIGHT][1] = 0;

		k_cells += g;

		abs_z++;
	}
	}
	}
	for (int i = 0; i < GRID_SIZE; i++)
	{
		kernel->grid[i][0] /= k_cells;
		kernel->real_k[i] /= k_cells;
		kernel->grid[i][1] = 0;
	}


	fftwf_execute_dft(forward_fft, kernel->grid, kernel->kg);

}


// transpose h and sum into fftwf complex
// forward fft h and sums (max 4 ffts)
// complex multiply h and sums
// reverse fft (max 4 ffts)
// normalize



void h_sum_blur(TD * t)
{

	if (BLUR_CENTRAL_DIFF)
	{
	int csection = (channels / SIM_THREADS) * t->id;
	int cend = (t->id == SIM_THREADS - 1 ? channels : csection + channels / SIM_THREADS);
	//int csection = 0; int cend = channels;

	for (int c = csection; c < cend; c++){
		int coffset = CPOFFSET(c);
		for (int z = 0; z < GRID_WIDTH; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = IDX(x, y, z);
			H_PREBLUR[c][idx][0] = H[coffset + PIDX(x, y, z)];
			H_PREBLUR[c][idx][1] = 0;
		}
	}

	if (t->id == 0){

		for (int z = 0; z < GRID_WIDTH; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = IDX(x, y, z);
			SUM_PREBLUR[idx][0] = SUM[PIDX(x, y, z)];
			SUM_PREBLUR[idx][1] = 0;
		}
	}

	for (int c = csection; c < cend; c++)
	{
		fftwf_execute_dft(forward_fft, H_PREBLUR[c], H_BLUR[c]);
	}
	if (t->id == 0)
	{
		fftwf_execute_dft(forward_fft, SUM_PREBLUR, SUM_BLUR);
	}

	for (int c = csection; c < cend; c++)
	{
		for (int i = 0; i < GRID_SIZE; i++)
		{
			
			TEMP_H_BLUR[c][i][0] = 
				H_BLUR[c][i][0] * blur_kernel.kg[i][0] - H_BLUR[c][i][1] * blur_kernel.kg[i][1];
			TEMP_H_BLUR[c][i][1] = 
				H_BLUR[c][i][0] * blur_kernel.kg[i][1] + H_BLUR[c][i][1] * blur_kernel.kg[i][0];
		}
		fftwf_execute_dft(reverse_fft, TEMP_H_BLUR[c], H_BLUR[c]);
	}
	if (t->id == 0)
	{
		for (int i = 0; i < GRID_SIZE; i++)
		{

			TEMP_BLUR[i][0] = 
				SUM_BLUR[i][0] * blur_kernel.kg[i][0] - SUM_BLUR[i][1] * blur_kernel.kg[i][1];
			TEMP_BLUR[i][1] = 
				SUM_BLUR[i][0] * blur_kernel.kg[i][1] + SUM_BLUR[i][1] * blur_kernel.kg[i][0];
		}
		fftwf_execute_dft(reverse_fft, TEMP_BLUR, SUM_BLUR);
	}

	for (int c = csection; c < cend; c++){
		int cpoffset = CPOFFSET(c);
		for (int z = 0; z < GRID_WIDTH; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int pidx = PIDX(x, y, z);
			int idx = IDX(x, y, z);
			U_S[cpoffset + pidx] = H_BLUR[c][idx][0] / GRID_SIZE;
		}
	}
	if (t->id == 0)
	{
		for (int z = 0; z < GRID_WIDTH; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = IDX(x, y, z);
			int pidx = PIDX(x, y, z);
			SUM_S[pidx] = SUM_BLUR[idx][0] / GRID_SIZE;
		}
	}
	}
	else
	{
		int section = (GRID_PADDED / SIM_THREADS) * t->id;
		int end = (t->id == SIM_THREADS - 1 ? GRID_PADDED : section + GRID_PADDED / SIM_THREADS);
		memcpy(&U_S[section], &H[section], sizeof(float) * (end - section) * MAX_CHANNELS);
		memcpy(&SUM_S[section], &SUM[section], sizeof(float) * (end - section) * MAX_CHANNELS);

	}

}


void h_sobel(int section, int end)
{
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		int coffset = COFFSET(c);
		for (int dz = -1; dz <= 1; dz++)
		for (int dy = -1; dy <= 1; dy++)
		for (int dx = -1; dx <= 1; dx++)
		{
			__m128 SX1 = _mm_set1_ps(sx[dx+1][dy+1][dz+1]);
			__m128 SY1 = _mm_set1_ps(sy[dx+1][dy+1][dz+1]);
			__m128 SZ1 = _mm_set1_ps(sz[dx+1][dy+1][dz+1]);

			for (int z = section; z < end; z++)
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x+=4)
			{
				int nidx = IDX(x, y, z);
		
				int ncidx = nidx + coffset;	

				int cidx = PIDX(x + dx, y + dy, z + dz) + cpoffset;

				__m128 h = _mm_loadu_ps(&H[cidx]);
				__m128 sumx = _mm_add_ps(_mm_load_ps(&NU_x[ncidx]), _mm_mul_ps(h, SX1));
				__m128 sumy = _mm_add_ps(_mm_load_ps(&NU_y[ncidx]), _mm_mul_ps(h, SY1));
				__m128 sumz = _mm_add_ps(_mm_load_ps(&NU_z[ncidx]), _mm_mul_ps(h, SZ1));
				_mm_store_ps(&NU_x[ncidx], sumx);
				_mm_store_ps(&NU_y[ncidx], sumy);
				_mm_store_ps(&NU_z[ncidx], sumz);


			}
		}
	}
}

void sum_sobel(int section, int end)
{
	for (int dz = -1; dz <= 1; dz++)
	for (int dy = -1; dy <= 1; dy++)
	for (int dx = -1; dx <= 1; dx++)
	{
		__m128 SX1 = _mm_set1_ps(sx[dx+1][dy+1][dz+1]);
		__m128 SY1 = _mm_set1_ps(sy[dx+1][dy+1][dz+1]);
		__m128 SZ1 = _mm_set1_ps(sz[dx+1][dy+1][dz+1]);

		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{
			int nidx = IDX(x, y, z);

			int cidx = PIDX(x + dx, y + dy, z + dz);

			__m128 sum = _mm_loadu_ps(&SUM[cidx]);
			__m128 sumx = _mm_add_ps(_mm_load_ps(&NA_x[nidx]), _mm_mul_ps(sum, SX1));
			__m128 sumy = _mm_add_ps(_mm_load_ps(&NA_y[nidx]), _mm_mul_ps(sum, SY1));
			__m128 sumz = _mm_add_ps(_mm_load_ps(&NA_z[nidx]), _mm_mul_ps(sum, SZ1));
			_mm_store_ps(&NA_x[nidx], sumx);
			_mm_store_ps(&NA_y[nidx], sumy);
			_mm_store_ps(&NA_z[nidx], sumz);


		}
	}
	
}
void h_central_diff(int section, int end)
{
	static const __m128 PFRG = _mm_set1_ps(0.5f);
	for (int c = 0; c < channels; c++)
	{
		int coffset = COFFSET(c);
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{

			int nidx = IDX(x, y, z);
			int mx = x - 1, ix = x + 1, my = y - 1, iy = y + 1, mz = z - 1, iz = z + 1;
			int cidx = nidx + coffset;

			__m128 MM = _mm_loadu_ps(&U_S[PIDX(mx, y, z) + cpoffset]);
			__m128 II = _mm_loadu_ps(&U_S[PIDX(ix, y, z) + cpoffset]);
			_mm_store_ps(&NU_x[cidx], _mm_mul_ps(_mm_sub_ps(II, MM), PFRG));

			MM = _mm_loadu_ps(&U_S[PIDX(x, my, z) + cpoffset]);
			II = _mm_loadu_ps(&U_S[PIDX(x, iy, z) + cpoffset]);
			_mm_store_ps(&NU_y[cidx], _mm_mul_ps(_mm_sub_ps(II, MM), PFRG));

			MM = _mm_loadu_ps(&U_S[PIDX(x, y, mz) + cpoffset]);
			II = _mm_loadu_ps(&U_S[PIDX(x, y, iz) + cpoffset]);
			_mm_store_ps(&NU_z[cidx], _mm_mul_ps(_mm_sub_ps(II, MM), PFRG));


		}
	}
}
void sum_central_diff(int section, int end)
{
	static const __m128 PFRG = _mm_set1_ps(0.5f);
	for (int z = section; z < end; z++)
	for (int y = 0; y < GRID_HEIGHT; y++)
	for (int x = 0; x < GRID_WIDTH; x+=4)
	{
		int nidx = IDX(x, y, z);
		int mx = x - 1, ix = x + 1, my = y - 1, iy = y + 1, mz = z - 1, iz = z + 1;

		__m128 MM = _mm_loadu_ps(&SUM_S[PIDX(mx, y, z)]);
		__m128 II = _mm_loadu_ps(&SUM_S[PIDX(ix, y, z)]);
		_mm_store_ps(&NA_x[nidx], _mm_mul_ps(_mm_sub_ps(II, MM), PFRG));

		MM = _mm_loadu_ps(&SUM_S[PIDX(x, my, z)]);
		II = _mm_loadu_ps(&SUM_S[PIDX(x, iy, z)]);
		_mm_store_ps(&NA_y[nidx], _mm_mul_ps(_mm_sub_ps(II, MM), PFRG));

		MM = _mm_loadu_ps(&SUM_S[PIDX(x, y, mz)]);
		II = _mm_loadu_ps(&SUM_S[PIDX(x, y, iz)]);
		_mm_store_ps(&NA_z[nidx], _mm_mul_ps(_mm_sub_ps(II, MM), PFRG));

	}
}





void compute_velocity_mu_wall(TD * t)
{

	float sigma = off_table.sigma;
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);

	float GWS = GRID_WIDTH - sigma;
	float GHS = GRID_HEIGHT - sigma;

	__m128 DAMP = _mm_set1_ps(V_DAMP);

	for (int c = 0; c < channels; c++)
	{
		int coffset = COFFSET(c);
		int cpoffset = CPOFFSET(c);

		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{

			int idx = IDX(x, y, z);
			int pidx = PIDX(x, y, z);

			int cidx = idx + coffset;
			int pcidx = pidx + cpoffset;
			compute_vel(pcidx, DAMP);
			{
				__m128 VX = _mm_load_ps(&V_x[pcidx]);
				__m128 VY = _mm_load_ps(&V_y[pcidx]);
				__m128 VZ = _mm_load_ps(&V_z[pcidx]);
				{
					__m128 TIMERG = _mm_set1_ps(timestep);
					VX = _mm_mul_ps(VX, TIMERG);
					VY = _mm_mul_ps(VY, TIMERG);
					VZ = _mm_mul_ps(VZ, TIMERG);
				}
				{
					__m128 NMARG = _mm_set1_ps(-off_table.ma);
					VX = _mm_max_ps(VX, NMARG);
					VY = _mm_max_ps(VY, NMARG);
					VZ = _mm_max_ps(VZ, NMARG);
				}
				{
					__m128 MARG = _mm_set1_ps(off_table.ma);
					VX = _mm_min_ps(VX, MARG);
					VY = _mm_min_ps(VY, MARG);
					VZ = _mm_min_ps(VZ, MARG);
				}

				VX = _mm_add_ps(_mm_load_ps(&POS_ARRAY_X[idx]), VX);
				VY = _mm_add_ps(_mm_load_ps(&POS_ARRAY_Y[idx]), VY);
				VZ = _mm_add_ps(_mm_load_ps(&POS_ARRAY_Z[idx]), VZ);

				{
					__m128 SIGMARG = _mm_set1_ps(sigma);
					VX = _mm_max_ps(VX, SIGMARG);
					VY = _mm_max_ps(VY, SIGMARG);
					VZ = _mm_max_ps(VZ, SIGMARG);
				}
				{
					__m128 SWRG = _mm_set1_ps(GWS);
					__m128 SHRG = _mm_set1_ps(GHS);
					VX = _mm_min_ps(VX, SWRG);
					VY = _mm_min_ps(VY, SHRG);
					VZ = _mm_min_ps(VZ, SWRG);
				}

				_mm_store_ps(&MU_x[pcidx], VX);
				_mm_store_ps(&MU_y[pcidx], VY);
				_mm_store_ps(&MU_z[pcidx], VZ);
			}
			
		}

	}

	
}

void compute_velocity_mu(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	__m128 DAMP = _mm_set1_ps(V_DAMP);

	for (int c = 0; c < channels; c++)
	{
		int coffset = COFFSET(c);
		int cpoffset = CPOFFSET(c);

		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{

			int idx = IDX(x, y, z);
			int pidx = PIDX(x, y, z);

	
			int cidx = idx + coffset;
			int pcidx = pidx + cpoffset;
			compute_vel(pcidx, DAMP);
			{
				__m128 VX = _mm_load_ps(&V_x[pcidx]);
				__m128 VY = _mm_load_ps(&V_y[pcidx]);
				__m128 VZ = _mm_load_ps(&V_z[pcidx]);
				{
					__m128 TIMERG = _mm_set1_ps(timestep);
					VX = _mm_mul_ps(VX, TIMERG);
					VY = _mm_mul_ps(VY, TIMERG);
					VZ = _mm_mul_ps(VZ, TIMERG);
				}
				{
					__m128 NMARG = _mm_set1_ps(-off_table.ma);
					VX = _mm_max_ps(VX, NMARG);
					VY = _mm_max_ps(VY, NMARG);
					VZ = _mm_max_ps(VZ, NMARG);
				}
				{
					__m128 MARG = _mm_set1_ps(off_table.ma);
					VX = _mm_min_ps(VX, MARG);
					VY = _mm_min_ps(VY, MARG);
					VZ = _mm_min_ps(VZ, MARG);
				}

				VX = _mm_add_ps(_mm_load_ps(&POS_ARRAY_X[idx]), VX);
				VY = _mm_add_ps(_mm_load_ps(&POS_ARRAY_Y[idx]), VY);
				VZ = _mm_add_ps(_mm_load_ps(&POS_ARRAY_Z[idx]), VZ);


				_mm_store_ps(&MU_x[pcidx], VX);
				_mm_store_ps(&MU_y[pcidx], VY);
				_mm_store_ps(&MU_z[pcidx], VZ);
			}
			
		}

	}
}

void compute_f_mu_wall(TD * t)
{
	float sigma = off_table.sigma;
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);

	float GWS = GRID_WIDTH - sigma;
	float GHS = GRID_HEIGHT - sigma;


	__m128 DAMP = _mm_set1_ps(V_DAMP);

	for (int c = 0; c < channels; c++)
	{
		int coffset = COFFSET(c);
		int cpoffset = CPOFFSET(c);

		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{

			int idx = IDX(x, y, z);
			int pidx = PIDX(x, y, z);

	
			int cidx = idx + coffset;
			int pcidx = pidx + cpoffset;
			{
				__m128 VX = _mm_load_ps(&F_x[pcidx]);
				__m128 VY = _mm_load_ps(&F_y[pcidx]);
				__m128 VZ = _mm_load_ps(&F_z[pcidx]);
				{
					__m128 TIMERG = _mm_set1_ps(timestep);
					VX = _mm_mul_ps(VX, TIMERG);
					VY = _mm_mul_ps(VY, TIMERG);
					VZ = _mm_mul_ps(VZ, TIMERG);
				}
				{
					__m128 NMARG = _mm_set1_ps(-off_table.ma);
					VX = _mm_max_ps(VX, NMARG);
					VY = _mm_max_ps(VY, NMARG);
					VZ = _mm_max_ps(VZ, NMARG);
				}
				{
					__m128 MARG = _mm_set1_ps(off_table.ma);
					VX = _mm_min_ps(VX, MARG);
					VY = _mm_min_ps(VY, MARG);
					VZ = _mm_min_ps(VZ, MARG);
				}

				VX = _mm_add_ps(_mm_load_ps(&POS_ARRAY_X[idx]), VX);
				VY = _mm_add_ps(_mm_load_ps(&POS_ARRAY_Y[idx]), VY);
				VZ = _mm_add_ps(_mm_load_ps(&POS_ARRAY_Z[idx]), VZ);


				{
					__m128 SIGMARG = _mm_set1_ps(sigma);
					VX = _mm_max_ps(VX, SIGMARG);
					VY = _mm_max_ps(VY, SIGMARG);
					VZ = _mm_max_ps(VZ, SIGMARG);
				}
				{
					__m128 SWRG = _mm_set1_ps(GWS);
					__m128 SHRG = _mm_set1_ps(GHS);
					VX = _mm_min_ps(VX, SWRG);
					VY = _mm_min_ps(VY, SHRG);
					VZ = _mm_min_ps(VZ, SWRG);
				}


				_mm_store_ps(&MU_x[pcidx], VX);
				_mm_store_ps(&MU_y[pcidx], VY);
				_mm_store_ps(&MU_z[pcidx], VZ);
			}
			compute_vel(pcidx, DAMP);
		}

	}
}

void compute_f_mu(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	__m128 DAMP = _mm_set1_ps(V_DAMP);

	for (int c = 0; c < channels; c++)
	{
		int coffset = COFFSET(c);
		int cpoffset = CPOFFSET(c);

		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{

			int idx = IDX(x, y, z);
			int pidx = PIDX(x, y, z);

	
			int cidx = idx + coffset;
			int pcidx = pidx + cpoffset;
			{
				__m128 VX = _mm_load_ps(&F_x[pcidx]);
				__m128 VY = _mm_load_ps(&F_y[pcidx]);
				__m128 VZ = _mm_load_ps(&F_z[pcidx]);
				{
					__m128 TIMERG = _mm_set1_ps(timestep);
					VX = _mm_mul_ps(VX, TIMERG);
					VY = _mm_mul_ps(VY, TIMERG);
					VZ = _mm_mul_ps(VZ, TIMERG);
				}
				{
					__m128 NMARG = _mm_set1_ps(-off_table.ma);
					VX = _mm_max_ps(VX, NMARG);
					VY = _mm_max_ps(VY, NMARG);
					VZ = _mm_max_ps(VZ, NMARG);
				}
				{
					__m128 MARG = _mm_set1_ps(off_table.ma);
					VX = _mm_min_ps(VX, MARG);
					VY = _mm_min_ps(VY, MARG);
					VZ = _mm_min_ps(VZ, MARG);
				}

				VX = _mm_add_ps(_mm_load_ps(&POS_ARRAY_X[idx]), VX);
				VY = _mm_add_ps(_mm_load_ps(&POS_ARRAY_Y[idx]), VY);
				VZ = _mm_add_ps(_mm_load_ps(&POS_ARRAY_Z[idx]), VZ);


				_mm_store_ps(&MU_x[pcidx], VX);
				_mm_store_ps(&MU_y[pcidx], VY);
				_mm_store_ps(&MU_z[pcidx], VZ);
			}
			compute_vel(pcidx, DAMP);
		}

	}
}




void srt(float * ARR, int tid, int section, int end, float clip, float sn)
{

	__m128 si = _mm_set1_ps(0.5f + off_table.sigma);
	__m128 zero = _mm_set1_ps(0.0f);
	__m128 cl = _mm_set1_ps(clip);
	__m128 sns = _mm_set1_ps(sn);
	for (int z = section; z < end; z++)
	for (int y = 0; y < GRID_HEIGHT; y++)
	for (int x = 0; x < GRID_WIDTH; x+=4)
	{

		int pidx = PIDX(x, y, z);
		int idx = IDX(x, y, z);
			
		int c = 0;
		int cidx = CPIDX(pidx, c);

		__m128 mx = abs_ps(_mm_sub_ps(_mm_load_ps(&POS_ARRAY_X[idx]), _mm_load_ps(&MU_x[cidx])));
		__m128 my = abs_ps(_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Y[idx]), _mm_load_ps(&MU_y[cidx])));
		__m128 mz = abs_ps(_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Z[idx]), _mm_load_ps(&MU_z[cidx])));
		mx = _mm_sub_ps(si, mx);
		my = _mm_sub_ps(si, my);
		mz = _mm_sub_ps(si, mz);
		mx = _mm_min_ps(_mm_max_ps(mx, zero), cl);
		my = _mm_min_ps(_mm_max_ps(my, zero), cl);
		mz = _mm_min_ps(_mm_max_ps(mz, zero), cl);
		__m128 w = _mm_mul_ps(_mm_load_ps(&ARR[cidx]), mx);
		w = _mm_mul_ps(w, my);
		w = _mm_mul_ps(w, mz);
		w = _mm_mul_ps(w, sns);
		w = _mm_add_ps(_mm_load_ps(&TT[cidx]), w);

		_mm_store_ps(&TT[cidx], w);

		if (channels == 1) continue;

		c++;
		cidx = CPIDX(pidx, c);

		mx = abs_ps(_mm_sub_ps(_mm_load_ps(&POS_ARRAY_X[idx]), _mm_load_ps(&MU_x[cidx])));
		my = abs_ps(_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Y[idx]), _mm_load_ps(&MU_y[cidx])));
		mz = abs_ps(_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Z[idx]), _mm_load_ps(&MU_z[cidx])));
		mx = _mm_sub_ps(si, mx);
		my = _mm_sub_ps(si, my);
		mz = _mm_sub_ps(si, mz);
		mx = _mm_min_ps(_mm_max_ps(mx, zero), cl);
		my = _mm_min_ps(_mm_max_ps(my, zero), cl);
		mz = _mm_min_ps(_mm_max_ps(mz, zero), cl);
		w = _mm_mul_ps(_mm_load_ps(&ARR[cidx]), mx);
		w = _mm_mul_ps(w, my);
		w = _mm_mul_ps(w, mz);
		w = _mm_mul_ps(w, sns);
		w = _mm_add_ps(_mm_load_ps(&TT[cidx]), w);

		_mm_store_ps(&TT[cidx], w);

		if (channels == 2) continue;

		c++;
		cidx = CPIDX(pidx, c);

		mx = abs_ps(_mm_sub_ps(_mm_load_ps(&POS_ARRAY_X[idx]), _mm_load_ps(&MU_x[cidx])));
		my = abs_ps(_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Y[idx]), _mm_load_ps(&MU_y[cidx])));
		mz = abs_ps(_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Z[idx]), _mm_load_ps(&MU_z[cidx])));
		mx = _mm_sub_ps(si, mx);
		my = _mm_sub_ps(si, my);
		mz = _mm_sub_ps(si, mz);
		mx = _mm_min_ps(_mm_max_ps(mx, zero), cl);
		my = _mm_min_ps(_mm_max_ps(my, zero), cl);
		mz = _mm_min_ps(_mm_max_ps(mz, zero), cl);
		w = _mm_mul_ps(_mm_load_ps(&ARR[cidx]), mx);
		w = _mm_mul_ps(w, my);
		w = _mm_mul_ps(w, mz);
		w = _mm_mul_ps(w, sns);
		w = _mm_add_ps(_mm_load_ps(&TT[cidx]), w);

		_mm_store_ps(&TT[cidx], w);
		
	}

}


void rt_border_torus(float * ARR, int tid, int section, int end, float clip, float sn)
{

	int dd = (int)off_table.dd;
	
	for (int dz = -dd; dz <= dd; dz++)
	for (int dy = -dd; dy <= dd; dy++)
	for (int dx = -dd; dx <= dd; dx++)
	{
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{

			int rx = x + dx;
			int ry = y + dy;
			int rz = z + dz;

			int idx = IDX(x, y, z);
			int dpidx = PIDX(x, y, z);

			__m128 pos_x = _mm_load_ps(&POS_ARRAY_X[idx]);
			__m128 pos_y = _mm_load_ps(&POS_ARRAY_Y[idx]);
			__m128 pos_z = _mm_load_ps(&POS_ARRAY_Z[idx]);
		

			int spidx = PIDX(rx, ry, rz);

			int c = 0;	
			int dcidx = CPIDX(dpidx, c);

			int scidx = CPIDX(spidx, c);

			__m128 mx = _mm_loadu_ps(&MU_x[scidx]);
			__m128 my = _mm_loadu_ps(&MU_y[scidx]);
			__m128 mz = _mm_loadu_ps(&MU_z[scidx]);

			mx = _mm_sub_ps(pos_x, mx);
			my = _mm_sub_ps(pos_y, my);
			mz = _mm_sub_ps(pos_z, mz);

			{
				__m128 cg = _mm_cmpgt_ps(mx, _mm_set1_ps(HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mx = _mm_sub_ps(mx, cg);

				cg = _mm_cmplt_ps(mx, _mm_set1_ps(-HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mx = _mm_add_ps(mx, cg);

				cg = _mm_cmpgt_ps(my, _mm_set1_ps(HALF_HEIGHT));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_HEIGHT));
				my = _mm_sub_ps(my, cg);

				cg = _mm_cmplt_ps(my, _mm_set1_ps(-HALF_HEIGHT));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_HEIGHT));
				my = _mm_add_ps(my, cg);

				cg = _mm_cmpgt_ps(mz, _mm_set1_ps(HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mz = _mm_sub_ps(mz, cg);

				cg = _mm_cmplt_ps(mz, _mm_set1_ps(-HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mz = _mm_add_ps(mz, cg);
			}

			mx = abs_ps(mx);
			my = abs_ps(my);
			mz = abs_ps(mz);

			{
				__m128 si = _mm_set1_ps(0.5f + off_table.sigma);
				mx = _mm_sub_ps(si, mx);
				my = _mm_sub_ps(si, my);
				mz = _mm_sub_ps(si, mz);
			}
			{
				__m128 ZRG = _mm_set1_ps(0.0f);
				mx = _mm_max_ps(mx, ZRG);
				my = _mm_max_ps(my, ZRG);
				mz = _mm_max_ps(mz, ZRG);
			}
			{
				__m128 CLRG = _mm_set1_ps(clip);
				mx = _mm_min_ps(mx, CLRG);
				my = _mm_min_ps(my, CLRG);
				mz = _mm_min_ps(mz, CLRG);
			}

			{
				__m128 w = _mm_mul_ps(_mm_loadu_ps(&ARR[scidx]), mx);
				w = _mm_mul_ps(w, my);
				w = _mm_mul_ps(w, mz);
				w = _mm_mul_ps(w, _mm_set1_ps(sn));


				_mm_store_ps(&TT[dcidx], _mm_add_ps(_mm_load_ps(&TT[dcidx]), w));
			}

			if (channels == 1) continue;

			c++;

			dcidx = CPIDX(dpidx, c);

			scidx = CPIDX(spidx, c);

			mx = _mm_loadu_ps(&MU_x[scidx]);
			my = _mm_loadu_ps(&MU_y[scidx]);
			mz = _mm_loadu_ps(&MU_z[scidx]);

			mx = _mm_sub_ps(pos_x, mx);
			my = _mm_sub_ps(pos_y, my);
			mz = _mm_sub_ps(pos_z, mz);

			{
				__m128 cg = _mm_cmpgt_ps(mx, _mm_set1_ps(HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mx = _mm_sub_ps(mx, cg);

				cg = _mm_cmplt_ps(mx, _mm_set1_ps(-HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mx = _mm_add_ps(mx, cg);

				cg = _mm_cmpgt_ps(my, _mm_set1_ps(HALF_HEIGHT));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_HEIGHT));
				my = _mm_sub_ps(my, cg);

				cg = _mm_cmplt_ps(my, _mm_set1_ps(-HALF_HEIGHT));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_HEIGHT));
				my = _mm_add_ps(my, cg);

				cg = _mm_cmpgt_ps(mz, _mm_set1_ps(HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mz = _mm_sub_ps(mz, cg);

				cg = _mm_cmplt_ps(mz, _mm_set1_ps(-HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mz = _mm_add_ps(mz, cg);
			}

			mx = abs_ps(mx);
			my = abs_ps(my);
			mz = abs_ps(mz);

			{
				__m128 si = _mm_set1_ps(0.5f + off_table.sigma);
				mx = _mm_sub_ps(si, mx);
				my = _mm_sub_ps(si, my);
				mz = _mm_sub_ps(si, mz);
			}
			{
				__m128 ZRG = _mm_set1_ps(0.0f);
				mx = _mm_max_ps(mx, ZRG);
				my = _mm_max_ps(my, ZRG);
				mz = _mm_max_ps(mz, ZRG);
			}
			{
				__m128 CLRG = _mm_set1_ps(clip);
				mx = _mm_min_ps(mx, CLRG);
				my = _mm_min_ps(my, CLRG);
				mz = _mm_min_ps(mz, CLRG);
			}

			{
				__m128 w = _mm_mul_ps(_mm_loadu_ps(&ARR[scidx]), mx);
				w = _mm_mul_ps(w, my);
				w = _mm_mul_ps(w, mz);
				w = _mm_mul_ps(w, _mm_set1_ps(sn));


				_mm_store_ps(&TT[dcidx], _mm_add_ps(_mm_load_ps(&TT[dcidx]), w));
			}

			if (channels == 2) continue;
			
			c++;
			dcidx = CPIDX(dpidx, c);

			scidx = CPIDX(spidx, c);
			mx = _mm_loadu_ps(&MU_x[scidx]);
			my = _mm_loadu_ps(&MU_y[scidx]);
			mz = _mm_loadu_ps(&MU_z[scidx]);

			mx = _mm_sub_ps(pos_x, mx);
			my = _mm_sub_ps(pos_y, my);
			mz = _mm_sub_ps(pos_z, mz);
			{
				__m128 cg = _mm_cmpgt_ps(mx, _mm_set1_ps(HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mx = _mm_sub_ps(mx, cg);

				cg = _mm_cmplt_ps(mx, _mm_set1_ps(-HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mx = _mm_add_ps(mx, cg);

				cg = _mm_cmpgt_ps(my, _mm_set1_ps(HALF_HEIGHT));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_HEIGHT));
				my = _mm_sub_ps(my, cg);

				cg = _mm_cmplt_ps(my, _mm_set1_ps(-HALF_HEIGHT));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_HEIGHT));
				my = _mm_add_ps(my, cg);

				cg = _mm_cmpgt_ps(mz, _mm_set1_ps(HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mz = _mm_sub_ps(mz, cg);

				cg = _mm_cmplt_ps(mz, _mm_set1_ps(-HALF_WIDTH));
				cg = _mm_and_ps(cg, _mm_set1_ps(GRID_WIDTH));
				mz = _mm_add_ps(mz, cg);
			}

			mx = abs_ps(mx);
			my = abs_ps(my);
			mz = abs_ps(mz);

			{
				__m128 si = _mm_set1_ps(0.5f + off_table.sigma);
				mx = _mm_sub_ps(si, mx);
				my = _mm_sub_ps(si, my);
				mz = _mm_sub_ps(si, mz);
			}
			{
				__m128 ZRG = _mm_set1_ps(0.0f);
				mx = _mm_max_ps(mx, ZRG);
				my = _mm_max_ps(my, ZRG);
				mz = _mm_max_ps(mz, ZRG);
			}
			{
				__m128 CLRG = _mm_set1_ps(clip);
				mx = _mm_min_ps(mx, CLRG);
				my = _mm_min_ps(my, CLRG);
				mz = _mm_min_ps(mz, CLRG);
			}

			{
				__m128 w = _mm_mul_ps(_mm_loadu_ps(&ARR[scidx]), mx);
				w = _mm_mul_ps(w, my);
				w = _mm_mul_ps(w, mz);
				w = _mm_mul_ps(w, _mm_set1_ps(sn));


				_mm_store_ps(&TT[dcidx], _mm_add_ps(_mm_load_ps(&TT[dcidx]), w));
			}

			
		}
	}
}

void rt_border_other(float * ARR, int tid, int section, int end, float clip, float sn)
{

	int dd = (int)off_table.dd;
	for (int dz = -dd; dz <= dd; dz++)
	for (int dy = -dd; dy <= dd; dy++)
	for (int dx = -dd; dx <= dd; dx++)
	{
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{

			int idx = IDX(x, y, z);
			int dpidx = PIDX(x, y, z);
			__m128 pos_x = _mm_load_ps(&POS_ARRAY_X[idx]);
			__m128 pos_y = _mm_load_ps(&POS_ARRAY_Y[idx]);
			__m128 pos_z = _mm_load_ps(&POS_ARRAY_Z[idx]);
		


			int rx = x + dx;
			int ry = y + dy;
			int rz = z + dz;

			int spidx = PIDX(rx, ry, rz);

			int c = 0;
			int dcidx = CPIDX(dpidx, c);

			int scidx = CPIDX(spidx, c);

			__m128 mx = _mm_loadu_ps(&MU_x[scidx]);
			__m128 my = _mm_loadu_ps(&MU_y[scidx]);
			__m128 mz = _mm_loadu_ps(&MU_z[scidx]);

			mx = _mm_sub_ps(pos_x, mx);
			my = _mm_sub_ps(pos_y, my);
			mz = _mm_sub_ps(pos_z, mz);
			
			mx = abs_ps(mx);
			my = abs_ps(my);
			mz = abs_ps(mz);
			{
				__m128 si = _mm_set1_ps(0.5f + off_table.sigma);
				mx = _mm_sub_ps(si, mx);
				my = _mm_sub_ps(si, my);
				mz = _mm_sub_ps(si, mz);
			}
			{
				__m128 ZRG = _mm_set1_ps(0.0f);
				mx = _mm_max_ps(mx, ZRG);
				my = _mm_max_ps(my, ZRG);
				mz = _mm_max_ps(mz, ZRG);
			}
			{
				__m128 CLRG = _mm_set1_ps(clip);
				mx = _mm_min_ps(mx, CLRG);
				my = _mm_min_ps(my, CLRG);
				mz = _mm_min_ps(mz, CLRG);
			}
			{
				__m128 w = _mm_mul_ps(_mm_loadu_ps(&ARR[scidx]), mx);
				w = _mm_mul_ps(w, my);
				w = _mm_mul_ps(w, mz);
				w = _mm_mul_ps(w, _mm_set1_ps(sn));


				_mm_store_ps(&TT[dcidx], _mm_add_ps(_mm_load_ps(&TT[dcidx]), w));
			}

			if (channels == 1) continue;

			c++;

			dcidx = CPIDX(dpidx, c);

			scidx = CPIDX(spidx, c);

			mx = _mm_loadu_ps(&MU_x[scidx]);
			my = _mm_loadu_ps(&MU_y[scidx]);
			mz = _mm_loadu_ps(&MU_z[scidx]);

			mx = _mm_sub_ps(pos_x, mx);
			my = _mm_sub_ps(pos_y, my);
			mz = _mm_sub_ps(pos_z, mz);
			
			mx = abs_ps(mx);
			my = abs_ps(my);
			mz = abs_ps(mz);
			{
				__m128 si = _mm_set1_ps(0.5f + off_table.sigma);
				mx = _mm_sub_ps(si, mx);
				my = _mm_sub_ps(si, my);
				mz = _mm_sub_ps(si, mz);
			}
			{
				__m128 ZRG = _mm_set1_ps(0.0f);
				mx = _mm_max_ps(mx, ZRG);
				my = _mm_max_ps(my, ZRG);
				mz = _mm_max_ps(mz, ZRG);
			}
			{
				__m128 CLRG = _mm_set1_ps(clip);
				mx = _mm_min_ps(mx, CLRG);
				my = _mm_min_ps(my, CLRG);
				mz = _mm_min_ps(mz, CLRG);
			}
			{
				__m128 w = _mm_mul_ps(_mm_loadu_ps(&ARR[scidx]), mx);
				w = _mm_mul_ps(w, my);
				w = _mm_mul_ps(w, mz);
				w = _mm_mul_ps(w, _mm_set1_ps(sn));


				_mm_store_ps(&TT[dcidx], _mm_add_ps(_mm_load_ps(&TT[dcidx]), w));
			}

			if (channels == 2) continue;

			c++;

			dcidx = CPIDX(dpidx, c);

			scidx = CPIDX(spidx, c);
			mx = _mm_loadu_ps(&MU_x[scidx]);
			my = _mm_loadu_ps(&MU_y[scidx]);
			mz = _mm_loadu_ps(&MU_z[scidx]);

			mx = _mm_sub_ps(pos_x, mx);
			my = _mm_sub_ps(pos_y, my);
			mz = _mm_sub_ps(pos_z, mz);
			
			mx = abs_ps(mx);
			my = abs_ps(my);
			mz = abs_ps(mz);
			{
				__m128 si = _mm_set1_ps(0.5f + off_table.sigma);
				mx = _mm_sub_ps(si, mx);
				my = _mm_sub_ps(si, my);
				mz = _mm_sub_ps(si, mz);
			}
			{
				__m128 ZRG = _mm_set1_ps(0.0f);
				mx = _mm_max_ps(mx, ZRG);
				my = _mm_max_ps(my, ZRG);
				mz = _mm_max_ps(mz, ZRG);
			}
			{
				__m128 CLRG = _mm_set1_ps(clip);
				mx = _mm_min_ps(mx, CLRG);
				my = _mm_min_ps(my, CLRG);
				mz = _mm_min_ps(mz, CLRG);
			}
			{
				__m128 w = _mm_mul_ps(_mm_loadu_ps(&ARR[scidx]), mx);
				w = _mm_mul_ps(w, my);
				w = _mm_mul_ps(w, mz);
				w = _mm_mul_ps(w, _mm_set1_ps(sn));


				_mm_store_ps(&TT[dcidx], _mm_add_ps(_mm_load_ps(&TT[dcidx]), w));
			}
			
		}
	}
}


void adv_vel(float * ARR, int tid, int section, int end, float clip, float sn)
{
	__declspec(align(16)) int X0[4];
	__declspec(align(16)) int Y0[4];
	__declspec(align(16)) int Z0[4];
	__declspec(align(16)) int X1[4];
	__declspec(align(16)) int Y1[4];
	__declspec(align(16)) int Z1[4];

	for (int c = 0; c < channels; c++)
	{
		int coffset = COFFSET(c);
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{

			int idx = IDX(x, y, z);
			int pidx = PIDX(x, y, z);
			int cidx = idx + coffset;
			int pcidx = pidx + cpoffset;
			
			
			__m128 X, Y, Z;
			{
				__m128 ts = _mm_set1_ps(timestep);
				__m128 nma = _mm_set1_ps(-off_table.ma);
				__m128 ma = _mm_set1_ps(off_table.ma);
				__m128 b = _mm_set1_ps(GRID_WIDTH);
				X = 
					_mm_sub_ps(
					_mm_sub_ps(_mm_load_ps(&POS_ARRAY_X[idx]), _mm_set1_ps(0.5f)), 
					_mm_min_ps(_mm_max_ps(
					_mm_mul_ps(_mm_load_ps(&V_x[pcidx]), _mm_set1_ps(timestep)), nma), ma));
				X = _mm_add_ps(X, b);
				X = _mm_sub_ps(X, _mm_mul_ps(simd_floor(_mm_div_ps(X, b)), b));

				b = _mm_set1_ps(GRID_HEIGHT);
				Y = 
					_mm_sub_ps(
					_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Y[idx]), _mm_set1_ps(0.5f)), 
					_mm_min_ps(_mm_max_ps(
					_mm_mul_ps(_mm_load_ps(&V_y[pcidx]), _mm_set1_ps(timestep)), nma), ma));
				Y = _mm_add_ps(Y, b);
				Y = _mm_sub_ps(Y, _mm_mul_ps(simd_floor(_mm_div_ps(Y, b)), b));

				b = _mm_set1_ps(GRID_WIDTH);
				Z = 
					_mm_sub_ps(
					_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Z[idx]), _mm_set1_ps(0.5f)), 
					_mm_min_ps(_mm_max_ps(
					_mm_mul_ps(_mm_load_ps(&V_z[pcidx]), _mm_set1_ps(timestep)), nma), ma));
				Z = _mm_add_ps(Z, b);
				Z = _mm_sub_ps(Z, _mm_mul_ps(simd_floor(_mm_div_ps(Z, b)), b));
			}

			__m128 FX, FY, FZ;
			{
				__m128 X0s = simd_floor(X);
				__m128 Y0s = simd_floor(Y);
				__m128 Z0s = simd_floor(Z);

				__m128i I = _mm_cvtps_epi32(X0s);
				__m128i ONERG = _mm_set1_epi32(1);

				_mm_store_si128((__m128i*)&X0[0], I);
				_mm_store_si128((__m128i*)&X1[0], _mm_add_epi32(I, ONERG));

				I = _mm_cvtps_epi32(Y0s);
				_mm_store_si128((__m128i*)&Y0[0], I);
				_mm_store_si128((__m128i*)&Y1[0], _mm_add_epi32(I, ONERG));

				I = _mm_cvtps_epi32(Z0s);
				_mm_store_si128((__m128i*)&Z0[0], I);
				_mm_store_si128((__m128i*)&Z1[0], _mm_add_epi32(I, ONERG));

				FX = _mm_sub_ps(X, X0s);
				FY = _mm_sub_ps(Y, Y0s);
				FZ = _mm_sub_ps(Z, Z0s);

			}


			__m128 NFX = _mm_sub_ps(_mm_set1_ps(1), FX);
			__m128 NFY = _mm_sub_ps(_mm_set1_ps(1), FY);
			__m128 NFZ = _mm_sub_ps(_mm_set1_ps(1), FZ);
			__m128 C = _mm_set_ps(
				ARR[PIDX(X0[3], Y0[3], Z0[3]) + cpoffset],
				ARR[PIDX(X0[2], Y0[2], Z0[2]) + cpoffset],
				ARR[PIDX(X0[1], Y0[1], Z0[1]) + cpoffset],
				ARR[PIDX(X0[0], Y0[0], Z0[0]) + cpoffset]);
			X = _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, NFX), NFY), NFZ);

			C = _mm_set_ps(
				ARR[PIDX(X1[3], Y0[3], Z0[3]) + cpoffset],
				ARR[PIDX(X1[2], Y0[2], Z0[2]) + cpoffset],
				ARR[PIDX(X1[1], Y0[1], Z0[1]) + cpoffset],
				ARR[PIDX(X1[0], Y0[0], Z0[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, FX), NFY), NFZ));

			C = _mm_set_ps(
				ARR[PIDX(X0[3], Y1[3], Z0[3]) + cpoffset],
				ARR[PIDX(X0[2], Y1[2], Z0[2]) + cpoffset],
				ARR[PIDX(X0[1], Y1[1], Z0[1]) + cpoffset],
				ARR[PIDX(X0[0], Y1[0], Z0[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, NFX), FY), NFZ));

			C = _mm_set_ps(
				ARR[PIDX(X1[3], Y1[3], Z0[3]) + cpoffset],
				ARR[PIDX(X1[2], Y1[2], Z0[2]) + cpoffset],
				ARR[PIDX(X1[1], Y1[1], Z0[1]) + cpoffset],
				ARR[PIDX(X1[0], Y1[0], Z0[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, FX), FY), NFZ));

			C = _mm_set_ps(
				ARR[PIDX(X0[3], Y0[3], Z1[3]) + cpoffset],
				ARR[PIDX(X0[2], Y0[2], Z1[2]) + cpoffset],
				ARR[PIDX(X0[1], Y0[1], Z1[1]) + cpoffset],
				ARR[PIDX(X0[0], Y0[0], Z1[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, NFX), NFY), FZ));

			C = _mm_set_ps(
				ARR[PIDX(X1[3], Y0[3], Z1[3]) + cpoffset],
				ARR[PIDX(X1[2], Y0[2], Z1[2]) + cpoffset],
				ARR[PIDX(X1[1], Y0[1], Z1[1]) + cpoffset],
				ARR[PIDX(X1[0], Y0[0], Z1[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, FX), NFY), FZ));

			C = _mm_set_ps(
				ARR[PIDX(X0[3], Y1[3], Z1[3]) + cpoffset],
				ARR[PIDX(X0[2], Y1[2], Z1[2]) + cpoffset],
				ARR[PIDX(X0[1], Y1[1], Z1[1]) + cpoffset],
				ARR[PIDX(X0[0], Y1[0], Z1[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, NFX), FY), FZ));

			C = _mm_set_ps(
				ARR[PIDX(X1[3], Y1[3], Z1[3]) + cpoffset],
				ARR[PIDX(X1[2], Y1[2], Z1[2]) + cpoffset],
				ARR[PIDX(X1[1], Y1[1], Z1[1]) + cpoffset],
				ARR[PIDX(X1[0], Y1[0], Z1[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, FX), FY), FZ));
			

			_mm_store_ps(&TT[pcidx], X);

		}
	}
}

void adv_f(float * ARR, int tid, int section, int end, float clip, float sn)
{
	__declspec(align(16)) int X0[4];
	__declspec(align(16)) int Y0[4];
	__declspec(align(16)) int Z0[4];
	__declspec(align(16)) int X1[4];
	__declspec(align(16)) int Y1[4];
	__declspec(align(16)) int Z1[4];

	for (int c = 0; c < channels; c++)
	{
		int coffset = COFFSET(c);
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{

			int idx = IDX(x, y, z);
			int pidx = PIDX(x, y, z);
			int cidx = idx + coffset;
			int pcidx = pidx + cpoffset;
			
			
			__m128 X, Y, Z;
			{
				__m128 ts = _mm_set1_ps(timestep);
				__m128 nma = _mm_set1_ps(-off_table.ma);
				__m128 ma = _mm_set1_ps(off_table.ma);
				__m128 b = _mm_set1_ps(GRID_WIDTH);
				X = 
					_mm_sub_ps(
					_mm_sub_ps(_mm_load_ps(&POS_ARRAY_X[idx]), _mm_set1_ps(0.5f)), 
					_mm_min_ps(_mm_max_ps(
					_mm_mul_ps(_mm_load_ps(&F_x[pcidx]), _mm_set1_ps(timestep)), nma), ma));
				X = _mm_add_ps(X, b);
				X = _mm_sub_ps(X, _mm_mul_ps(simd_floor(_mm_div_ps(X, b)), b));

				b = _mm_set1_ps(GRID_HEIGHT);
				Y = 
					_mm_sub_ps(
					_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Y[idx]), _mm_set1_ps(0.5f)), 
					_mm_min_ps(_mm_max_ps(
					_mm_mul_ps(_mm_load_ps(&F_y[pcidx]), _mm_set1_ps(timestep)), nma), ma));
				Y = _mm_add_ps(Y, b);
				Y = _mm_sub_ps(Y, _mm_mul_ps(simd_floor(_mm_div_ps(Y, b)), b));

				b = _mm_set1_ps(GRID_WIDTH);
				Z = 
					_mm_sub_ps(
					_mm_sub_ps(_mm_load_ps(&POS_ARRAY_Z[idx]), _mm_set1_ps(0.5f)), 
					_mm_min_ps(_mm_max_ps(
					_mm_mul_ps(_mm_load_ps(&F_z[pcidx]), _mm_set1_ps(timestep)), nma), ma));
				Z = _mm_add_ps(Z, b);
				Z = _mm_sub_ps(Z, _mm_mul_ps(simd_floor(_mm_div_ps(Z, b)), b));
			}

			__m128 FX, FY, FZ;
			{
				__m128 X0s = simd_floor(X);
				__m128 Y0s = simd_floor(Y);
				__m128 Z0s = simd_floor(Z);

				__m128i I = _mm_cvtps_epi32(X0s);
				__m128i ONERG = _mm_set1_epi32(1);

				_mm_store_si128((__m128i*)&X0[0], I);
				_mm_store_si128((__m128i*)&X1[0], _mm_add_epi32(I, ONERG));

				I = _mm_cvtps_epi32(Y0s);
				_mm_store_si128((__m128i*)&Y0[0], I);
				_mm_store_si128((__m128i*)&Y1[0], _mm_add_epi32(I, ONERG));

				I = _mm_cvtps_epi32(Z0s);
				_mm_store_si128((__m128i*)&Z0[0], I);
				_mm_store_si128((__m128i*)&Z1[0], _mm_add_epi32(I, ONERG));

				FX = _mm_sub_ps(X, X0s);
				FY = _mm_sub_ps(Y, Y0s);
				FZ = _mm_sub_ps(Z, Z0s);

			}


			__m128 NFX = _mm_sub_ps(_mm_set1_ps(1), FX);
			__m128 NFY = _mm_sub_ps(_mm_set1_ps(1), FY);
			__m128 NFZ = _mm_sub_ps(_mm_set1_ps(1), FZ);
			__m128 C = _mm_set_ps(
				ARR[PIDX(X0[3], Y0[3], Z0[3]) + cpoffset],
				ARR[PIDX(X0[2], Y0[2], Z0[2]) + cpoffset],
				ARR[PIDX(X0[1], Y0[1], Z0[1]) + cpoffset],
				ARR[PIDX(X0[0], Y0[0], Z0[0]) + cpoffset]);
			X = _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, NFX), NFY), NFZ);

			C = _mm_set_ps(
				ARR[PIDX(X1[3], Y0[3], Z0[3]) + cpoffset],
				ARR[PIDX(X1[2], Y0[2], Z0[2]) + cpoffset],
				ARR[PIDX(X1[1], Y0[1], Z0[1]) + cpoffset],
				ARR[PIDX(X1[0], Y0[0], Z0[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, FX), NFY), NFZ));

			C = _mm_set_ps(
				ARR[PIDX(X0[3], Y1[3], Z0[3]) + cpoffset],
				ARR[PIDX(X0[2], Y1[2], Z0[2]) + cpoffset],
				ARR[PIDX(X0[1], Y1[1], Z0[1]) + cpoffset],
				ARR[PIDX(X0[0], Y1[0], Z0[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, NFX), FY), NFZ));

			C = _mm_set_ps(
				ARR[PIDX(X1[3], Y1[3], Z0[3]) + cpoffset],
				ARR[PIDX(X1[2], Y1[2], Z0[2]) + cpoffset],
				ARR[PIDX(X1[1], Y1[1], Z0[1]) + cpoffset],
				ARR[PIDX(X1[0], Y1[0], Z0[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, FX), FY), NFZ));

			C = _mm_set_ps(
				ARR[PIDX(X0[3], Y0[3], Z1[3]) + cpoffset],
				ARR[PIDX(X0[2], Y0[2], Z1[2]) + cpoffset],
				ARR[PIDX(X0[1], Y0[1], Z1[1]) + cpoffset],
				ARR[PIDX(X0[0], Y0[0], Z1[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, NFX), NFY), FZ));

			C = _mm_set_ps(
				ARR[PIDX(X1[3], Y0[3], Z1[3]) + cpoffset],
				ARR[PIDX(X1[2], Y0[2], Z1[2]) + cpoffset],
				ARR[PIDX(X1[1], Y0[1], Z1[1]) + cpoffset],
				ARR[PIDX(X1[0], Y0[0], Z1[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, FX), NFY), FZ));

			C = _mm_set_ps(
				ARR[PIDX(X0[3], Y1[3], Z1[3]) + cpoffset],
				ARR[PIDX(X0[2], Y1[2], Z1[2]) + cpoffset],
				ARR[PIDX(X0[1], Y1[1], Z1[1]) + cpoffset],
				ARR[PIDX(X0[0], Y1[0], Z1[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, NFX), FY), FZ));

			C = _mm_set_ps(
				ARR[PIDX(X1[3], Y1[3], Z1[3]) + cpoffset],
				ARR[PIDX(X1[2], Y1[2], Z1[2]) + cpoffset],
				ARR[PIDX(X1[1], Y1[1], Z1[1]) + cpoffset],
				ARR[PIDX(X1[0], Y1[0], Z1[0]) + cpoffset]);
			X = _mm_add_ps(X, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(C, FX), FY), FZ));
			

			_mm_store_ps(&TT[pcidx], X);

		}
	}
}