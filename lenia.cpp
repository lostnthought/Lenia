#include "render.h"
#include "lenia.h"
#include <assert.h>
#include <stdio.h>
#include <xmmintrin.h>
#include "reintegration.h"


// THE ARRAYS

// Grid 
__declspec(align(16)) fftwf_complex val[MAX_CHANNELS][GRID_SIZE];

// Grid transposed from fftw complex
__declspec(align(16)) float GRID[GRID_PADDED * MAX_CHANNELS];

// Grid with growth applied early (used for RT options)
__declspec(align(16)) float BACK[GRID_PADDED * MAX_CHANNELS];

// Forward FFT of grid
__declspec(align(16)) fftwf_complex T[MAX_CHANNELS][GRID_SIZE];

// Temp buffer for complex mul
__declspec(align(16)) fftwf_complex TEMP[MAX_KERNELS][GRID_SIZE];

// Weighted Neighbor Sums
__declspec(align(16)) fftwf_complex U[MAX_KERNELS][GRID_SIZE];

// Pos array of offsets indexed by cell IDX()
__declspec(align(16)) float POS_ARRAY_X[GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH];
__declspec(align(16)) float POS_ARRAY_Y[GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH];
__declspec(align(16)) float POS_ARRAY_Z[GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH];

// Moving average of grid
__declspec(align(16)) float GA[GRID_SIZE * MAX_CHANNELS];

// abs(current grid - ga)
__declspec(align(16)) float ALPHA[GRID_PADDED * MAX_CHANNELS];

// alpha used for F computation
__declspec(align(16)) float PREALPHA[GRID_SIZE * MAX_CHANNELS];

// Sum of grid across all channels
__declspec(align(16)) float SUM[GRID_PADDED];

// Sum of grid across all channels blurred (or unblurred) prepped for central difference calc
__declspec(align(16)) float SUM_S[GRID_PADDED];

// Weighted neighbor sums blurred (or unblurred) prepped for central difference calc
__declspec(align(16)) float U_S[GRID_PADDED * MAX_CHANNELS];

// Sum of weighted neibor sums used for mesh gen
__declspec(align(16)) float U_SUM_RB[GRID_PADDED * MAX_CHANNELS];
__declspec(align(16)) float U_SUM[GRID_PADDED * MAX_CHANNELS];

// Back buffer of grid used for mesh gen
__declspec(align(16)) float RBB[GRID_PADDED * MAX_CHANNELS];

// growth and back buffer for mesh gen
__declspec(align(16)) float H_RB[GRID_PADDED * MAX_CHANNELS];
__declspec(align(16)) float H[GRID_PADDED * MAX_CHANNELS];

// Nabla a = gradient of grid sum across all channels
__declspec(align(16)) float NA_x[GRID_SIZE];
__declspec(align(16)) float NA_y[GRID_SIZE];
__declspec(align(16)) float NA_z[GRID_SIZE];

// Nabla u = gradient of growth
__declspec(align(16)) float NU_x[GRID_SIZE * MAX_CHANNELS];
__declspec(align(16)) float NU_y[GRID_SIZE * MAX_CHANNELS];
__declspec(align(16)) float NU_z[GRID_SIZE * MAX_CHANNELS];

// F = Nabla u * (1-alpha) - Nabla a * alpha
__declspec(align(16)) float F_x[GRID_PADDED * MAX_CHANNELS];
__declspec(align(16)) float F_y[GRID_PADDED * MAX_CHANNELS];
__declspec(align(16)) float F_z[GRID_PADDED * MAX_CHANNELS];

// Velocity = V + F * timestep - gamma * V 
__declspec(align(16)) float V_x[GRID_PADDED * MAX_CHANNELS];
__declspec(align(16)) float V_y[GRID_PADDED * MAX_CHANNELS];
__declspec(align(16)) float V_z[GRID_PADDED * MAX_CHANNELS];

// Displacement for each cell calculated by V or F * timestep and clamped by ma = sigma - dd
__declspec(align(16)) float MU_x[GRID_PADDED * MAX_CHANNELS];
__declspec(align(16)) float MU_y[GRID_PADDED * MAX_CHANNELS];
__declspec(align(16)) float MU_z[GRID_PADDED * MAX_CHANNELS];

// Accumulated reintegration array 
__declspec(align(16)) float TT[GRID_PADDED * MAX_CHANNELS];

// Growth lookup table
float BELL_LUT[BELL_LUT_SIZE];
void init_bell_lut()
{
	for (int u = 0; u < 256; u++)
	for (int m = 0; m < 256; m++)
	for (int s = 0; s < 256; s++)
	{
		BELL_LUT[m * 256 * 256 + s * 256 + u] = bell((float)u / 255, (float)m / 255, (float)s / 255);	
	}
}

// Old optional kernel for some gradient blurring (not in use)
__declspec(align(16)) float gaussian_kernel[3][3][3];
void init_gaussian_kernel()
{
	const float s = 1.0f;
	float sum = 0.0f;
	for (float dx = -1; dx <= 1; dx++)
	for (float dy = -1; dy <= 1; dy++)
	for (float dz = -1; dz <= 1; dz++)
	{
		float d2 = dx * dx + dy * dy + dz * dz;
		float w = expf(-d2 / (s * s * s));
		gaussian_kernel[(int)(dx + 1)][(int)(dy + 1)][(int)(dz + 1)] = w;
		sum += w;
	}
	for (float dx = -1; dx <= 1; dx++)
	for (float dy = -1; dy <= 1; dy++)
	for (float dz = -1; dz <= 1; dz++)
	{
		gaussian_kernel[(int)(dx + 1)][(int)(dy + 1)][(int)(dz + 1)] /= sum;
	}
}

void init_pos_array()
{
	for (int x = 0; x < GRID_WIDTH; x++)
	for (int y = 0; y < GRID_HEIGHT; y++)
	for (int z = 0; z < GRID_WIDTH; z++)
	{
		POS_ARRAY_X[z * GRID_WIDTH * GRID_HEIGHT + y * GRID_WIDTH + x] = x + 0.5f;
		POS_ARRAY_Y[z * GRID_WIDTH * GRID_HEIGHT + y * GRID_WIDTH + x] = y + 0.5f;
		POS_ARRAY_Z[z * GRID_WIDTH * GRID_HEIGHT + y * GRID_WIDTH + x] = z + 0.5f;
	}
}

float SIGMA[3][3][3];
void init_offset_table(Button * b)
{
	int s = 1;
	off_table.row_size = (2 * off_table.dd + 1) * (2 * off_table.dd + 1) * (2 * off_table.dd + 1);
	off_table.ma = off_table.dd - off_table.sigma;

	if (off_table.X) free(off_table.X);
	off_table.X = NULL;
	if (off_table.Y) free(off_table.Y);
	off_table.Y = NULL;
	if (off_table.Z) free(off_table.Z);
	off_table.Z = NULL;

	off_table.X = (int *)calloc(off_table.row_size, sizeof(int));
	off_table.Y = (int *)calloc(off_table.row_size, sizeof(int));
	off_table.Z = (int *)calloc(off_table.row_size, sizeof(int));


	int index = 0;
	for (int dx = -off_table.dd; dx <= off_table.dd; dx++)
	for (int dy = -off_table.dd; dy <= off_table.dd; dy++)
	for (int dz = -off_table.dd; dz <= off_table.dd; dz++)
	{
		off_table.X[index] = dx;
		off_table.Y[index] = dy;
		off_table.Z[index] = dz;
		
		index++;
	}
	float sum = 0.0f;

	for (float dx = -1; dx <= 1; dx++)
	for (float dy = -1; dy <= 1; dy++)
	for (float dz = -1; dz <= 1; dz++)
	{
		float d2 = dx * dx + dy * dy + dz * dz;
		float w = expf(-d2 / (2 * off_table.sigma * off_table.sigma));
		SIGMA[(int)(dx + 1)][(int)(dy + 1)][(int)(dz + 1)] = w;
		sum += w;
		
	}
	for (float dx = -1; dx <= 1; dx++)
	for (float dy = -1; dy <= 1; dy++)
	for (float dz = -1; dz <= 1; dz++)
	{
		SIGMA[(int)(dx + 1)][(int)(dy + 1)][(int)(dz + 1)] /= sum;
	}

}

void init_fftw_plans()
{


	reverse_fft = fftwf_plan_dft_3d(GRID_WIDTH, GRID_HEIGHT, GRID_WIDTH, NULL, NULL, FFTW_BACKWARD, FFTW_ESTIMATE);
	forward_fft = fftwf_plan_dft_3d(GRID_WIDTH, GRID_HEIGHT, GRID_WIDTH, NULL, NULL, FFTW_FORWARD, FFTW_ESTIMATE);

}

// init a new kernel and forward fft it
void init_kernel(Kernel * kernel)
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
	

	for (int r_x = -offset_w; r_x < offset_w; r_x++)
	for (int r_y = -offset_h; r_y < offset_h; r_y++)
	for (int r_z = -offset_w; r_z < offset_w; r_z++)
	{
		float x = fabsf((float)r_x);
		float y = fabsf((float)r_y);
		float z = fabsf((float)r_z);
		
		float dist = sqrt(x * x + y * y + z * z);

		//if (dist < r)
		//{
			 dist = dist / r * beta_size / kernel->r_a;
			
			if (dist < beta_size)
			{
				int real_x = emod(r_x, GRID_WIDTH);
				int real_y = emod(r_y, GRID_HEIGHT);
				int real_z = emod(r_z, GRID_WIDTH);


				int index = (int)MIN((int)floorf(dist), (int)beta_size - 1); 

				float final = kernel->b[index] * (float)bell(fmodf((float)dist, 1.0f), 0.5f, 0.15f);

				kernel->grid[(real_x) + (real_y) * GRID_WIDTH + (real_z) * GRID_WIDTH * GRID_HEIGHT][0] = final;
				kernel->real_k[r_x + offset_w + (r_y + offset_h) * GRID_WIDTH + (r_z + offset_w) * GRID_WIDTH * GRID_HEIGHT] = final; 

				k_cells += final;

			}

		//}

	}
	for (int x = 0; x < GRID_SIZE; x++)
	{

		kernel->grid[x][0] /= k_cells;
		kernel->real_k[x] /= k_cells;
		kernel->grid[x][1] = 0;

	}

	fftwf_execute_dft(forward_fft, kernel->grid, kernel->kg);

}

// functions used for changing the field during runtime

void clear_buffers()
{
	for (int z = 0; z < GRID_WIDTH; z++)
	for (int y = 0; y < GRID_HEIGHT; y++)
	for (int x = 0; x < GRID_WIDTH; x++)
	{
		for (int c = 0; c < MAX_CHANNELS; c++)
		{
			int cidx = CIDX(IDX(x, y, z), c);
			int pcidx = CPIDX(PIDX(x, y, z), c);
			GA[cidx] = 0;
			ALPHA[pcidx] = 0;
			V_x[pcidx] = 0;
			V_y[pcidx] = 0;
			V_z[pcidx] = 0;
		}
	}
}

// randomize field
void init_conway(){

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = 0; z < GRID_WIDTH; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			float a = (float)(rand() % 1000) / 1000;
			
			val[c][IDX(x, y, z)][0] = a; 
			val[c][IDX(x, y, z)][1] = 0; 
			GRID[PIDX(x, y, z) + cpoffset] = a;

		}
	}
}

// seeds grid with blobs based on a neighborhood around kernel size
void seed_grid_gaussian(int r)
{
	for (int i = 0; i < 5; i++)
	{
		int cx = rand() % GRID_WIDTH;
		int cy = rand() % GRID_HEIGHT;
		int cz = rand() % GRID_WIDTH;

		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int z = 0; z < GRID_WIDTH; z++)
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int pidx = PIDX(x, y, z);
				int dx = x - cx;
				int dy = y - cy;
				int dz = z - cz;
				float dist = dx * dx + dy * dy + dz * dz;
				float v = exp(-dist / (2 * r * r));
				//float val = sqrt(dist);
				int idx = x + y * GRID_WIDTH + z * GRID_WIDTH * GRID_HEIGHT;
				val[c][idx][0] += v;
				GRID[pidx + cpoffset] += v;
				val[c][idx][1] = 0;
			}
		}
	}
}

void generate_neighborhood_random(fftwf_complex * cll, float * ksum)
{

	int real_r = 0;
	int offset_w = GRID_WIDTH / 2;
	int offset_h = GRID_HEIGHT / 2;
	for (int row = -offset_w; row < offset_w; row++)
	{
		int real_c = 0;
		for (int col = -offset_h; col < offset_h; col++)
		{
			int real_d = 0;
			for (int dep = -offset_w; dep < offset_w; dep++)
			{
				float x = (float)abs(row);
				float y = (float)abs(col);
				float z = (float)abs(dep);
				float r = sqrt(x * x + y * y + z * z);
				if (r < neighborhood)
				{
					int offset_r = emod(row, GRID_WIDTH);
					int offset_c = emod(col, GRID_HEIGHT);
					int offset_d = emod(dep, GRID_WIDTH);
					int random_offset = rand() % 20 - 10;
					cll[real_r + real_c * (offset_w * 2) + real_d * (offset_w * 2) * (offset_h * 2)][0] =
						MIN(ksum[(int)MAX(
						offset_d * GRID_WIDTH * GRID_HEIGHT + offset_c * GRID_WIDTH + offset_r, 0)] / kernel_count, 1);
				}
			}
		}
	}
}


void generate_random()
{

	neighborhood = 1 + (rand() % 8);
	kernel_count = rand() % 10 + 1;

	for (int i = 0; i < MAX_CHANNELS; i++)
	{
		ch[i].asymptotic = 0;
		ch[i].soft_clip = 0;
	}

	float ksum[GRID_SIZE];
	for (int i = 0; i < kernel_count; i++)
	{

		Kernel * k = &kernels[i];

		k->beta_size = rand() % 2 + 1;

		k->radius = neighborhood /*+ rand() % (int)(2) - 1*/;

		for (int e = 0; e < k->beta_size; e++)
		{
			k->b[e] = (float)(rand() % 100)/100;
		}
		k->m = (float)(100 + rand() % 500) / 1000;

		k->s = (k->m / 10) + (float)(rand() % 200) / 1000;
	
		if (channels == 1 || i == 0)
		{
			k->c0 = 0; 
			k->c1 = 0;
		} else 
		{
			k->c0 = rand() % channels;
			k->c1 = rand() % channels;
		}
		k->h = 1;
		k->r_a = 1;
		init_kernel(k);
	}
	
	Vertex * v;
	kernel_verts->Lock(0, 0, (void **) &v, D3DLOCK_DISCARD);
	memset(v, 0, GRID_SIZE * 6 * MAX_KERNELS * sizeof(Vertex));

	rebuild_kernel_mesh(v);

	kernel_verts->Unlock();
	memset(val, 0, sizeof(fftwf_complex) * GRID_SIZE * MAX_CHANNELS);
	memset(GRID, 0, sizeof(float) * GRID_PADDED * MAX_CHANNELS);

	clear_buffers();

	seed_grid_gaussian(neighborhood);
	redraw_config = true;

}


// sets up arrays at the beginning of a new sim step 
// complex multiplies
// reverse ffts to get U
void fft_step(TD * t){
	
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	for (int c = 0; c < channels; c++)
	{
		int coffset = COFFSET(c);
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = IDX(x, y, z);
			int cidx = coffset + idx;
			
			int pidx = PIDX(x, y, z);
			int pcidx = pidx + cpoffset;

			U_SUM_RB[pcidx] = U_SUM[pcidx];
			U_SUM[pcidx] = 0;

			H_RB[pcidx] = H[pcidx];
			H[pcidx] = 0;
			
			NA_x[idx] = 0;
			NA_y[idx] = 0;
			NA_z[idx] = 0;

			NU_x[cidx] = 0;
			NU_y[cidx] = 0;
			NU_z[cidx] = 0;
		}
	}	

	int ksection = (kernel_count / SIM_THREADS) * t->id;
	int kend = (t->id == SIM_THREADS - 1 ? kernel_count : ksection + kernel_count / SIM_THREADS);

	for (int k = ksection; k < kend; k++)
	{
		Kernel * kernel = &kernels[k];

		int c = kernel->c0;

		for (int i = 0; i < GRID_SIZE; i++)
		{

			TEMP[k][i][0] = 
				T[c][i][0] * kernel->kg[i][0] - T[c][i][1] * kernel->kg[i][1];
			TEMP[k][i][1] = 
				T[c][i][0] * kernel->kg[i][1] + T[c][i][1] * kernel->kg[i][0];

		}

		fftwf_execute_dft(reverse_fft, TEMP[k], U[k]);
	}

}

// normalizes U
// calculates growth
// applies it to back optionally
void mt_fft_step(TD * t)
{

	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	for (int k = 0; k < kernel_count; k++)
	{
		Kernel * kernel = &kernels[k];
		float m = kernel->m; float s = kernel->s; float h = kernel->h;
		int c = kernels[k].c1;
		if (ch[c].asymptotic)
		{
			for (int z = section; z < end; z++)
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{

				int i = IDX(x, y, z);
				int cidx = CIDX(i, c);
				int cpidx = CPIDX(PIDX(x, y, z), c);

				float u = U[k][i][0] / GRID_SIZE;

				U_SUM[cpidx] += u;

				H[cpidx] += (bell_lookup(u, m, s) - GRID[cpidx]) * h;
			}
		}
		else
		{
			for (int z = section; z < end; z++)
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int i = IDX(x, y, z);
				int cpidx = CPIDX(PIDX(x, y, z), c);
				float u = U[k][i][0] / GRID_SIZE;
				U_SUM[cpidx] += u;
				H[cpidx] += bell_lookup(u, m, s) * h;
				
			}
		}
	}


	
	if (GROWTH_TYPE != GROWTH_NOGROWTH && GROWTH_TYPE != GROWTH_TRANSPOSITION)
	{
		int csection = (channels / SIM_THREADS) * t->id;
		int cend = (t->id == SIM_THREADS - 1 ? channels : csection + channels / SIM_THREADS);
		__m128 t = _mm_set1_ps(timestep);
		for (int c = 0; c < channels; c++)
		{
			if (ch[c].soft_clip)
			{
				static const __m128 ZERORG = _mm_set1_ps(0);
				static const __m128 ONERG = _mm_set1_ps(1);
				static const __m128 ZEROFIVERG = _mm_set1_ps(0.5f);
				static const __m128 NEGFOURRG = _mm_set1_ps(-4.0f);
				int cpoffset = CPOFFSET(c);

				for (int z = section; z < end; z++)
				for (int y = 0; y < GRID_HEIGHT; y++)
				for (int x = 0; x < GRID_WIDTH; x+=4)
				{
					int cpidx = PIDX(x, y, z) + cpoffset;

					__m128 g = _mm_load_ps(&GRID[cpidx]);
					__m128 h = _mm_load_ps(&H[cpidx]);
					h = _mm_add_ps(g, _mm_mul_ps(h, t));
					h = _mm_mul_ps(NEGFOURRG, _mm_sub_ps(h, ZEROFIVERG));
					h = fastexp(h);
					h = _mm_add_ps(ONERG, h);
					h = _mm_div_ps(ONERG, h);
					_mm_store_ps(&BACK[cpidx], _mm_min_ps(_mm_max_ps(h, ZERORG), ONERG));			
				}
			}
			else
			{
				
				static const __m128 ZERORG = _mm_set1_ps(0);
				static const __m128 ONERG = _mm_set1_ps(1);
				int cpoffset = CPOFFSET(c);

				for (int z = section; z < end; z++)
				for (int y = 0; y < GRID_HEIGHT; y++)
				for (int x = 0; x < GRID_WIDTH; x+=4)
				{
					int cpidx = PIDX(x, y, z) + cpoffset;

					__m128 g = _mm_load_ps(&GRID[cpidx]);
					__m128 h = _mm_load_ps(&H[cpidx]);
					_mm_store_ps(&BACK[cpidx], _mm_min_ps(_mm_max_ps(_mm_add_ps(g, _mm_mul_ps(h, t)), ZERORG), ONERG));		
				}
			}
		}
	}
	
}

// calculates the gradient for sum
void sobel_convolve(TD * t)
{

	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);

	if (GRADIENT_MODE == GRADIENT_SOBEL)
	{
		sum_sobel(section, end);
	}
	else if (GRADIENT_MODE == GRADIENT_CENTRAL_DIFF)
	{
		sum_central_diff(section, end);
	}
	t->is_done = true;
}

// calculates gradient for growth
// precalculates alpha so we can vectorize f calculation
void sobel_convolve_h(TD * t)
{
	
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);

	if (GRADIENT_MODE == GRADIENT_SOBEL)
	{
		h_sobel(section, end);
	}
	else if (GRADIENT_MODE == GRADIENT_CENTRAL_DIFF)
	{
		h_central_diff(section, end);
	}

	float * gr;
	switch(ALPHA_MODE)
	{
		case ALPHA_GRID:
			gr = GRID;
			break;
		case ALPHA_GRID_GROWTH:
			gr = BACK;
			break;
		case ALPHA_GA:
			gr = GA;
			break;
		case ALPHA_GA_DIFF:
			gr = ALPHA;
			break;
	}
	__m128 ZERORG = _mm_set1_ps(0);
	__m128 ONERG = _mm_set1_ps(1);
	__m128 TA = _mm_set1_ps(THETA_A);
	for (int c = 0; c < channels; c++)
	{
		int coffset = CPOFFSET(c);

		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{
			int pidx = PIDX(x, y, z);
			int pcidx = CPIDX(pidx, c);
			int cidx = coffset + IDX(x, y, z);
			float sq = fabsf(GRID[pcidx] - GA[cidx]);
			__m128 gr = _mm_load_ps(&GRID[pcidx]);
			__m128 ga = _mm_load_ps(&GA[cidx]);
			_mm_store_ps(&ALPHA[pcidx], abs_ps(_mm_sub_ps(gr, ga)));

			__m128 pg = _mm_div_ps(gr, TA);
			// this disregards n in favor of vectorization. technically possible to respect n as a float
			// but impractical in terms of performance
			pg = _mm_mul_ps(pg, pg);
			_mm_store_ps(&PREALPHA[cidx], _mm_min_ps(_mm_max_ps(pg, ZERORG), ONERG));
		}
	}

	
	t->is_done = true;
}

// computes f what else would you expect
void compute_f(TD * t)
{

	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	
	for (int c = 0; c < channels; c++)
	{
		int coffset = CPOFFSET(c);

		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{
			int index = IDX(x, y, z);
			int pidx = PIDX(x, y, z);
			
			int cidx = CIDX(index, c);

			int pcidx = pidx + coffset;

			__m128 nu_x = _mm_load_ps(&NU_x[cidx]);
			__m128 nu_y = _mm_load_ps(&NU_y[cidx]);
			__m128 nu_z = _mm_load_ps(&NU_z[cidx]);

			__m128 na_x = _mm_load_ps(&NA_x[index]);
			__m128 na_y = _mm_load_ps(&NA_y[index]);
			__m128 na_z = _mm_load_ps(&NA_z[index]);

			__m128 alpha = _mm_load_ps(&PREALPHA[cidx]);
			__m128 alpha1 = _mm_sub_ps(_mm_set1_ps(1.0f), alpha);
			_mm_store_ps(&F_x[pcidx], _mm_sub_ps(_mm_mul_ps(nu_x, alpha1), _mm_mul_ps(na_x, alpha)));
			_mm_store_ps(&F_y[pcidx], _mm_sub_ps(_mm_mul_ps(nu_y, alpha1), _mm_mul_ps(na_y, alpha)));
			_mm_store_ps(&F_z[pcidx], _mm_sub_ps(_mm_mul_ps(nu_z, alpha1), _mm_mul_ps(na_z, alpha)));


		}

	}

}


// for mu computation see reintegration.c



// reintegrates the mass or applies advection. this is the transport step
void reintegration_step(TD * t)
{
	// interpolates between 2d or 3d normalization
	float sigma = off_table.sigma;
	float s2 = 1 / (4 * sigma * sigma);
	float s23 = 1 / (2 * sigma * sigma * sigma);
	float sn = s2 * (1 - S_NORM) + (S_NORM * s23);

	// maximum mass per x, y, z as seen in the flow lenia formula
	float clip = MIN(2.0f * off_table.sigma, 1.0f);

	int tid = t->id;
	int section = tid * (GRID_WIDTH / SIM_THREADS);
	int end = (tid == SIM_THREADS - 1 ? GRID_WIDTH : section + (GRID_WIDTH / SIM_THREADS));

	switch(MODE_TYPE)
	{
		// singular reintegration mode (no transport (very fast), amplify noise)
		case MODE_SRT:
			if (GROWTH_TYPE == GROWTH_NOGROWTH || GROWTH_TYPE == GROWTH_TRANSPOSITION){
				srt(GRID, tid, section, end, clip, sn);
			} else {
				srt(BACK, tid, section, end, clip, sn);
			}

			break;
		// main mode (transport the mass within dd radius of destination cells)
		case MODE_RT:
			{
				if (BORDER == BORDER_TORUS)
				{
					if (GROWTH_TYPE == GROWTH_NOGROWTH || GROWTH_TYPE == GROWTH_TRANSPOSITION)
					{
						rt_border_torus(GRID, tid, section, end, clip, sn);
					}
					else
					{
						rt_border_torus(BACK, tid, section, end, clip, sn);
					}

				}
				else
				{
					if (GROWTH_TYPE == GROWTH_NOGROWTH || GROWTH_TYPE == GROWTH_TRANSPOSITION)
					{
						rt_border_other(GRID, tid, section, end, clip, sn);
					}
					else
					{
						rt_border_other(BACK, tid, section, end, clip, sn);
					}

				}

			}
			break;
		// lagranian advection
		case MODE_ADV:
			{
				float * ARR = GROWTH_TYPE == GROWTH_NOGROWTH || GROWTH_TYPE == GROWTH_TRANSPOSITION ? GRID : BACK;

				if (ADV_VEL)
				{
					adv_vel(ARR, tid, section, end, clip, sn);
				}
				else
				{
					adv_f(ARR, tid, section, end, clip, sn);
				}
			
			}
			break;
	}

	// optional osciallators
	// the only part of the codebase that isn't optimized, very optional toys that are super fast anyway
	if (OSC_LA || OSC_FLOW || OSC_DIV_FLOW)
	{

		for (int c = 0; c < channels; c++)
		{
			int coffset = COFFSET(c);
			int cpoffset = CPOFFSET(c);

			for (int z = section; z < end; z++)
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = IDX(x, y, z);
				int pidx = PIDX(x, y, z);

				int cidx = idx + coffset;
				int cpidx = pidx + cpoffset;

				int mx1 = PIDX((x + 1), y, z);
				int mx0 = PIDX((x - 1), + y, z);
				int my1 = PIDX(x, (y + 1), z);
				int my0 = PIDX(x, (y - 1), z);
				int mz1 = PIDX(x, y, (z + 1));
				int mz0 = PIDX(x, y, (z - 1));
				
				float fx1 = GRID[CPIDX(mx1, c)];
				float fx0 = GRID[CPIDX(mx0, c)];

				float fy1 = GRID[CPIDX(my1, c)];
				float fy0 = GRID[CPIDX(my0, c)];

				float fz1 = GRID[CPIDX(mz1, c)];
				float fz0 = GRID[CPIDX(mz0, c)];

				float diff_x = fx1 - fx0;
				float diff_y = fy1 - fy0;
				float diff_z = fz1 - fz0;

				if (OSC_LA)
				{
					float la = fx1 + fx0 + fy1 + fy0 + fz1 + fz0 - 6 * GRID[cpidx];
					
					TT[cpidx] += la * OSC_LA_AMT;
					
				}
				if (OSC_FLOW)
				{
					if (OSC_VEL)
					{
						float flow = 
							V_x[cpidx] * diff_x +
							V_y[cpidx] * diff_y +
							V_z[cpidx] * diff_z;
						TT[cpidx] += flow * OSC_FLOW_AMT;
					}
					else
					{
						float flow = 
							F_x[cpidx] * diff_x +
							F_y[cpidx] * diff_y +
							F_z[cpidx] * diff_z;
						TT[cpidx] += flow * OSC_FLOW_AMT;
					}
					
				}
				if (OSC_DIV_FLOW)
				{
					float jx1 = 0, jy1 = 0, jz1 = 0;
					int mx0p = CPIDX(mx0, c);
					int mx1p = CPIDX(mx1, c);
					int my0p = CPIDX(my0, c);
					int my1p = CPIDX(my1, c);
					int mz0p = CPIDX(mz0, c);
					int mz1p = CPIDX(mz1, c);

					if (OSC_VEL)
					{
						jx1 = 
							V_x[mx1p] * ALPHA[mx1p] -
							V_x[mx0p] * ALPHA[mx0p];
						jy1 = 
							V_y[my1p] * ALPHA[my1p] -
							V_y[my0p] * ALPHA[my0p];
						jz1 = 
							V_z[mz1p] * ALPHA[mz1p] -
							V_z[mz0p] * ALPHA[mz0p];
					}
					else 
					{
						jx1 = 
							F_x[mx1p] * ALPHA[mx1p] -
							F_x[mx0p] * ALPHA[mx0p];
						jy1 = 
							F_y[my1p] * ALPHA[my1p] -
							F_y[my0p] * ALPHA[my0p];
						jz1 = 
							F_z[mz1p] * ALPHA[mz1p] -
							F_z[mz0p] * ALPHA[mz0p];

					}
					float divj = jx1 + jy1 + jz1;

					TT[cpidx] += divj * OSC_DIV_FLOW_AMT;
				}

			}
		}

	}
	if (tid == 0)
	{
		memset(SUM, 0, sizeof(SUM));
	}
	
	t->is_done = true;
}


// the cleanup for the sim step. transposes to fftw complex so we can fft again, also applies growth if tr
// growth is set
// computes moving averages for grid and sum
void transpose_to_row_major(TD * t)
{

	int section = t->id * (GRID_WIDTH / SIM_THREADS);
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + (GRID_WIDTH / SIM_THREADS));
	
	if (GROWTH_TYPE == GROWTH_TRANSPOSITION)
	{
		__m128 t = _mm_set1_ps(timestep);
		for (int c = 0; c < channels; c++)
		{
			if (ch[c].soft_clip)
			{
				static const __m128 ZERORG = _mm_set1_ps(0);
				static const __m128 ONERG = _mm_set1_ps(1);
				static const __m128 ZEROFIVERG = _mm_set1_ps(0.5f);
				static const __m128 NEGFOURRG = _mm_set1_ps(-4.0f);
				int cpoffset = CPOFFSET(c);

				for (int z = section; z < end; z++)
				for (int y = 0; y < GRID_HEIGHT; y++)
				for (int x = 0; x < GRID_WIDTH; x+=4)
				{
					int cpidx = PIDX(x, y, z) + cpoffset;

					__m128 g = _mm_load_ps(&GRID[cpidx]);
					__m128 h = _mm_load_ps(&H[cpidx]);
					h = _mm_add_ps(g, _mm_mul_ps(h, t));
					h = _mm_mul_ps(NEGFOURRG, _mm_sub_ps(h, ZEROFIVERG));
					h = fastexp(h);
					h = _mm_add_ps(ONERG, h);
					h = _mm_div_ps(ONERG, h);
					_mm_store_ps(&BACK[cpidx], _mm_min_ps(_mm_max_ps(h, ZERORG), ONERG));			
				}
			}
			else
			{
				static const __m128 ZERORG = _mm_set1_ps(0);
				static const __m128 ONERG = _mm_set1_ps(1);
				int cpoffset = CPOFFSET(c);

				for (int z = section; z < end; z++)
				for (int y = 0; y < GRID_HEIGHT; y++)
				for (int x = 0; x < GRID_WIDTH; x+=4)
				{
					int cpidx = PIDX(x, y, z) + cpoffset;

					__m128 g = _mm_load_ps(&TT[cpidx]);
					__m128 h = _mm_load_ps(&H[cpidx]);
					_mm_store_ps(&TT[cpidx], _mm_min_ps(_mm_max_ps(_mm_add_ps(g, _mm_mul_ps(h, t)), ZERORG), ONERG));		
				}
			}
		}
	}
	else
	{
		static const __m128 ZERORG = _mm_set1_ps(0.0f);
		static const __m128 ONERG = _mm_set1_ps(1.0f);
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int z = section; z < end; z++)
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x+=4)
			{
				int pidx = PIDX(x, y, z) + cpoffset;
				_mm_store_ps(&TT[pidx], _mm_min_ps(_mm_max_ps(_mm_load_ps(&TT[pidx]), ZERORG), ONERG));
			}
		}
	}

	{
		__m128 STAM = _mm_set1_ps(1.0f - STA);
		__m128 STAS = _mm_set1_ps(STA);
		__m128 GTAS = _mm_set1_ps(GTA);
		__m128 GTAM = _mm_set1_ps(1.0f - GTA);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		for (int x = 0; x < GRID_WIDTH; x+=4)
		{
			int i = IDX(x, y, z);
			int pidx = PIDX(x, y, z);
			__m128 ttsum = _mm_set1_ps(0.0f);

			for (int c = 0; c < channels; c++)
			{
			
				int cidx = CIDX(i, c);
				int pcidx = CPIDX(pidx, c);

				__m128 tt = _mm_load_ps(&TT[pcidx]);
				ttsum = _mm_add_ps(ttsum, tt);

				_mm_store_ps(&TT[pcidx], _mm_set1_ps(0.0f));

				_mm_store_ps(&GA[cidx],
					_mm_add_ps(
					_mm_mul_ps(GTAM, _mm_load_ps(&GA[cidx])), 
					_mm_mul_ps(tt, GTAS)));

				_mm_store_ps(&GRID[pcidx], tt);

				_mm_store_ps(&val[c][i][0], _mm_unpacklo_ps(tt, _mm_set1_ps(0.0f)));
				_mm_store_ps(&val[c][i+2][0], _mm_unpackhi_ps(tt, _mm_set1_ps(0.0f)));

			}

			_mm_store_ps(&SUM[pidx],
					_mm_add_ps(
					_mm_mul_ps(STAM, _mm_load_ps(&SUM[pidx])), 
					_mm_mul_ps(ttsum, STAS)));
			
		}
	}

	int csection = (channels / SIM_THREADS) * t->id;
	int cend = (t->id == SIM_THREADS - 1 ? channels : csection + channels / SIM_THREADS);
	for (int c = csection; c < cend; c++)
	{	
		fftwf_execute_dft(forward_fft, val[c], T[c]);
	}

}


