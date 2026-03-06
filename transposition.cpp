#include "render.h"
#include "transposition.h"


// torus
void transpose_central_diff_borders(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(-1, y, z);
			int edx = PIDX(GRID_WIDTH - 1, y, z);

			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = SUM_S[eidx];
			U_S[cidx] = U_S[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(GRID_WIDTH, y, z);
			int edx = PIDX(0, y, z);
		

			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = SUM_S[eidx];
			U_S[cidx] = U_S[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, -1, z);
			int edx = PIDX(x, GRID_HEIGHT - 1, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = SUM_S[eidx];
			U_S[cidx] = U_S[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section ; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, GRID_HEIGHT, z);
			int edx = PIDX(x, 0, z);
			
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = SUM_S[eidx];
			U_S[cidx] = U_S[eidx];
		}
	}

	if (t->id == 0 || SIM_THREADS == 1)
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, -1);
				int edx = PIDX(x, y, GRID_WIDTH-1);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM_S[cidx] = SUM_S[eidx];
				U_S[cidx] = U_S[eidx];
			}
			
		}
	}

	if (t->id == 1 || SIM_THREADS == 1) 
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, GRID_WIDTH);
				int edx = PIDX(x, y, 0);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM_S[cidx] = SUM_S[eidx];
				U_S[cidx] = U_S[eidx];
			}
		}
	}
	
}

void wall_central_diff_borders(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(-1, y, z);
			int edx = PIDX(0, y, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = SUM_S[eidx];
			U_S[cidx] = U_S[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(GRID_WIDTH, y, z);
			int edx = PIDX(GRID_WIDTH-1, y, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = SUM_S[eidx];
			U_S[cidx] = U_S[eidx];
		}
	}
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, -1, z);
			int edx = PIDX(x, 0, z);
			
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = SUM_S[eidx];
			U_S[cidx] = U_S[eidx];
		}
	}
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section ; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, GRID_HEIGHT, z);
			int edx = PIDX(x, GRID_HEIGHT-1, z);

			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = SUM_S[eidx];
			U_S[cidx] = U_S[eidx];
		}
	}

	if (t->id == 0 || SIM_THREADS == 1)
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, -1);
				int edx = PIDX(x, y, 0);

				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM_S[cidx] = SUM_S[eidx];
				U_S[cidx] = U_S[eidx];
			}
			
		}
	}

	if (t->id == 1 || SIM_THREADS == 1) 
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, GRID_WIDTH);
				int edx = PIDX(x, y, GRID_WIDTH-1);

				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM_S[cidx] = SUM_S[eidx];
				U_S[cidx] = U_S[eidx];
			}
		}
	}
	
}
void zero_central_diff_borders(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(-1, y, z);
			int edx = PIDX(GRID_WIDTH - 1, y, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = 0;
			U_S[cidx] = 0;
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(GRID_WIDTH, y, z);
			int edx = PIDX(0, y, z);
			
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = 0;
			U_S[cidx] = 0;
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, -1, z);
			int edx = PIDX(x, GRID_HEIGHT - 1, z);

			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = 0;
			U_S[cidx] = 0;
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section ; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, GRID_HEIGHT, z);
			int edx = PIDX(x, 0, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM_S[cidx] = 0;
			U_S[cidx] = 0;
		}
	}

	if (t->id == 0 || SIM_THREADS == 1)
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, -1);
				int edx = PIDX(x, y, GRID_WIDTH-1);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM_S[cidx] = 0;
				U_S[cidx] = 0;
			}
			
		}
	}

	if (t->id == 1 || SIM_THREADS == 1) 
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, GRID_WIDTH);
				int edx = PIDX(x, y, 0);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM_S[cidx] = 0;
				U_S[cidx] = 0;
			}
		}
	}
	
}





// torus
void transpose_sobel_borders(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(-1, y, z);
			int edx = PIDX(GRID_WIDTH - 1, y, z);

			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = SUM[eidx];
			H[cidx] = H[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(GRID_WIDTH, y, z);
			int edx = PIDX(0, y, z);
		

			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = SUM[eidx];
			H[cidx] = H[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, -1, z);
			int edx = PIDX(x, GRID_HEIGHT - 1, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = SUM[eidx];
			H[cidx] = H[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section ; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, GRID_HEIGHT, z);
			int edx = PIDX(x, 0, z);
			
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = SUM[eidx];
			H[cidx] = H[eidx];
		}
	}

	if (t->id == 0 || SIM_THREADS == 1)
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, -1);
				int edx = PIDX(x, y, GRID_WIDTH-1);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM[cidx] = SUM[eidx];
				H[cidx] = H[eidx];
			}
			
		}
	}

	if (t->id == 1 || SIM_THREADS == 1) 
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, GRID_WIDTH);
				int edx = PIDX(x, y, 0);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM[cidx] = SUM[eidx];
				H[cidx] = H[eidx];
			}
		}
	}
	
}



// wall
void wall_sobel_borders(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(-1, y, z);
			int edx = PIDX(0, y, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = SUM[eidx];
			H[cidx] = H[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(GRID_WIDTH, y, z);
			int edx = PIDX(GRID_WIDTH-1, y, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = SUM[eidx];
			H[cidx] = H[eidx];
		}
	}
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, -1, z);
			int edx = PIDX(x, 0, z);
			
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = SUM[eidx];
			H[cidx] = H[eidx];
		}
	}
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section ; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, GRID_HEIGHT, z);
			int edx = PIDX(x, GRID_HEIGHT-1, z);

			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = SUM[eidx];
			H[cidx] = H[eidx];
		}
	}

	if (t->id == 0 || SIM_THREADS == 1)
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, -1);
				int edx = PIDX(x, y, 0);

				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM[cidx] = SUM[eidx];
				H[cidx] = H[eidx];
			}
			
		}
	}

	if (t->id == 1 || SIM_THREADS == 1) 
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, GRID_WIDTH);
				int edx = PIDX(x, y, GRID_WIDTH-1);

				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM[cidx] = SUM[eidx];
				H[cidx] = H[eidx];
			}
		}
	}
	
}
// none
void zero_sobel_borders(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);
	
	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(-1, y, z);
			int edx = PIDX(GRID_WIDTH - 1, y, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = 0;
			H[cidx] = 0;
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(GRID_WIDTH, y, z);
			int edx = PIDX(0, y, z);
			
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = 0;
			H[cidx] = 0;
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, -1, z);
			int edx = PIDX(x, GRID_HEIGHT - 1, z);

			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = 0;
			H[cidx] = 0;
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section ; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, GRID_HEIGHT, z);
			int edx = PIDX(x, 0, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			SUM[cidx] = 0;
			H[cidx] = 0;
		}
	}

	if (t->id == 0 || SIM_THREADS == 1)
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, -1);
				int edx = PIDX(x, y, GRID_WIDTH-1);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM[cidx] = 0;
				H[cidx] = 0;
			}
			
		}
	}

	if (t->id == 1 || SIM_THREADS == 1) 
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, GRID_WIDTH);
				int edx = PIDX(x, y, 0);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				SUM[cidx] = 0;
				H[cidx] = 0;
			}
		}
	}
	
}





void transpose_borders(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int dd = 0; dd < (int)off_table.dd; dd++)
		{
			int wi = -1 - dd;
			int we = GRID_WIDTH - 1 - dd;
			for (int z = section; z < end; z++)
			for (int y = 0; y < GRID_HEIGHT; y++)
			{
				int idx = PIDX(wi, y, z);
				int edx = PIDX(we, y, z);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;

				MU_x[cidx] = MU_x[eidx];
				MU_y[cidx] = MU_y[eidx];
				MU_z[cidx] = MU_z[eidx];
				ALPHA[cidx] = ALPHA[eidx];
				V_x[cidx] = V_x[eidx];
				V_y[cidx] = V_y[eidx];
				V_z[cidx] = V_z[eidx];
				F_x[cidx] = F_x[eidx];
				F_y[cidx] = F_y[eidx];
				F_z[cidx] = F_z[eidx];
				GRID[cidx] = GRID[eidx];
				BACK[cidx] = BACK[eidx];
			}
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int dd = 0; dd < (int)off_table.dd; dd++)
		{
			int wi = GRID_WIDTH + dd;
			int we = 0 + dd;
			for (int z = section; z < end; z++)
			for (int y = 0; y < GRID_HEIGHT; y++)
			{
				int idx = PIDX(wi, y, z);
				int edx = PIDX(we, y, z);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;

				MU_x[cidx] = MU_x[eidx];
				MU_y[cidx] = MU_y[eidx];
				MU_z[cidx] = MU_z[eidx];
				ALPHA[cidx] = ALPHA[eidx];
				V_x[cidx] = V_x[eidx];
				V_y[cidx] = V_y[eidx];
				V_z[cidx] = V_z[eidx];
				F_x[cidx] = F_x[eidx];
				F_y[cidx] = F_y[eidx];
				F_z[cidx] = F_z[eidx];
				GRID[cidx] = GRID[eidx];
				BACK[cidx] = BACK[eidx];
			}
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int dd = 0; dd < (int)off_table.dd; dd++)
		{
			int wi = -1 - dd;
			int we = GRID_HEIGHT - 1 - dd;
			for (int z = section; z < end; z++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, wi, z);
				int edx = PIDX(x, we, z);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				MU_x[cidx] = MU_x[eidx];
				MU_y[cidx] = MU_y[eidx];
				MU_z[cidx] = MU_z[eidx];
				ALPHA[cidx] = ALPHA[eidx];
				V_x[cidx] = V_x[eidx];
				V_y[cidx] = V_y[eidx];
				V_z[cidx] = V_z[eidx];
				F_x[cidx] = F_x[eidx];
				F_y[cidx] = F_y[eidx];
				F_z[cidx] = F_z[eidx];
				GRID[cidx] = GRID[eidx];
				BACK[cidx] = BACK[eidx];
			}
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);

		for (int dd = 0; dd < (int)off_table.dd; dd++)
		{
			int wi = GRID_HEIGHT + dd;
			int we = 0 + dd;
			for (int z = section; z < end; z++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, wi, z);
				int edx = PIDX(x, we, z);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				MU_x[cidx] = MU_x[eidx];
				MU_y[cidx] = MU_y[eidx];
				MU_z[cidx] = MU_z[eidx];
				ALPHA[cidx] = ALPHA[eidx];
				V_x[cidx] = V_x[eidx];
				V_y[cidx] = V_y[eidx];
				V_z[cidx] = V_z[eidx];
				F_x[cidx] = F_x[eidx];
				F_y[cidx] = F_y[eidx];
				F_z[cidx] = F_z[eidx];
				GRID[cidx] = GRID[eidx];
				BACK[cidx] = BACK[eidx];
			}
		}
	}

	if (t->id == 0 || SIM_THREADS == 1)
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int dd = 0; dd < (int)off_table.dd; dd++)
			{
				
				int wi = -1 - dd;
				int we = GRID_WIDTH - 1 - dd;

				for (int y = 0; y < GRID_HEIGHT; y++)
				for (int x = 0; x < GRID_WIDTH; x++)
				{
					int idx = PIDX(x, y, wi);
					int edx = PIDX(x, y, we);
				
					int cidx = idx + cpoffset;
					int eidx = edx + cpoffset;
					MU_x[cidx] = MU_x[eidx];
					MU_y[cidx] = MU_y[eidx];
					MU_z[cidx] = MU_z[eidx];
					ALPHA[cidx] = ALPHA[eidx];
					V_x[cidx] = V_x[eidx];
					V_y[cidx] = V_y[eidx];
					V_z[cidx] = V_z[eidx];
					F_x[cidx] = F_x[eidx];
					F_y[cidx] = F_y[eidx];
					F_z[cidx] = F_z[eidx];
					GRID[cidx] = GRID[eidx];
					BACK[cidx] = BACK[eidx];
				}
			}
			
		}
	}

	if (t->id == 1 || SIM_THREADS == 1) 
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);

			for (int dd = 0; dd < (int)off_table.dd; dd++)
			{
				int wi = GRID_WIDTH + dd;
				int we = 0 + dd;
				for (int y = 0; y < GRID_HEIGHT; y++)
				for (int x = 0; x < GRID_WIDTH; x++)
				{
					
					int idx = PIDX(x, y, wi);
					int edx = PIDX(x, y, we);
					
					int cidx = idx + cpoffset;
					int eidx = edx + cpoffset;
					MU_x[cidx] = MU_x[eidx];
					MU_y[cidx] = MU_y[eidx];
					MU_z[cidx] = MU_z[eidx];
					ALPHA[cidx] = ALPHA[eidx];
					V_x[cidx] = V_x[eidx];
					V_y[cidx] = V_y[eidx];
					V_z[cidx] = V_z[eidx];
					F_x[cidx] = F_x[eidx];
					F_y[cidx] = F_y[eidx];
					F_z[cidx] = F_z[eidx];
					GRID[cidx] = GRID[eidx];
					BACK[cidx] = BACK[eidx];
				}
			}
		}
	}
	
}




void transpose_mesh_borders(TD * t)
{
	int section = (GRID_WIDTH / SIM_THREADS) * t->id;
	int end = (t->id == SIM_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / SIM_THREADS);

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(-1, y, z);
			int edx = PIDX(GRID_WIDTH - 1, y, z);
			
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			RBB[cidx] = GRID[eidx];
			H_RB[cidx] = H[eidx];
			U_SUM_RB[cidx] = U_SUM[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int y = 0; y < GRID_HEIGHT; y++)
		{
			int idx = PIDX(GRID_WIDTH, y, z);
			int edx = PIDX(0, y, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			RBB[cidx] = GRID[eidx];
			H_RB[cidx] = H[eidx];
			U_SUM_RB[cidx] = U_SUM[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, -1, z);
			int edx = PIDX(x, GRID_HEIGHT - 1, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			RBB[cidx] = GRID[eidx];
			H_RB[cidx] = H[eidx];
			U_SUM_RB[cidx] = U_SUM[eidx];
		}
	}

	for (int c = 0; c < channels; c++)
	{
		int cpoffset = CPOFFSET(c);
		for (int z = section ; z < end; z++)
		for (int x = 0; x < GRID_WIDTH; x++)
		{
			int idx = PIDX(x, GRID_HEIGHT, z);
			int edx = PIDX(x, 0, z);
		
			int cidx = idx + cpoffset;
			int eidx = edx + cpoffset;
			RBB[cidx] = GRID[eidx];
			H_RB[cidx] = H[eidx];
			U_SUM_RB[cidx] = U_SUM[eidx];
		}
	}

	if (t->id == 0 || SIM_THREADS == 1)
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, -1);
				int edx = PIDX(x, y, GRID_WIDTH-1);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				RBB[cidx] = GRID[eidx];
				H_RB[cidx] = H[eidx];
				U_SUM_RB[cidx] = U_SUM[eidx];
			}
			
		}
	}

	if (t->id == 1 || SIM_THREADS == 1) 
	{
		for (int c = 0; c < channels; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x++)
			{
				int idx = PIDX(x, y, GRID_WIDTH);
				int edx = PIDX(x, y, 0);
			
				int cidx = idx + cpoffset;
				int eidx = edx + cpoffset;
				RBB[cidx] = GRID[eidx];
				H_RB[cidx] = H[eidx];
				U_SUM_RB[cidx] = U_SUM[eidx];
			}
		}
	}
	
}