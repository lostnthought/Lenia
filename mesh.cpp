#include "render.h"
#include "transposition.h"
#include "reintegration.h"
#include "camera.h"
#include "config.h"
#include "lenia.h"
#include "marching_cubes.h"
float mesh_mean;

int V1_EDGES[MESH_THREADS * EDGE_GRID_SIZE * AXIS_MAX];
int V2_EDGES[MESH_THREADS * EDGE_GRID_SIZE * AXIS_MAX];
int V3_EDGES[MESH_THREADS * EDGE_GRID_SIZE * AXIS_MAX];

int LB[MESH_THREADS * EDGE_GRID_SIZE * AXIS_MAX];
int LBI[MESH_THREADS * EDGE_GRID_SIZE * AXIS_MAX];
int UB[MESH_THREADS * EDGE_GRID_SIZE * AXIS_MAX];
int UBI[MESH_THREADS * EDGE_GRID_SIZE * AXIS_MAX];

const __declspec(align(16)) float DELTA_X[] = 
{
	0, 1, 1, 0, 0, 1, 1, 0
};
const __declspec(align(16)) float DELTA_Y[] = 
{
	0, 0, 0, 0, 1, 1, 1, 1
};


const __declspec(align(16)) float DELTA_Z[] = 
{
	0, 0, 1, 1, 0, 0, 1, 1
};

void generate_mesh_v1(TD * t)
{

	D3DXCOLOR color;

	Timer time;
	start_timer(&time);

	int section = (GRID_WIDTH / MESH_THREADS) * t->id;
	int end = (t->id == MESH_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / MESH_THREADS);

	__declspec(align(16)) float POINTS_X[8];
	__declspec(align(16)) float POINTS_Y[8];
	__declspec(align(16)) float POINTS_Z[8];
	__declspec(align(16)) float POINTS_INDICES[8];

	int V1_VAL;
	__declspec(align(16)) float V1_VALS[8];
	__declspec(align(16)) float V1_VERTS_X[12];
	__declspec(align(16)) float V1_VERTS_Y[12];
	__declspec(align(16)) float V1_VERTS_Z[12];

	__declspec(align(16)) float V1_SVERTSA_X[5];
	__declspec(align(16)) float V1_SVERTSA_Y[5];
	__declspec(align(16)) float V1_SVERTSA_Z[5];

	__declspec(align(16)) float V1_SVERTSB_X[5];
	__declspec(align(16)) float V1_SVERTSB_Y[5];
	__declspec(align(16)) float V1_SVERTSB_Z[5];

	__declspec(align(16)) float V1_SVERTSC_X[5];
	__declspec(align(16)) float V1_SVERTSC_Y[5];
	__declspec(align(16)) float V1_SVERTSC_Z[5];

	__declspec(align(16)) float V1_NORM_X[5];
	__declspec(align(16)) float V1_NORM_Y[5];
	__declspec(align(16)) float V1_NORM_Z[5];


	int tid = t->id;
	for (int c = 0; c < channels; c++)
	{	
		color = ch[c].color;

		for (int gz = section; gz < end; gz++){

			for (int gy = 0; gy < GRID_HEIGHT; gy++){

				for (int gx = 0; gx < GRID_WIDTH; gx++)
				{
		
					V1_VAL = 0;

					// MAKE POINTS FROM DELTAS
					__m128 mmx = _mm_add_ps(_mm_set1_ps(gx), _mm_load_ps(&DELTA_X[0]));
					
					__m128 mmy = _mm_add_ps(_mm_set1_ps(gy), _mm_load_ps(&DELTA_Y[0]));
					
					__m128 mmz = _mm_add_ps(_mm_set1_ps(gz), _mm_load_ps(&DELTA_Z[0]));

					

					// MAKE PADDED INDICES
					__m128 ii;
					__m128 chr = _mm_mul_ps(_mm_set1_ps(c), _mm_set1_ps(GRID_PADDED));
					{
						__m128 PADRG = _mm_set1_ps(PADDING);


						ii = _mm_add_ps(
							 _mm_add_ps(
							 _mm_add_ps(mmx, PADRG), 
							 _mm_mul_ps(_mm_add_ps(mmy, PADRG), _mm_set1_ps(PW))), 
							 _mm_mul_ps(_mm_add_ps(mmz, PADRG), _mm_set1_ps(PZOFFSET)));

					}

					ii = _mm_add_ps(ii, chr);

					_mm_store_ps(&POINTS_INDICES[0], ii);
					{
						__m128 CS = _mm_set1_ps(CELL_SIZE);
						mmx = _mm_mul_ps(mmx, CS);
						mmy = _mm_mul_ps(mmy, CS);
						mmz = _mm_mul_ps(mmz, CS);
					}

					_mm_store_ps(&POINTS_X[0], mmx);
					_mm_store_ps(&POINTS_Y[0], mmy);
					_mm_store_ps(&POINTS_Z[0], mmz);

					_mm_store_ps(&V1_VERTS_X[0], mmx);
					_mm_store_ps(&V1_VERTS_Y[0], mmy);
					_mm_store_ps(&V1_VERTS_Z[0], mmz);



					// MAKE POINTS FROM DELTAS
					mmx = _mm_add_ps(_mm_set1_ps(gx), _mm_load_ps(&DELTA_X[4]));
					
					mmy = _mm_add_ps(_mm_set1_ps(gy), _mm_load_ps(&DELTA_Y[4]));
				
					mmz = _mm_add_ps(_mm_set1_ps(gz), _mm_load_ps(&DELTA_Z[4]));

					

					// MAKE PADDED INDICES
					{
						__m128 PADRG = _mm_set1_ps(PADDING);


						ii = _mm_add_ps(
							 _mm_add_ps(
							 _mm_add_ps(mmx, PADRG), 
							 _mm_mul_ps(_mm_add_ps(mmy, PADRG), _mm_set1_ps(PW))), 
							 _mm_mul_ps(_mm_add_ps(mmz, PADRG), _mm_set1_ps(PZOFFSET)));

					}

					ii = _mm_add_ps(ii, chr);

					_mm_store_ps(&POINTS_INDICES[4], ii);

					{
						__m128 CS = _mm_set1_ps(CELL_SIZE);
						mmx = _mm_mul_ps(mmx, CS);
						mmy = _mm_mul_ps(mmy, CS);
						mmz = _mm_mul_ps(mmz, CS);
					}
					
					_mm_store_ps(&POINTS_X[4], mmx);
					_mm_store_ps(&POINTS_Y[4], mmy);
					_mm_store_ps(&POINTS_Z[4], mmz);

					_mm_store_ps(&V1_VERTS_X[4], mmx);
					_mm_store_ps(&V1_VERTS_Y[4], mmy);
					_mm_store_ps(&V1_VERTS_Z[4], mmz);


					for (int i = 0; i < 8; i++){

						V1_VALS[i] = RBB[(int)POINTS_INDICES[i]];
						if (V1_VALS[i] > ISO_VALUE) 
						{
							V1_VAL |= 1 << i;
						}
					}




					int edge = EDGE_TABLE[V1_VAL];
					
				
					for (int i = 0; i < 12; i++)
					{

						if (edge & (1 << i))
						{
							int a = EDGE_TO_VERTICES[i][0];
							int b = EDGE_TO_VERTICES[i][1];

							float T = (ISO_VALUE - V1_VALS[a]) / (V1_VALS[b] - V1_VALS[a]);
							V1_VERTS_X[i] = POINTS_X[a] + (POINTS_X[b] - POINTS_X[a]) * T;
							V1_VERTS_Y[i] = POINTS_Y[a] + (POINTS_Y[b] - POINTS_Y[a]) * T;
							V1_VERTS_Z[i] = POINTS_Z[a] + (POINTS_Z[b] - POINTS_Z[a]) * T;
							

						}
					}

					int j = 0;
					int tris1 = 0;
					int nv = TRIANGLE_TABLE[V1_VAL][j++];
					while (nv != -1)
					{
						V1_SVERTSA_X[tris1] = V1_VERTS_X[nv];
						V1_SVERTSA_Y[tris1] = V1_VERTS_Y[nv];
						V1_SVERTSA_Z[tris1] = V1_VERTS_Z[nv];
						nv = TRIANGLE_TABLE[V1_VAL][j++];

						V1_SVERTSB_X[tris1] = V1_VERTS_X[nv];
						V1_SVERTSB_Y[tris1] = V1_VERTS_Y[nv];
						V1_SVERTSB_Z[tris1] = V1_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V1_VAL][j++];

						V1_SVERTSC_X[tris1] = V1_VERTS_X[nv];
						V1_SVERTSC_Y[tris1] = V1_VERTS_Y[nv];
						V1_SVERTSC_Z[tris1] = V1_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V1_VAL][j++];
						tris1++;
					}

					if (tris1 >= 2) {

						__m128 sva = _mm_load_ps(&V1_SVERTSA_X[0]);
						__m128 svb = _mm_load_ps(&V1_SVERTSB_X[0]);

						__m128 ux = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V1_SVERTSB_Y[0]);

						__m128 uy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V1_SVERTSB_Z[0]);

						__m128 uz = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_X[0]);
						svb = _mm_load_ps(&V1_SVERTSC_X[0]);

						__m128 vx = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V1_SVERTSC_Y[0]);

						__m128 vy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V1_SVERTSC_Z[0]);

						__m128 vz = _mm_sub_ps(svb, sva);

						_mm_store_ps(&V1_NORM_X[0], _mm_sub_ps(_mm_mul_ps(uy, vz), _mm_mul_ps(uz, vy)));
						_mm_store_ps(&V1_NORM_Y[0], _mm_sub_ps(_mm_mul_ps(uz, vx), _mm_mul_ps(ux, vz)));
						_mm_store_ps(&V1_NORM_Z[0], _mm_sub_ps(_mm_mul_ps(ux, vy), _mm_mul_ps(uy, vx)));

						if (tris1 == 5)
						{
							float ux = V1_SVERTSB_X[4] - V1_SVERTSA_X[4];
							float uy = V1_SVERTSB_Y[4] - V1_SVERTSA_Y[4];
							float uz = V1_SVERTSB_Z[4] - V1_SVERTSA_Z[4];

							float vx = V1_SVERTSC_X[4] - V1_SVERTSA_X[4];
							float vy = V1_SVERTSC_Y[4] - V1_SVERTSA_Y[4];
							float vz = V1_SVERTSC_Z[4] - V1_SVERTSA_Z[4];
							V1_NORM_X[4] = uy * vz - uz * vy;
							V1_NORM_Y[4] = uz * vx - ux * vz;
							V1_NORM_Z[4] = ux * vy - uy * vx;
						}
					} 
					else 
					{
						for (int tri = 0; tri < tris1; tri++)
						{
							float ux = V1_SVERTSB_X[tri] - V1_SVERTSA_X[tri];
							float uy = V1_SVERTSB_Y[tri] - V1_SVERTSA_Y[tri];
							float uz = V1_SVERTSB_Z[tri] - V1_SVERTSA_Z[tri];

							float vx = V1_SVERTSC_X[tri] - V1_SVERTSA_X[tri];
							float vy = V1_SVERTSC_Y[tri] - V1_SVERTSA_Y[tri];
							float vz = V1_SVERTSC_Z[tri] - V1_SVERTSA_Z[tri];
							V1_NORM_X[tri] = uy * vz - uz * vy;
							V1_NORM_Y[tri] = uz * vx - ux * vz;
							V1_NORM_Z[tri] = ux * vy - uy * vx;

						}
					}

					int tcnt = 0;
					int next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];
					for (int tri = 0; tri < tris1; tri++)
					{
						
						int p0x = gx + EI[next_vertex].offset.x;
						int p0y = gy + EI[next_vertex].offset.y;
						int p0z = gz + EI[next_vertex].offset.z;
						int eidx = ETIDX(BIDX(p0x, p0y, p0z), EI[next_vertex].axis, tid);
						int oidx = V1_EDGES[eidx];

						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V1_SVERTSA_X[tri], V1_SVERTSA_Y[tri], V1_SVERTSA_Z[tri], color);
							V1_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V1_NORM_X[tri];
						t->v[oidx]._ny += V1_NORM_Y[tri];
						t->v[oidx]._nz += V1_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p0z == section)
						{
							int idx = BIDX(p0x, p0y, p0z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
							//LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;

						}
						else if (p0z == end)
						{
							int idx = BIDX(p0x, p0y, p0z);
							//UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];

						int p1x = gx + EI[next_vertex].offset.x;
						int p1y = gy + EI[next_vertex].offset.y;
						int p1z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p1x, p1y, p1z), EI[next_vertex].axis, tid);

						oidx = V1_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V1_SVERTSB_X[tri], V1_SVERTSB_Y[tri], V1_SVERTSB_Z[tri], color);
							V1_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V1_NORM_X[tri];
						t->v[oidx]._ny += V1_NORM_Y[tri];
						t->v[oidx]._nz += V1_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p1z == section)
						{
							int idx = BIDX(p1x, p1y, p1z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
						}
						else if (p1z == end)
						{
							int idx = BIDX(p1x, p1y, p1z);
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];


						int p2x = gx + EI[next_vertex].offset.x;
						int p2y = gy + EI[next_vertex].offset.y;
						int p2z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p2x, p2y, p2z), EI[next_vertex].axis, tid);
						oidx = V1_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V1_SVERTSC_X[tri], V1_SVERTSC_Y[tri], V1_SVERTSC_Z[tri], color);
							V1_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V1_NORM_X[tri];
						t->v[oidx]._ny += V1_NORM_Y[tri];
						t->v[oidx]._nz += V1_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p2z == section)
						{
							int idx = BIDX(p2x, p2y, p2z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
						}
						else if (p2z == end)
						{
							int idx = BIDX(p2x, p2y, p2z);
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];

					}

				}

			}
		

			

		}

	}

	float diff = end_timer(&time);
	char str[400];
	sprintf_s(str, 399, "THREAD FINISH: %.4f\n", diff);
	OutputDebugStringA(str);
	if (t->id == 0) mesh_mean = (0.9f) * mesh_mean + (1 - 0.9f) * diff;
	sprintf_s(str, 399, "MEAN MESH: %.4f\n", mesh_mean);
	OutputDebugStringA(str);

	t->is_done = true;
}


void generate_mesh_v2(TD * t)
{

	D3DXCOLOR color;

	Timer time;
	start_timer(&time);

	int section = (GRID_WIDTH / MESH_THREADS) * t->id;
	int end = (t->id == MESH_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / MESH_THREADS);

	__declspec(align(16)) float POINTS_X[8];
	__declspec(align(16)) float POINTS_Y[8];
	__declspec(align(16)) float POINTS_Z[8];
	__declspec(align(16)) float POINTS_INDICES[8];

	int V1_VAL;
	__declspec(align(16)) float V1_VALS[8];
	__declspec(align(16)) float V1_VERTS_X[12];
	__declspec(align(16)) float V1_VERTS_Y[12];
	__declspec(align(16)) float V1_VERTS_Z[12];

	__declspec(align(16)) float V1_SVERTSA_X[5];
	__declspec(align(16)) float V1_SVERTSA_Y[5];
	__declspec(align(16)) float V1_SVERTSA_Z[5];

	__declspec(align(16)) float V1_SVERTSB_X[5];
	__declspec(align(16)) float V1_SVERTSB_Y[5];
	__declspec(align(16)) float V1_SVERTSB_Z[5];

	__declspec(align(16)) float V1_SVERTSC_X[5];
	__declspec(align(16)) float V1_SVERTSC_Y[5];
	__declspec(align(16)) float V1_SVERTSC_Z[5];

	__declspec(align(16)) float V1_NORM_X[5];
	__declspec(align(16)) float V1_NORM_Y[5];
	__declspec(align(16)) float V1_NORM_Z[5];

	int V2_VAL;
	__declspec(align(16)) float V2_VALS[8];
	__declspec(align(16)) float V2_VERTS_X[12];
	__declspec(align(16)) float V2_VERTS_Y[12];
	__declspec(align(16)) float V2_VERTS_Z[12];

	__declspec(align(16)) float V2_SVERTSA_X[5];
	__declspec(align(16)) float V2_SVERTSA_Y[5];
	__declspec(align(16)) float V2_SVERTSA_Z[5];

	__declspec(align(16)) float V2_SVERTSB_X[5];
	__declspec(align(16)) float V2_SVERTSB_Y[5];
	__declspec(align(16)) float V2_SVERTSB_Z[5];

	__declspec(align(16)) float V2_SVERTSC_X[5];
	__declspec(align(16)) float V2_SVERTSC_Y[5];
	__declspec(align(16)) float V2_SVERTSC_Z[5];

	__declspec(align(16)) float V2_NORM_X[5];
	__declspec(align(16)) float V2_NORM_Y[5];
	__declspec(align(16)) float V2_NORM_Z[5];

	int tid = t->id;
	for (int c = 0; c < channels; c++)
	{	
		color = ch[c].color;

		for (int gz = section; gz < end; gz++){

			for (int gy = 0; gy < GRID_HEIGHT; gy++){

				for (int gx = 0; gx < GRID_WIDTH; gx++)
				{
		
					V1_VAL = 0; V2_VAL = 0;

					// MAKE POINTS FROM DELTAS
					__m128 mmx = _mm_add_ps(_mm_set1_ps(gx), _mm_load_ps(&DELTA_X[0]));
					
					__m128 mmy = _mm_add_ps(_mm_set1_ps(gy), _mm_load_ps(&DELTA_Y[0]));
					
					__m128 mmz = _mm_add_ps(_mm_set1_ps(gz), _mm_load_ps(&DELTA_Z[0]));

					

					// MAKE PADDED INDICES
					__m128 ii;
					__m128 chr = _mm_mul_ps(_mm_set1_ps(c), _mm_set1_ps(GRID_PADDED));
					{
						__m128 PADRG = _mm_set1_ps(PADDING);


						ii = _mm_add_ps(
							 _mm_add_ps(
							 _mm_add_ps(mmx, PADRG), 
							 _mm_mul_ps(_mm_add_ps(mmy, PADRG), _mm_set1_ps(PW))), 
							 _mm_mul_ps(_mm_add_ps(mmz, PADRG), _mm_set1_ps(PZOFFSET)));

					}

					ii = _mm_add_ps(ii, chr);

					_mm_store_ps(&POINTS_INDICES[0], ii);
					{
						__m128 CS = _mm_set1_ps(CELL_SIZE);
						mmx = _mm_mul_ps(mmx, CS);
						mmy = _mm_mul_ps(mmy, CS);
						mmz = _mm_mul_ps(mmz, CS);
					}

					_mm_store_ps(&POINTS_X[0], mmx);
					_mm_store_ps(&POINTS_Y[0], mmy);
					_mm_store_ps(&POINTS_Z[0], mmz);

					_mm_store_ps(&V1_VERTS_X[0], mmx);
					_mm_store_ps(&V1_VERTS_Y[0], mmy);
					_mm_store_ps(&V1_VERTS_Z[0], mmz);

					_mm_store_ps(&V2_VERTS_X[0], mmx);
					_mm_store_ps(&V2_VERTS_Y[0], mmy);
					_mm_store_ps(&V2_VERTS_Z[0], mmz);


					// MAKE POINTS FROM DELTAS
					mmx = _mm_add_ps(_mm_set1_ps(gx), _mm_load_ps(&DELTA_X[4]));
					
					mmy = _mm_add_ps(_mm_set1_ps(gy), _mm_load_ps(&DELTA_Y[4]));
				
					mmz = _mm_add_ps(_mm_set1_ps(gz), _mm_load_ps(&DELTA_Z[4]));

					

					// MAKE PADDED INDICES
					{
						__m128 PADRG = _mm_set1_ps(PADDING);


						ii = _mm_add_ps(
							 _mm_add_ps(
							 _mm_add_ps(mmx, PADRG), 
							 _mm_mul_ps(_mm_add_ps(mmy, PADRG), _mm_set1_ps(PW))), 
							 _mm_mul_ps(_mm_add_ps(mmz, PADRG), _mm_set1_ps(PZOFFSET)));

					}

					ii = _mm_add_ps(ii, chr);

					_mm_store_ps(&POINTS_INDICES[4], ii);

					{
						__m128 CS = _mm_set1_ps(CELL_SIZE);
						mmx = _mm_mul_ps(mmx, CS);
						mmy = _mm_mul_ps(mmy, CS);
						mmz = _mm_mul_ps(mmz, CS);
					}
					
					_mm_store_ps(&POINTS_X[4], mmx);
					_mm_store_ps(&POINTS_Y[4], mmy);
					_mm_store_ps(&POINTS_Z[4], mmz);

					_mm_store_ps(&V1_VERTS_X[4], mmx);
					_mm_store_ps(&V1_VERTS_Y[4], mmy);
					_mm_store_ps(&V1_VERTS_Z[4], mmz);

					_mm_store_ps(&V2_VERTS_X[4], mmx);
					_mm_store_ps(&V2_VERTS_Y[4], mmy);
					_mm_store_ps(&V2_VERTS_Z[4], mmz);


					for (int i = 0; i < 8; i++){

						V1_VALS[i] = RBB[(int)POINTS_INDICES[i]];
						if (V1_VALS[i] > ISO_VALUE) 
						{
							V1_VAL |= 1 << i;
						}
						
						V2_VALS[i] = U_SUM_RB[(int)POINTS_INDICES[i]];
						if (V2_VALS[i] > ISO_VALUE) 
						{
							V2_VAL |= 1 << i;
						}

					}


					int edge = EDGE_TABLE[V1_VAL];
					
				
					for (int i = 0; i < 12; i++)
					{

						if (edge & (1 << i))
						{
							int a = EDGE_TO_VERTICES[i][0];
							int b = EDGE_TO_VERTICES[i][1];

							float T = (ISO_VALUE - V1_VALS[a]) / (V1_VALS[b] - V1_VALS[a]);
							V1_VERTS_X[i] = POINTS_X[a] + (POINTS_X[b] - POINTS_X[a]) * T;
							V1_VERTS_Y[i] = POINTS_Y[a] + (POINTS_Y[b] - POINTS_Y[a]) * T;
							V1_VERTS_Z[i] = POINTS_Z[a] + (POINTS_Z[b] - POINTS_Z[a]) * T;
							

						}
					}

					int edge2 = EDGE_TABLE[V2_VAL];
					for (int i = 0; i < 12; i++)
					{
						if (edge2 & (1 << i))
						{
							int a = EDGE_TO_VERTICES[i][0];
							int b = EDGE_TO_VERTICES[i][1];
							float T = (ISO_VALUE - V2_VALS[a]) / (V2_VALS[b] - V2_VALS[a]);
							V2_VERTS_X[i] = POINTS_X[a] + (POINTS_X[b] - POINTS_X[a]) * T;
							V2_VERTS_Y[i] = POINTS_Y[a] + (POINTS_Y[b] - POINTS_Y[a]) * T;
							V2_VERTS_Z[i] = POINTS_Z[a] + (POINTS_Z[b] - POINTS_Z[a]) * T;

						}
					}
					__m128 off2_x = _mm_set1_ps(GRID_WIDTH * CELL_SIZE);
					_mm_storeu_ps(&V2_VERTS_X[0], _mm_add_ps(_mm_loadu_ps(&V2_VERTS_X[0]), off2_x));
					_mm_storeu_ps(&V2_VERTS_X[4], _mm_add_ps(_mm_loadu_ps(&V2_VERTS_X[4]), off2_x));
					_mm_storeu_ps(&V2_VERTS_X[8], _mm_add_ps(_mm_loadu_ps(&V2_VERTS_X[8]), off2_x));

					int j = 0;
					int tris1 = 0;
					int nv = TRIANGLE_TABLE[V1_VAL][j++];
					while (nv != -1)
					{
						V1_SVERTSA_X[tris1] = V1_VERTS_X[nv];
						V1_SVERTSA_Y[tris1] = V1_VERTS_Y[nv];
						V1_SVERTSA_Z[tris1] = V1_VERTS_Z[nv];
						nv = TRIANGLE_TABLE[V1_VAL][j++];

						V1_SVERTSB_X[tris1] = V1_VERTS_X[nv];
						V1_SVERTSB_Y[tris1] = V1_VERTS_Y[nv];
						V1_SVERTSB_Z[tris1] = V1_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V1_VAL][j++];

						V1_SVERTSC_X[tris1] = V1_VERTS_X[nv];
						V1_SVERTSC_Y[tris1] = V1_VERTS_Y[nv];
						V1_SVERTSC_Z[tris1] = V1_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V1_VAL][j++];
						tris1++;
					}

					j = 0;
					int tris2 = 0;
					nv = TRIANGLE_TABLE[V2_VAL][j++];
					while (nv != -1)
					{
						V2_SVERTSA_X[tris2] = V2_VERTS_X[nv];
						V2_SVERTSA_Y[tris2] = V2_VERTS_Y[nv];
						V2_SVERTSA_Z[tris2] = V2_VERTS_Z[nv];
						nv = TRIANGLE_TABLE[V2_VAL][j++];

						V2_SVERTSB_X[tris2] = V2_VERTS_X[nv];
						V2_SVERTSB_Y[tris2] = V2_VERTS_Y[nv];
						V2_SVERTSB_Z[tris2] = V2_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V2_VAL][j++];

						V2_SVERTSC_X[tris2] = V2_VERTS_X[nv];
						V2_SVERTSC_Y[tris2] = V2_VERTS_Y[nv];
						V2_SVERTSC_Z[tris2] = V2_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V2_VAL][j++];
						tris2++;
					}

					if (tris1 >= 2) {

						__m128 sva = _mm_load_ps(&V1_SVERTSA_X[0]);
						__m128 svb = _mm_load_ps(&V1_SVERTSB_X[0]);

						__m128 ux = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V1_SVERTSB_Y[0]);

						__m128 uy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V1_SVERTSB_Z[0]);

						__m128 uz = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_X[0]);
						svb = _mm_load_ps(&V1_SVERTSC_X[0]);

						__m128 vx = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V1_SVERTSC_Y[0]);

						__m128 vy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V1_SVERTSC_Z[0]);

						__m128 vz = _mm_sub_ps(svb, sva);

						_mm_store_ps(&V1_NORM_X[0], _mm_sub_ps(_mm_mul_ps(uy, vz), _mm_mul_ps(uz, vy)));
						_mm_store_ps(&V1_NORM_Y[0], _mm_sub_ps(_mm_mul_ps(uz, vx), _mm_mul_ps(ux, vz)));
						_mm_store_ps(&V1_NORM_Z[0], _mm_sub_ps(_mm_mul_ps(ux, vy), _mm_mul_ps(uy, vx)));

						if (tris1 == 5)
						{
							float ux = V1_SVERTSB_X[4] - V1_SVERTSA_X[4];
							float uy = V1_SVERTSB_Y[4] - V1_SVERTSA_Y[4];
							float uz = V1_SVERTSB_Z[4] - V1_SVERTSA_Z[4];

							float vx = V1_SVERTSC_X[4] - V1_SVERTSA_X[4];
							float vy = V1_SVERTSC_Y[4] - V1_SVERTSA_Y[4];
							float vz = V1_SVERTSC_Z[4] - V1_SVERTSA_Z[4];
							V1_NORM_X[4] = uy * vz - uz * vy;
							V1_NORM_Y[4] = uz * vx - ux * vz;
							V1_NORM_Z[4] = ux * vy - uy * vx;
						}
					} 
					else 
					{
						for (int tri = 0; tri < tris1; tri++)
						{
							float ux = V1_SVERTSB_X[tri] - V1_SVERTSA_X[tri];
							float uy = V1_SVERTSB_Y[tri] - V1_SVERTSA_Y[tri];
							float uz = V1_SVERTSB_Z[tri] - V1_SVERTSA_Z[tri];

							float vx = V1_SVERTSC_X[tri] - V1_SVERTSA_X[tri];
							float vy = V1_SVERTSC_Y[tri] - V1_SVERTSA_Y[tri];
							float vz = V1_SVERTSC_Z[tri] - V1_SVERTSA_Z[tri];
							V1_NORM_X[tri] = uy * vz - uz * vy;
							V1_NORM_Y[tri] = uz * vx - ux * vz;
							V1_NORM_Z[tri] = ux * vy - uy * vx;

						}
					}
					if (tris2 >= 2) {

						__m128 sva = _mm_load_ps(&V2_SVERTSA_X[0]);
						__m128 svb = _mm_load_ps(&V2_SVERTSB_X[0]);

						__m128 ux = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V2_SVERTSB_Y[0]);

						__m128 uy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V2_SVERTSB_Z[0]);

						__m128 uz = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_X[0]);
						svb = _mm_load_ps(&V2_SVERTSC_X[0]);

						__m128 vx = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V2_SVERTSC_Y[0]);

						__m128 vy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V2_SVERTSC_Z[0]);

						__m128 vz = _mm_sub_ps(svb, sva);

						_mm_store_ps(&V2_NORM_X[0], _mm_sub_ps(_mm_mul_ps(uy, vz), _mm_mul_ps(uz, vy)));
						_mm_store_ps(&V2_NORM_Y[0], _mm_sub_ps(_mm_mul_ps(uz, vx), _mm_mul_ps(ux, vz)));
						_mm_store_ps(&V2_NORM_Z[0], _mm_sub_ps(_mm_mul_ps(ux, vy), _mm_mul_ps(uy, vx)));

						if (tris1 == 5)
						{
							float ux = V2_SVERTSB_X[4] - V2_SVERTSA_X[4];
							float uy = V2_SVERTSB_Y[4] - V2_SVERTSA_Y[4];
							float uz = V2_SVERTSB_Z[4] - V2_SVERTSA_Z[4];

							float vx = V2_SVERTSC_X[4] - V2_SVERTSA_X[4];
							float vy = V2_SVERTSC_Y[4] - V2_SVERTSA_Y[4];
							float vz = V2_SVERTSC_Z[4] - V2_SVERTSA_Z[4];
							V2_NORM_X[4] = uy * vz - uz * vy;
							V2_NORM_Y[4] = uz * vx - ux * vz;
							V2_NORM_Z[4] = ux * vy - uy * vx;
						}
					} 
					else 
					{
						for (int tri = 0; tri < tris2; tri++)
						{
							float ux = V2_SVERTSB_X[tri] - V2_SVERTSA_X[tri];
							float uy = V2_SVERTSB_Y[tri] - V2_SVERTSA_Y[tri];
							float uz = V2_SVERTSB_Z[tri] - V2_SVERTSA_Z[tri];

							float vx = V2_SVERTSC_X[tri] - V2_SVERTSA_X[tri];
							float vy = V2_SVERTSC_Y[tri] - V2_SVERTSA_Y[tri];
							float vz = V2_SVERTSC_Z[tri] - V2_SVERTSA_Z[tri];
							V2_NORM_X[tri] = uy * vz - uz * vy;
							V2_NORM_Y[tri] = uz * vx - ux * vz;
							V2_NORM_Z[tri] = ux * vy - uy * vx;

						}
					}

					int tcnt = 0;
					int next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];
					for (int tri = 0; tri < tris1; tri++)
					{
						
						int p0x = gx + EI[next_vertex].offset.x;
						int p0y = gy + EI[next_vertex].offset.y;
						int p0z = gz + EI[next_vertex].offset.z;
						int eidx = ETIDX(BIDX(p0x, p0y, p0z), EI[next_vertex].axis, tid);
						int oidx = V1_EDGES[eidx];

						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V1_SVERTSA_X[tri], V1_SVERTSA_Y[tri], V1_SVERTSA_Z[tri], color);
							V1_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V1_NORM_X[tri];
						t->v[oidx]._ny += V1_NORM_Y[tri];
						t->v[oidx]._nz += V1_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p0z == section)
						{
							int idx = BIDX(p0x, p0y, p0z);
							//LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//LBI[ETIDX(t->lower_stitches, EI[next_vertex].axis, tid)] = idx;
							//t->lower_stitches++;
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;

						}
						else if (p0z == end)
						{
							int idx = BIDX(p0x, p0y, p0z);
							//UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//UBI[ETIDX(t->upper_stitches, EI[next_vertex].axis, tid)] = idx;
							//t->upper_stitches++;
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];

						int p1x = gx + EI[next_vertex].offset.x;
						int p1y = gy + EI[next_vertex].offset.y;
						int p1z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p1x, p1y, p1z), EI[next_vertex].axis, tid);

						oidx = V1_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V1_SVERTSB_X[tri], V1_SVERTSB_Y[tri], V1_SVERTSB_Z[tri], color);
							V1_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V1_NORM_X[tri];
						t->v[oidx]._ny += V1_NORM_Y[tri];
						t->v[oidx]._nz += V1_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p1z == section)
						{
							int idx = BIDX(p1x, p1y, p1z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
							//LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//LBI[ETIDX(t->lower_stitches, EI[next_vertex].axis, tid)] = idx;
							//t->lower_stitches++;
						}
						else if (p1z == end)
						{
							int idx = BIDX(p1x, p1y, p1z);
							//UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//UBI[ETIDX(t->upper_stitches, EI[next_vertex].axis, tid)] = idx;
							//t->upper_stitches++;
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];


						int p2x = gx + EI[next_vertex].offset.x;
						int p2y = gy + EI[next_vertex].offset.y;
						int p2z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p2x, p2y, p2z), EI[next_vertex].axis, tid);
						oidx = V1_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V1_SVERTSC_X[tri], V1_SVERTSC_Y[tri], V1_SVERTSC_Z[tri], color);
							V1_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V1_NORM_X[tri];
						t->v[oidx]._ny += V1_NORM_Y[tri];
						t->v[oidx]._nz += V1_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p2z == section)
						{
							int idx = BIDX(p2x, p2y, p2z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
							//LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//LBI[ETIDX(t->lower_stitches, EI[next_vertex].axis, tid)] = idx;
							t->lower_stitches++;
						}
						else if (p2z == end)
						{
							int idx = BIDX(p2x, p2y, p2z);
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
							//UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//UBI[ETIDX(t->upper_stitches, EI[next_vertex].axis, tid)] = idx;
							t->upper_stitches++;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];

					}

					tcnt = 0;
					next_vertex = TRIANGLE_TABLE[V2_VAL][tcnt++];
					for (int tri = 0; tri < tris2; tri++)
					{
						
						int p0x = gx + EI[next_vertex].offset.x;
						int p0y = gy + EI[next_vertex].offset.y;
						int p0z = gz + EI[next_vertex].offset.z;
						int eidx = ETIDX(BIDX(p0x, p0y, p0z), EI[next_vertex].axis, tid);
						int oidx = V2_EDGES[eidx];

						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V2_SVERTSA_X[tri], V2_SVERTSA_Y[tri], V2_SVERTSA_Z[tri], color);
							V2_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V2_NORM_X[tri];
						t->v[oidx]._ny += V2_NORM_Y[tri];
						t->v[oidx]._nz += V2_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p0z == section)
						{
							int idx = BIDX(p0x, p0y, p0z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
							//LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;

						}
						else if (p0z == end)
						{
							int idx = BIDX(p0x, p0y, p0z);
							//UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V2_VAL][tcnt++];

						int p1x = gx + EI[next_vertex].offset.x;
						int p1y = gy + EI[next_vertex].offset.y;
						int p1z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p1x, p1y, p1z), EI[next_vertex].axis, tid);

						oidx = V2_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V2_SVERTSB_X[tri], V2_SVERTSB_Y[tri], V2_SVERTSB_Z[tri], color);
							V2_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V2_NORM_X[tri];
						t->v[oidx]._ny += V2_NORM_Y[tri];
						t->v[oidx]._nz += V2_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p1z == section)
						{
							int idx = BIDX(p1x, p1y, p1z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
							//LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						else if (p1z == end)
						{
							int idx = BIDX(p1x, p1y, p1z);
							//UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V2_VAL][tcnt++];


						int p2x = gx + EI[next_vertex].offset.x;
						int p2y = gy + EI[next_vertex].offset.y;
						int p2z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p2x, p2y, p2z), EI[next_vertex].axis, tid);
						oidx = V2_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V2_SVERTSC_X[tri], V2_SVERTSC_Y[tri], V2_SVERTSC_Z[tri], color);
							V2_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V2_NORM_X[tri];
						t->v[oidx]._ny += V2_NORM_Y[tri];
						t->v[oidx]._nz += V2_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p2z == section)
						{
							int idx = BIDX(p2x, p2y, p2z);
							//LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
						}
						else if (p2z == end)
						{
							int idx = BIDX(p2x, p2y, p2z);
							//UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V2_VAL][tcnt++];

					}
				}

			}
		

			

		}

	}

	float diff = end_timer(&time);
	char str[400];
	sprintf_s(str, 399, "THREAD FINISH: %.4f\n", diff);
	OutputDebugStringA(str);
	if (t->id == 0) mesh_mean = (0.9) * mesh_mean + (1 - 0.9) * diff;
	sprintf_s(str, 399, "MEAN MESH: %.4f\n", mesh_mean);
	OutputDebugStringA(str);

	t->is_done = true;
}





void generate_mesh_v3(TD * t)
{

	D3DXCOLOR color;

	Timer time;
	start_timer(&time);

	int section = (GRID_WIDTH / MESH_THREADS) * t->id;
	int end = (t->id == MESH_THREADS - 1 ? GRID_WIDTH : section + GRID_WIDTH / MESH_THREADS);

	__declspec(align(16)) float POINTS_X[8];
	__declspec(align(16)) float POINTS_Y[8];
	__declspec(align(16)) float POINTS_Z[8];
	__declspec(align(16)) float POINTS_INDICES[8];

	int V1_VAL;
	__declspec(align(16)) float V1_VALS[8];
	__declspec(align(16)) float V1_VERTS_X[12];
	__declspec(align(16)) float V1_VERTS_Y[12];
	__declspec(align(16)) float V1_VERTS_Z[12];

	__declspec(align(16)) float V1_SVERTSA_X[5];
	__declspec(align(16)) float V1_SVERTSA_Y[5];
	__declspec(align(16)) float V1_SVERTSA_Z[5];

	__declspec(align(16)) float V1_SVERTSB_X[5];
	__declspec(align(16)) float V1_SVERTSB_Y[5];
	__declspec(align(16)) float V1_SVERTSB_Z[5];

	__declspec(align(16)) float V1_SVERTSC_X[5];
	__declspec(align(16)) float V1_SVERTSC_Y[5];
	__declspec(align(16)) float V1_SVERTSC_Z[5];

	__declspec(align(16)) float V1_NORM_X[5];
	__declspec(align(16)) float V1_NORM_Y[5];
	__declspec(align(16)) float V1_NORM_Z[5];

	int V2_VAL;
	__declspec(align(16)) float V2_VALS[8];
	__declspec(align(16)) float V2_VERTS_X[12];
	__declspec(align(16)) float V2_VERTS_Y[12];
	__declspec(align(16)) float V2_VERTS_Z[12];

	__declspec(align(16)) float V2_SVERTSA_X[5];
	__declspec(align(16)) float V2_SVERTSA_Y[5];
	__declspec(align(16)) float V2_SVERTSA_Z[5];

	__declspec(align(16)) float V2_SVERTSB_X[5];
	__declspec(align(16)) float V2_SVERTSB_Y[5];
	__declspec(align(16)) float V2_SVERTSB_Z[5];

	__declspec(align(16)) float V2_SVERTSC_X[5];
	__declspec(align(16)) float V2_SVERTSC_Y[5];
	__declspec(align(16)) float V2_SVERTSC_Z[5];

	__declspec(align(16)) float V2_NORM_X[5];
	__declspec(align(16)) float V2_NORM_Y[5];
	__declspec(align(16)) float V2_NORM_Z[5];

	int V3_VAL;
	__declspec(align(16)) float V3_VALS[8];
	__declspec(align(16)) float V3_VERTS_X[12];
	__declspec(align(16)) float V3_VERTS_Y[12];
	__declspec(align(16)) float V3_VERTS_Z[12];

	__declspec(align(16)) float V3_SVERTSA_X[5];
	__declspec(align(16)) float V3_SVERTSA_Y[5];
	__declspec(align(16)) float V3_SVERTSA_Z[5];

	__declspec(align(16)) float V3_SVERTSB_X[5];
	__declspec(align(16)) float V3_SVERTSB_Y[5];
	__declspec(align(16)) float V3_SVERTSB_Z[5];

	__declspec(align(16)) float V3_SVERTSC_X[5];
	__declspec(align(16)) float V3_SVERTSC_Y[5];
	__declspec(align(16)) float V3_SVERTSC_Z[5];

	__declspec(align(16)) float V3_NORM_X[5];
	__declspec(align(16)) float V3_NORM_Y[5];
	__declspec(align(16)) float V3_NORM_Z[5];
	int tid = t->id;
	for (int c = 0; c < channels; c++)
	{	
		color = ch[c].color;

		for (int gz = section; gz < end; gz++){

			for (int gy = 0; gy < GRID_HEIGHT; gy++){

				for (int gx = 0; gx < GRID_WIDTH; gx++)
				{
		
					V1_VAL = 0; V2_VAL = 0; V3_VAL = 0;

					// MAKE POINTS FROM DELTAS
					__m128 mmx = _mm_add_ps(_mm_set1_ps(gx), _mm_load_ps(&DELTA_X[0]));
					
					__m128 mmy = _mm_add_ps(_mm_set1_ps(gy), _mm_load_ps(&DELTA_Y[0]));
					
					__m128 mmz = _mm_add_ps(_mm_set1_ps(gz), _mm_load_ps(&DELTA_Z[0]));

					

					// MAKE PADDED INDICES
					__m128 ii;
					__m128 chr = _mm_mul_ps(_mm_set1_ps(c), _mm_set1_ps(GRID_PADDED));
					{
						__m128 PADRG = _mm_set1_ps(PADDING);


						ii = _mm_add_ps(
							 _mm_add_ps(
							 _mm_add_ps(mmx, PADRG), 
							 _mm_mul_ps(_mm_add_ps(mmy, PADRG), _mm_set1_ps(PW))), 
							 _mm_mul_ps(_mm_add_ps(mmz, PADRG), _mm_set1_ps(PZOFFSET)));

					}

					ii = _mm_add_ps(ii, chr);

					_mm_store_ps(&POINTS_INDICES[0], ii);
					{
						__m128 CS = _mm_set1_ps(CELL_SIZE);
						mmx = _mm_mul_ps(mmx, CS);
						mmy = _mm_mul_ps(mmy, CS);
						mmz = _mm_mul_ps(mmz, CS);
					}

					_mm_store_ps(&POINTS_X[0], mmx);
					_mm_store_ps(&POINTS_Y[0], mmy);
					_mm_store_ps(&POINTS_Z[0], mmz);

					_mm_store_ps(&V1_VERTS_X[0], mmx);
					_mm_store_ps(&V1_VERTS_Y[0], mmy);
					_mm_store_ps(&V1_VERTS_Z[0], mmz);

					_mm_store_ps(&V2_VERTS_X[0], mmx);
					_mm_store_ps(&V2_VERTS_Y[0], mmy);
					_mm_store_ps(&V2_VERTS_Z[0], mmz);

					_mm_store_ps(&V3_VERTS_X[0], mmx);
					_mm_store_ps(&V3_VERTS_Y[0], mmy);
					_mm_store_ps(&V3_VERTS_Z[0], mmz);


					// MAKE POINTS FROM DELTAS
					mmx = _mm_add_ps(_mm_set1_ps(gx), _mm_load_ps(&DELTA_X[4]));
					
					mmy = _mm_add_ps(_mm_set1_ps(gy), _mm_load_ps(&DELTA_Y[4]));
				
					mmz = _mm_add_ps(_mm_set1_ps(gz), _mm_load_ps(&DELTA_Z[4]));

					

					// MAKE PADDED INDICES
					{
						__m128 PADRG = _mm_set1_ps(PADDING);


						ii = _mm_add_ps(
							 _mm_add_ps(
							 _mm_add_ps(mmx, PADRG), 
							 _mm_mul_ps(_mm_add_ps(mmy, PADRG), _mm_set1_ps(PW))), 
							 _mm_mul_ps(_mm_add_ps(mmz, PADRG), _mm_set1_ps(PZOFFSET)));

					}

					ii = _mm_add_ps(ii, chr);

					_mm_store_ps(&POINTS_INDICES[4], ii);

					{
						__m128 CS = _mm_set1_ps(CELL_SIZE);
						mmx = _mm_mul_ps(mmx, CS);
						mmy = _mm_mul_ps(mmy, CS);
						mmz = _mm_mul_ps(mmz, CS);
					}
					
					_mm_store_ps(&POINTS_X[4], mmx);
					_mm_store_ps(&POINTS_Y[4], mmy);
					_mm_store_ps(&POINTS_Z[4], mmz);

					_mm_store_ps(&V1_VERTS_X[4], mmx);
					_mm_store_ps(&V1_VERTS_Y[4], mmy);
					_mm_store_ps(&V1_VERTS_Z[4], mmz);

					_mm_store_ps(&V2_VERTS_X[4], mmx);
					_mm_store_ps(&V2_VERTS_Y[4], mmy);
					_mm_store_ps(&V2_VERTS_Z[4], mmz);

					_mm_store_ps(&V3_VERTS_X[4], mmx);
					_mm_store_ps(&V3_VERTS_Y[4], mmy);
					_mm_store_ps(&V3_VERTS_Z[4], mmz);


					for (int i = 0; i < 8; i++){

						V1_VALS[i] = RBB[(int)POINTS_INDICES[i]];
						if (V1_VALS[i] > ISO_VALUE) 
						{
							V1_VAL |= 1 << i;
						}
						
						V2_VALS[i] = U_SUM_RB[(int)POINTS_INDICES[i]];
						if (V2_VALS[i] > ISO_VALUE) 
						{
							V2_VAL |= 1 << i;
						}

						V3_VALS[i] = H_RB[(int)POINTS_INDICES[i]];
						if (V3_VALS[i] > ISO_VALUE) 
						{
							V3_VAL |= 1 << i;
						}

					}
					//int edge = EDGE_TABLE[V1_VAL];
					//

					//__m128 e1 = _mm_set_ps(edge & 1 << 0, edge & 1 << 1, edge & 1 << 2, edge & 1 << 3);
					//__m128 e2 = _mm_set_ps(edge & 1 << 4, edge & 1 << 5, edge & 1 << 6, edge & 1 << 7);
					//__m128 e3 = _mm_set_ps(edge & 1 << 8, edge & 1 << 9, edge & 1 << 10, edge & 1 << 11);



					int edge = EDGE_TABLE[V1_VAL];
					
				
					for (int i = 0; i < 12; i++)
					{

						if (edge & (1 << i))
						{
							int a = EDGE_TO_VERTICES[i][0];
							int b = EDGE_TO_VERTICES[i][1];

							float T = (ISO_VALUE - V1_VALS[a]) / (V1_VALS[b] - V1_VALS[a]);
							V1_VERTS_X[i] = POINTS_X[a] + (POINTS_X[b] - POINTS_X[a]) * T;
							V1_VERTS_Y[i] = POINTS_Y[a] + (POINTS_Y[b] - POINTS_Y[a]) * T;
							V1_VERTS_Z[i] = POINTS_Z[a] + (POINTS_Z[b] - POINTS_Z[a]) * T;
							

						}
					}

					int edge2 = EDGE_TABLE[V2_VAL];
					for (int i = 0; i < 12; i++)
					{
						if (edge2 & (1 << i))
						{
							int a = EDGE_TO_VERTICES[i][0];
							int b = EDGE_TO_VERTICES[i][1];
							float T = (ISO_VALUE - V2_VALS[a]) / (V2_VALS[b] - V2_VALS[a]);
							V2_VERTS_X[i] = POINTS_X[a] + (POINTS_X[b] - POINTS_X[a]) * T;
							V2_VERTS_Y[i] = POINTS_Y[a] + (POINTS_Y[b] - POINTS_Y[a]) * T;
							V2_VERTS_Z[i] = POINTS_Z[a] + (POINTS_Z[b] - POINTS_Z[a]) * T;

						}
					}
					__m128 off2_x = _mm_set1_ps(GRID_WIDTH * CELL_SIZE);
					_mm_storeu_ps(&V2_VERTS_X[0], _mm_add_ps(_mm_loadu_ps(&V2_VERTS_X[0]), off2_x));
					_mm_storeu_ps(&V2_VERTS_X[4], _mm_add_ps(_mm_loadu_ps(&V2_VERTS_X[4]), off2_x));
					_mm_storeu_ps(&V2_VERTS_X[8], _mm_add_ps(_mm_loadu_ps(&V2_VERTS_X[8]), off2_x));

					int edge3 = EDGE_TABLE[V3_VAL];
					for (int i = 0; i < 12; i++)
					{
						if (edge3 & (1 << i))
						{
							int a = EDGE_TO_VERTICES[i][0];
							int b = EDGE_TO_VERTICES[i][1];
							float T = (ISO_VALUE - V3_VALS[a]) / (V3_VALS[b] - V3_VALS[a]);
							V3_VERTS_X[i] = POINTS_X[a] + (POINTS_X[b] - POINTS_X[a]) * T;
							V3_VERTS_Y[i] = POINTS_Y[a] + (POINTS_Y[b] - POINTS_Y[a]) * T;
							V3_VERTS_Z[i] = POINTS_Z[a] + (POINTS_Z[b] - POINTS_Z[a]) * T;

						}
					}

					_mm_storeu_ps(&V3_VERTS_Z[0], _mm_add_ps(_mm_loadu_ps(&V3_VERTS_Z[0]), off2_x));
					_mm_storeu_ps(&V3_VERTS_Z[4], _mm_add_ps(_mm_loadu_ps(&V3_VERTS_Z[4]), off2_x));
					_mm_storeu_ps(&V3_VERTS_Z[8], _mm_add_ps(_mm_loadu_ps(&V3_VERTS_Z[8]), off2_x));



					int j = 0;
					int tris1 = 0;
					int nv = TRIANGLE_TABLE[V1_VAL][j++];
					while (nv != -1)
					{
						V1_SVERTSA_X[tris1] = V1_VERTS_X[nv];
						V1_SVERTSA_Y[tris1] = V1_VERTS_Y[nv];
						V1_SVERTSA_Z[tris1] = V1_VERTS_Z[nv];
						nv = TRIANGLE_TABLE[V1_VAL][j++];

						V1_SVERTSB_X[tris1] = V1_VERTS_X[nv];
						V1_SVERTSB_Y[tris1] = V1_VERTS_Y[nv];
						V1_SVERTSB_Z[tris1] = V1_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V1_VAL][j++];

						V1_SVERTSC_X[tris1] = V1_VERTS_X[nv];
						V1_SVERTSC_Y[tris1] = V1_VERTS_Y[nv];
						V1_SVERTSC_Z[tris1] = V1_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V1_VAL][j++];
						tris1++;
					}

					j = 0;
					int tris2 = 0;
					nv = TRIANGLE_TABLE[V2_VAL][j++];
					while (nv != -1)
					{
						V2_SVERTSA_X[tris2] = V2_VERTS_X[nv];
						V2_SVERTSA_Y[tris2] = V2_VERTS_Y[nv];
						V2_SVERTSA_Z[tris2] = V2_VERTS_Z[nv];
						nv = TRIANGLE_TABLE[V2_VAL][j++];

						V2_SVERTSB_X[tris2] = V2_VERTS_X[nv];
						V2_SVERTSB_Y[tris2] = V2_VERTS_Y[nv];
						V2_SVERTSB_Z[tris2] = V2_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V2_VAL][j++];

						V2_SVERTSC_X[tris2] = V2_VERTS_X[nv];
						V2_SVERTSC_Y[tris2] = V2_VERTS_Y[nv];
						V2_SVERTSC_Z[tris2] = V2_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V2_VAL][j++];
						tris2++;
					}

					j = 0;
					int tris3 = 0;
					nv = TRIANGLE_TABLE[V3_VAL][j++];
					while (nv != -1)
					{
						V3_SVERTSA_X[tris3] = V3_VERTS_X[nv];
						V3_SVERTSA_Y[tris3] = V3_VERTS_Y[nv];
						V3_SVERTSA_Z[tris3] = V3_VERTS_Z[nv];
						nv = TRIANGLE_TABLE[V3_VAL][j++];

						V3_SVERTSB_X[tris3] = V3_VERTS_X[nv];
						V3_SVERTSB_Y[tris3] = V3_VERTS_Y[nv];
						V3_SVERTSB_Z[tris3] = V3_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V3_VAL][j++];

						V3_SVERTSC_X[tris3] = V3_VERTS_X[nv];
						V3_SVERTSC_Y[tris3] = V3_VERTS_Y[nv];
						V3_SVERTSC_Z[tris3] = V3_VERTS_Z[nv];

						nv = TRIANGLE_TABLE[V3_VAL][j++];
						tris3++;
					}

					if (tris1 >= 2) {

						__m128 sva = _mm_load_ps(&V1_SVERTSA_X[0]);
						__m128 svb = _mm_load_ps(&V1_SVERTSB_X[0]);

						__m128 ux = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V1_SVERTSB_Y[0]);

						__m128 uy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V1_SVERTSB_Z[0]);

						__m128 uz = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_X[0]);
						svb = _mm_load_ps(&V1_SVERTSC_X[0]);

						__m128 vx = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V1_SVERTSC_Y[0]);

						__m128 vy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V1_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V1_SVERTSC_Z[0]);

						__m128 vz = _mm_sub_ps(svb, sva);

						_mm_store_ps(&V1_NORM_X[0], _mm_sub_ps(_mm_mul_ps(uy, vz), _mm_mul_ps(uz, vy)));
						_mm_store_ps(&V1_NORM_Y[0], _mm_sub_ps(_mm_mul_ps(uz, vx), _mm_mul_ps(ux, vz)));
						_mm_store_ps(&V1_NORM_Z[0], _mm_sub_ps(_mm_mul_ps(ux, vy), _mm_mul_ps(uy, vx)));

						if (tris1 == 5)
						{
							float ux = V1_SVERTSB_X[4] - V1_SVERTSA_X[4];
							float uy = V1_SVERTSB_Y[4] - V1_SVERTSA_Y[4];
							float uz = V1_SVERTSB_Z[4] - V1_SVERTSA_Z[4];

							float vx = V1_SVERTSC_X[4] - V1_SVERTSA_X[4];
							float vy = V1_SVERTSC_Y[4] - V1_SVERTSA_Y[4];
							float vz = V1_SVERTSC_Z[4] - V1_SVERTSA_Z[4];
							V1_NORM_X[4] = uy * vz - uz * vy;
							V1_NORM_Y[4] = uz * vx - ux * vz;
							V1_NORM_Z[4] = ux * vy - uy * vx;
						}
					} 
					else 
					{
						for (int tri = 0; tri < tris1; tri++)
						{
							float ux = V1_SVERTSB_X[tri] - V1_SVERTSA_X[tri];
							float uy = V1_SVERTSB_Y[tri] - V1_SVERTSA_Y[tri];
							float uz = V1_SVERTSB_Z[tri] - V1_SVERTSA_Z[tri];

							float vx = V1_SVERTSC_X[tri] - V1_SVERTSA_X[tri];
							float vy = V1_SVERTSC_Y[tri] - V1_SVERTSA_Y[tri];
							float vz = V1_SVERTSC_Z[tri] - V1_SVERTSA_Z[tri];
							V1_NORM_X[tri] = uy * vz - uz * vy;
							V1_NORM_Y[tri] = uz * vx - ux * vz;
							V1_NORM_Z[tri] = ux * vy - uy * vx;

						}
					}
					if (tris2 >= 2) {

						__m128 sva = _mm_load_ps(&V2_SVERTSA_X[0]);
						__m128 svb = _mm_load_ps(&V2_SVERTSB_X[0]);

						__m128 ux = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V2_SVERTSB_Y[0]);

						__m128 uy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V2_SVERTSB_Z[0]);

						__m128 uz = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_X[0]);
						svb = _mm_load_ps(&V2_SVERTSC_X[0]);

						__m128 vx = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V2_SVERTSC_Y[0]);

						__m128 vy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V2_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V2_SVERTSC_Z[0]);

						__m128 vz = _mm_sub_ps(svb, sva);

						_mm_store_ps(&V2_NORM_X[0], _mm_sub_ps(_mm_mul_ps(uy, vz), _mm_mul_ps(uz, vy)));
						_mm_store_ps(&V2_NORM_Y[0], _mm_sub_ps(_mm_mul_ps(uz, vx), _mm_mul_ps(ux, vz)));
						_mm_store_ps(&V2_NORM_Z[0], _mm_sub_ps(_mm_mul_ps(ux, vy), _mm_mul_ps(uy, vx)));

						if (tris1 == 5)
						{
							float ux = V2_SVERTSB_X[4] - V2_SVERTSA_X[4];
							float uy = V2_SVERTSB_Y[4] - V2_SVERTSA_Y[4];
							float uz = V2_SVERTSB_Z[4] - V2_SVERTSA_Z[4];

							float vx = V2_SVERTSC_X[4] - V2_SVERTSA_X[4];
							float vy = V2_SVERTSC_Y[4] - V2_SVERTSA_Y[4];
							float vz = V2_SVERTSC_Z[4] - V2_SVERTSA_Z[4];
							V2_NORM_X[4] = uy * vz - uz * vy;
							V2_NORM_Y[4] = uz * vx - ux * vz;
							V2_NORM_Z[4] = ux * vy - uy * vx;
						}
					} 
					else 
					{
						for (int tri = 0; tri < tris2; tri++)
						{
							float ux = V2_SVERTSB_X[tri] - V2_SVERTSA_X[tri];
							float uy = V2_SVERTSB_Y[tri] - V2_SVERTSA_Y[tri];
							float uz = V2_SVERTSB_Z[tri] - V2_SVERTSA_Z[tri];

							float vx = V2_SVERTSC_X[tri] - V2_SVERTSA_X[tri];
							float vy = V2_SVERTSC_Y[tri] - V2_SVERTSA_Y[tri];
							float vz = V2_SVERTSC_Z[tri] - V2_SVERTSA_Z[tri];
							V2_NORM_X[tri] = uy * vz - uz * vy;
							V2_NORM_Y[tri] = uz * vx - ux * vz;
							V2_NORM_Z[tri] = ux * vy - uy * vx;

						}
					}

					if (tris3 >= 2) {

						__m128 sva = _mm_load_ps(&V3_SVERTSA_X[0]);
						__m128 svb = _mm_load_ps(&V3_SVERTSB_X[0]);

						__m128 ux = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V3_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V3_SVERTSB_Y[0]);

						__m128 uy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V3_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V3_SVERTSB_Z[0]);

						__m128 uz = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V3_SVERTSA_X[0]);
						svb = _mm_load_ps(&V3_SVERTSC_X[0]);

						__m128 vx = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V3_SVERTSA_Y[0]);
						svb = _mm_load_ps(&V3_SVERTSC_Y[0]);

						__m128 vy = _mm_sub_ps(svb, sva);

						sva = _mm_load_ps(&V3_SVERTSA_Z[0]);
						svb = _mm_load_ps(&V3_SVERTSC_Z[0]);

						__m128 vz = _mm_sub_ps(svb, sva);

						_mm_store_ps(&V3_NORM_X[0], _mm_sub_ps(_mm_mul_ps(uy, vz), _mm_mul_ps(uz, vy)));
						_mm_store_ps(&V3_NORM_Y[0], _mm_sub_ps(_mm_mul_ps(uz, vx), _mm_mul_ps(ux, vz)));
						_mm_store_ps(&V3_NORM_Z[0], _mm_sub_ps(_mm_mul_ps(ux, vy), _mm_mul_ps(uy, vx)));

						if (tris1 == 5)
						{
							float ux = V3_SVERTSB_X[4] - V3_SVERTSA_X[4];
							float uy = V3_SVERTSB_Y[4] - V3_SVERTSA_Y[4];
							float uz = V3_SVERTSB_Z[4] - V3_SVERTSA_Z[4];

							float vx = V3_SVERTSC_X[4] - V3_SVERTSA_X[4];
							float vy = V3_SVERTSC_Y[4] - V3_SVERTSA_Y[4];
							float vz = V3_SVERTSC_Z[4] - V3_SVERTSA_Z[4];
							V3_NORM_X[4] = uy * vz - uz * vy;
							V3_NORM_Y[4] = uz * vx - ux * vz;
							V3_NORM_Z[4] = ux * vy - uy * vx;
						}
					} 
					else 
					{
						for (int tri = 0; tri < tris2; tri++)
						{
							float ux = V3_SVERTSB_X[tri] - V3_SVERTSA_X[tri];
							float uy = V3_SVERTSB_Y[tri] - V3_SVERTSA_Y[tri];
							float uz = V3_SVERTSB_Z[tri] - V3_SVERTSA_Z[tri];

							float vx = V3_SVERTSC_X[tri] - V3_SVERTSA_X[tri];
							float vy = V3_SVERTSC_Y[tri] - V3_SVERTSA_Y[tri];
							float vz = V3_SVERTSC_Z[tri] - V3_SVERTSA_Z[tri];
							V3_NORM_X[tri] = uy * vz - uz * vy;
							V3_NORM_Y[tri] = uz * vx - ux * vz;
							V3_NORM_Z[tri] = ux * vy - uy * vx;

						}
					}

					


					


					int tcnt = 0;
					int next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];
					for (int tri = 0; tri < tris1; tri++)
					{
						
						int p0x = gx + EI[next_vertex].offset.x;
						int p0y = gy + EI[next_vertex].offset.y;
						int p0z = gz + EI[next_vertex].offset.z;
						int eidx = ETIDX(BIDX(p0x, p0y, p0z), EI[next_vertex].axis, tid);
						int oidx = V1_EDGES[eidx];

						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V1_SVERTSA_X[tri], V1_SVERTSA_Y[tri], V1_SVERTSA_Z[tri], color);
							V1_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V1_NORM_X[tri];
						t->v[oidx]._ny += V1_NORM_Y[tri];
						t->v[oidx]._nz += V1_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p0z == section)
						{
							int idx = BIDX(p0x, p0y, p0z);
							LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;

						}
						else if (p0z == end)
						{
							int idx = BIDX(p0x, p0y, p0z);
							UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];

						int p1x = gx + EI[next_vertex].offset.x;
						int p1y = gy + EI[next_vertex].offset.y;
						int p1z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p1x, p1y, p1z), EI[next_vertex].axis, tid);

						oidx = V1_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V1_SVERTSB_X[tri], V1_SVERTSB_Y[tri], V1_SVERTSB_Z[tri], color);
							V1_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V1_NORM_X[tri];
						t->v[oidx]._ny += V1_NORM_Y[tri];
						t->v[oidx]._nz += V1_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p1z == section)
						{
							int idx = BIDX(p1x, p1y, p1z);
							LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						else if (p1z == end)
						{
							int idx = BIDX(p1x, p1y, p1z);
							UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];


						int p2x = gx + EI[next_vertex].offset.x;
						int p2y = gy + EI[next_vertex].offset.y;
						int p2z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p2x, p2y, p2z), EI[next_vertex].axis, tid);
						oidx = V1_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V1_SVERTSC_X[tri], V1_SVERTSC_Y[tri], V1_SVERTSC_Z[tri], color);
							V1_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V1_NORM_X[tri];
						t->v[oidx]._ny += V1_NORM_Y[tri];
						t->v[oidx]._nz += V1_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p2z == section)
						{
							int idx = BIDX(p2x, p2y, p2z);
							LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						else if (p2z == end)
						{
							int idx = BIDX(p2x, p2y, p2z);
							UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V1_VAL][tcnt++];

					}

					tcnt = 0;
					next_vertex = TRIANGLE_TABLE[V2_VAL][tcnt++];
					for (int tri = 0; tri < tris2; tri++)
					{
						
						int p0x = gx + EI[next_vertex].offset.x;
						int p0y = gy + EI[next_vertex].offset.y;
						int p0z = gz + EI[next_vertex].offset.z;
						int eidx = ETIDX(BIDX(p0x, p0y, p0z), EI[next_vertex].axis, tid);
						int oidx = V2_EDGES[eidx];

						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V2_SVERTSA_X[tri], V2_SVERTSA_Y[tri], V2_SVERTSA_Z[tri], color);
							V2_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V2_NORM_X[tri];
						t->v[oidx]._ny += V2_NORM_Y[tri];
						t->v[oidx]._nz += V2_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p0z == section)
						{
							int idx = BIDX(p0x, p0y, p0z);
							LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;

						}
						else if (p0z == end)
						{
							int idx = BIDX(p0x, p0y, p0z);
							UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V2_VAL][tcnt++];

						int p1x = gx + EI[next_vertex].offset.x;
						int p1y = gy + EI[next_vertex].offset.y;
						int p1z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p1x, p1y, p1z), EI[next_vertex].axis, tid);

						oidx = V2_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V2_SVERTSB_X[tri], V2_SVERTSB_Y[tri], V2_SVERTSB_Z[tri], color);
							V2_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V2_NORM_X[tri];
						t->v[oidx]._ny += V2_NORM_Y[tri];
						t->v[oidx]._nz += V2_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p1z == section)
						{
							int idx = BIDX(p1x, p1y, p1z);
							LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						else if (p1z == end)
						{
							int idx = BIDX(p1x, p1y, p1z);
							UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V2_VAL][tcnt++];


						int p2x = gx + EI[next_vertex].offset.x;
						int p2y = gy + EI[next_vertex].offset.y;
						int p2z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p2x, p2y, p2z), EI[next_vertex].axis, tid);
						oidx = V2_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V2_SVERTSC_X[tri], V2_SVERTSC_Y[tri], V2_SVERTSC_Z[tri], color);
							V2_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V2_NORM_X[tri];
						t->v[oidx]._ny += V2_NORM_Y[tri];
						t->v[oidx]._nz += V2_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p2z == section)
						{
							int idx = BIDX(p2x, p2y, p2z);
							LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						else if (p2z == end)
						{
							int idx = BIDX(p2x, p2y, p2z);
							UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V2_VAL][tcnt++];

					}

					tcnt = 0;
					next_vertex = TRIANGLE_TABLE[V3_VAL][tcnt++];
					for (int tri = 0; tri < tris3; tri++)
					{
						
						int p0x = gx + EI[next_vertex].offset.x;
						int p0y = gy + EI[next_vertex].offset.y;
						int p0z = gz + EI[next_vertex].offset.z;
						int eidx = ETIDX(BIDX(p0x, p0y, p0z), EI[next_vertex].axis, tid);
						int oidx = V3_EDGES[eidx];

						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V3_SVERTSA_X[tri], V3_SVERTSA_Y[tri], V3_SVERTSA_Z[tri], color);
							V3_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V3_NORM_X[tri];
						t->v[oidx]._ny += V3_NORM_Y[tri];
						t->v[oidx]._nz += V3_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p0z == section)
						{
							int idx = BIDX(p0x, p0y, p0z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
							//LB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//LBI[ETIDX(t->lower_stitches++, EI[next_vertex].axis, tid)] = idx;

						}
						else if (p0z == end)
						{
							int idx = BIDX(p0x, p0y, p0z);
							//UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V3_VAL][tcnt++];

						int p1x = gx + EI[next_vertex].offset.x;
						int p1y = gy + EI[next_vertex].offset.y;
						int p1z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p1x, p1y, p1z), EI[next_vertex].axis, tid);

						oidx = V3_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V3_SVERTSB_X[tri], V3_SVERTSB_Y[tri], V3_SVERTSB_Z[tri], color);
							V3_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V3_NORM_X[tri];
						t->v[oidx]._ny += V3_NORM_Y[tri];
						t->v[oidx]._nz += V3_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p1z == section)
						{
							int idx = BIDX(p1x, p1y, p1z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
						}
						else if (p1z == end)
						{
							int idx = BIDX(p1x, p1y, p1z);
							//UB[ETIDX(idx, EI[next_vertex].axis, tid)] = oidx;
							//UBI[ETIDX(t->upper_stitches++, EI[next_vertex].axis, tid)] = idx;
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V3_VAL][tcnt++];


						int p2x = gx + EI[next_vertex].offset.x;
						int p2y = gy + EI[next_vertex].offset.y;
						int p2z = gz + EI[next_vertex].offset.z;
						eidx = ETIDX(BIDX(p2x, p2y, p2z), EI[next_vertex].axis, tid);
						oidx = V3_EDGES[eidx];
						if (oidx == -1)
						{
							oidx = t->vidx[cur_buf];
							t->v[oidx] = Vertex(V3_SVERTSC_X[tri], V3_SVERTSC_Y[tri], V3_SVERTSC_Z[tri], color);
							V3_EDGES[eidx] = oidx;
							t->vidx[cur_buf]++;
						} 
						t->v[oidx]._nx += V3_NORM_X[tri];
						t->v[oidx]._ny += V3_NORM_Y[tri];
						t->v[oidx]._nz += V3_NORM_Z[tri];
						t->ib[t->iidx[cur_buf]] = oidx;

						if (p2z == section)
						{
							int idx = BIDX(p2x, p2y, p2z);
							t->lower_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->lower_indices[EI[next_vertex].axis].edges[t->lower_stitches++] = idx;
						}
						else if (p2z == end)
						{
							int idx = BIDX(p2x, p2y, p2z);
							t->upper_border[EI[next_vertex].axis].edges[idx] = oidx; 
							t->upper_indices[EI[next_vertex].axis].edges[t->upper_stitches++] = idx;
						}
						
						t->iidx[cur_buf]++;
						next_vertex = TRIANGLE_TABLE[V3_VAL][tcnt++];

					}
				}

			}
		

			

		}

	}

	float diff = end_timer(&time);
	char str[400];
	sprintf_s(str, 399, "THREAD FINISH: %.4f\n", diff);
	OutputDebugStringA(str);
	if (t->id == 0) mesh_mean = (0.9) * mesh_mean + (1 - 0.9) * diff;
	sprintf_s(str, 399, "MEAN MESH: %.4f\n", mesh_mean);
	OutputDebugStringA(str);

	t->is_done = true;
}