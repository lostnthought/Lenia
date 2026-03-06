#include "render.h"
#include "transposition.h"
#include "reintegration.h"
#include "camera.h"
#include "config.h"
#include "mesh.h"
#include "lenia.h"
#include "marching_cubes.h"
#include "string.h"
#include "stdio.h"
#include <intrin.h>
//#include <xmmintrin.h>



vec3 calc_normal(vec3 a, vec3 b, vec3 c)
{
	vec3 u = vec3(b.x - a.x, b.y - a.y, b.z - a.z);
	vec3 v = vec3(c.x - a.x, c.y - a.y, c.z - a.z);
	return vec3(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
}
vec3 normalize(vec3 a)
{
	float m = sqrt(a.x * a.x + a.y * a.y + a.z * a.z); 

	a.x /= m;
	a.y /= m;
	a.z /= m;
	return a;
}

vec3 centroid(vec3 p0, vec3 p1, vec3 p2)
{
	
	return vec3((p0.x + p1.x + p2.x) / 3,
				(p0.y + p1.y + p2.y) / 3,
				(p0.z + p1.z + p2.z) / 3);
}

vec3 avg_normals(vec3 n[3], int cnt)
{
	vec3 a = vec3(0, 0, 0);
	for (int i = 0 ; i < cnt; i++)
	{
		a.x += n[i].x;
		a.y += n[i].y;
		a.z += n[i].z;
	}
	a.x /= cnt;
	a.y /= cnt;
	a.z /= cnt;
	return a;
}


vec3 central_difference(fftwf_complex * c, vec3 p)
{

	/*int offset = (int)GRID_WIDTH / 2 * CELL_SIZE;
	int idx = grid_index(p);
	int x0i = grid_index(round_vec3(vec3(p.x - 1, p.y, p.z)));
	if (!bounds_check(x0i)) x0i = idx;
	int x1i = grid_index(round_vec3(vec3(p.x + 1, p.y, p.z)));
	if (!bounds_check(x1i)) x1i = idx;
	int y0i = grid_index(round_vec3(vec3(p.x, p.y - 1, p.z)));
	if (!bounds_check(y0i)) y0i = idx;
	int y1i = grid_index(round_vec3(vec3(p.x, p.y + 1, p.z)));
	if (!bounds_check(y1i)) y1i = idx;
	int z0i = grid_index(round_vec3(vec3(p.x, p.y, p.z - 1)));
	if (!bounds_check(z0i)) z0i = idx;
	int z1i = grid_index(round_vec3(vec3(p.x, p.y, p.z + 1)));
	if (!bounds_check(z1i)) z1i = idx;
	
	float mx = cells[MAX(x1i, 0)][0] - cells[MAX(x0i, 0)][0];
	float my = cells[MAX(y1i, 0)][0] - cells[MAX(y0i, 0)][0];
	float mz = cells[MAX(z1i, 0)][0] - cells[MAX(z0i, 0)][0];
	
	vec3 a = vec3(mx, my, mz);
	a = normalize(a);*/
	
	return vec3(0, 0, 0);

}

EdgeInfo EI[12];

void precompute_edge_info()
{
	EI[0].axis = AXIS_X;
	EI[0].corner_a = 0;
	EI[0].corner_b = 1;
	EI[0].offset = ivec3(0, 0, 0);

	EI[1].axis = AXIS_Z;
	EI[1].corner_a = 1;
	EI[1].corner_b = 2;
	EI[1].offset = ivec3(1, 0, 0);

	EI[2].axis = AXIS_X;
	EI[2].corner_a = 2;
	EI[2].corner_b = 3;
	EI[2].offset = ivec3(0, 0, 1);

	EI[3].axis = AXIS_Z;
	EI[3].corner_a = 3;
	EI[3].corner_b = 0;
	EI[3].offset = ivec3(0, 0, 0);

	EI[4].axis = AXIS_X;
	EI[4].corner_a = 4;
	EI[4].corner_b = 5;
	EI[4].offset = ivec3(0, 1, 0);

	EI[5].axis = AXIS_Z;
	EI[5].corner_a = 5;
	EI[5].corner_b = 6;
	EI[5].offset = ivec3(1, 1, 0);

	EI[6].axis = AXIS_X;
	EI[6].corner_a = 6;
	EI[6].corner_b = 7;
	EI[6].offset = ivec3(0, 1, 1);

	EI[7].axis = AXIS_Z;
	EI[7].corner_a = 4;
	EI[7].corner_b = 7;
	EI[7].offset = ivec3(0, 1, 0);

	EI[8].axis = AXIS_Y;
	EI[8].corner_a = 0;
	EI[8].corner_b = 4;
	EI[8].offset = ivec3(0, 0, 0);

	EI[9].axis = AXIS_Y;
	EI[9].corner_a = 1;
	EI[9].corner_b = 5;
	EI[9].offset = ivec3(1, 0, 0);

	EI[10].axis = AXIS_Y;
	EI[10].corner_a = 2;
	EI[10].corner_b = 6;
	EI[10].offset = ivec3(1, 0, 1);

	EI[11].axis = AXIS_Y;
	EI[11].corner_a = 3;
	EI[11].corner_b = 7;
	EI[11].offset = ivec3(0, 0, 1);


}


void reset_mesh_generation()
{

	Timer time;
	start_timer(&time);
	int old = cur_buf;

	stitch_borders();
	Vertex * v;
	DWORD * ib;
	cell_verts[old]->Lock(0, 0, (void **) &v, D3DLOCK_DISCARD);
	cvi[old]->Lock(0, 0, (void **) &ib, D3DLOCK_DISCARD);
	voffset[old] = 0;
	ioffset[old] = 0;
	for (int i = 0; i < MESH_THREADS; i++)
	{
		memcpy(v + voffset[old], mesh_td[i].v, sizeof(Vertex) * mesh_td[i].vidx[old]);
 		for (int j = 0; j < mesh_td[i].iidx[old]; j++)
		{
			ib[j + ioffset[old]] = mesh_td[i].ib[j] + voffset[old];
		}
		voffset[old] += mesh_td[i].vidx[old];
		ioffset[old] += mesh_td[i].iidx[old];
	}

	cell_verts[old]->Unlock();
	cvi[old]->Unlock();
	cur_buf ^= 1;
	
	memset(V1_EDGES, -1, sizeof(int) * EDGE_ARR_SIZE);
	memset(V2_EDGES, -1, sizeof(int) * EDGE_ARR_SIZE);
	memset(V3_EDGES, -1, sizeof(int) * EDGE_ARR_SIZE);
	//memset(LB, -1, sizeof(int) * EDGE_ARR_SIZE);
	//memset(LBI, 0, sizeof(int) * EDGE_ARR_SIZE);
	//memset(UB, -1, sizeof(int) * EDGE_ARR_SIZE);
	//memset(UBI, 0, sizeof(int) * EDGE_ARR_SIZE);
	for (int i = 0; i < MESH_THREADS; i++)
	{
		memset(mesh_td[i].lower_border, -1, sizeof(int) * 3 * (GRID_WIDTH + 1) * (GRID_HEIGHT + 1) * (GRID_WIDTH + 1));
		memset(mesh_td[i].lower_indices, 0, sizeof(int) * 3 * (GRID_WIDTH + 1) * (GRID_HEIGHT + 1) * (GRID_WIDTH + 1));
		memset(mesh_td[i].upper_indices, 0, sizeof(int) * 3 * (GRID_WIDTH + 1) * (GRID_HEIGHT + 1) * (GRID_WIDTH + 1));
		memset(mesh_td[i].upper_border, -1, sizeof(int) * 3 * (GRID_WIDTH + 1) * (GRID_HEIGHT + 1) * (GRID_WIDTH + 1));
		
		mesh_td[i].is_done = false;
		mesh_td[i].lower_stitches = 0;
		mesh_td[i].upper_stitches = 0;


		mesh_td[i].vidx[cur_buf] = 0;
		mesh_td[i].iidx[cur_buf] = 0;


	}
	float diff = end_timer(&time);
	char str[400];
	sprintf_s(str, 399, "RESET FINISH: %.4f\n", diff);
	OutputDebugStringA(str);

	
}



void rebuild_kernel_mesh(Vertex * v)
{

	const int delta[][3] = 
	{

		{0, 0, 0},
		{1, 0, 0},
		{1, 0, 1},
		{0, 0, 1},
		{0, 1, 0},
		{1, 1, 0},
		{1, 1, 1},
		{0, 1, 1}
	};
	float offset_w = GRID_WIDTH / 2 * CELL_SIZE;
	float offset_h = GRID_HEIGHT / 2 * CELL_SIZE;

	vec3 offset = vec3(GRID_WIDTH, 0, GRID_WIDTH);

	const D3DXCOLOR colors[MAX_KERNELS] =
	{
		D3DXCOLOR(1.0, 1.0, 1.0, 1.0),
		D3DXCOLOR(1.0, 0.0, 0.0, 1.0),
		D3DXCOLOR(0.0, 1.0, 0.0, 1.0),
		D3DXCOLOR(0.0, 0.0, 1.0, 1.0),
		D3DXCOLOR(0.0, 1.0, 1.0, 1.0),
		D3DXCOLOR(1.0, 1.0, 0.0, 1.0),
		D3DXCOLOR(1.0, 1.0, 1.0, 1.0),
		D3DXCOLOR(0.0, 0.0, 0.0, 1.0),
		D3DXCOLOR(0.0, 0.0, 0.0, 1.0),
		D3DXCOLOR(0.0, 0.0, 0.0, 1.0),
		D3DXCOLOR(0.0, 0.0, 0.0, 1.0),
	};

	float start_x = -GRID_WIDTH / 2 * CELL_SIZE;
	float start_y = -GRID_HEIGHT / 2 * CELL_SIZE;
	float start_z = -GRID_WIDTH / 2 * CELL_SIZE;

	float x = start_x;
	float y = start_y;
	float z = start_z;
	const float K_ISO = 0.00005;
	kernel_vcnt = 0;
	for (int gx = 0; gx < GRID_WIDTH; gx++){
		y = start_y;
		for (int gy = 0; gy < GRID_HEIGHT; gy++){
			z = start_z;
			for (int gz = 0; gz < GRID_WIDTH; gz++)
			{

				if (gx >= GRID_WIDTH / 2 && gz >= GRID_WIDTH / 2)
				{
					continue;
				}
				for (int k = 0; k < kernel_count; k++)
				{
				
					D3DXCOLOR color = colors[k];
					CubeData cube_data;
					cube_data.val = 0;

					for (int i = 0; i < 8; i++){
						
						vec3 p = vec3(clamp(floorf(x + delta[i][0] + offset_w + 0.5f), 0, GRID_WIDTH - 1), 
							clamp(floorf(y + delta[i][1] + offset_h + 0.5f), 0, GRID_HEIGHT - 1), 
							clamp(floorf(z + delta[i][2] + offset_w + 0.5f), 0, GRID_WIDTH - 1));

						cube_data.vals[i] = kernels[k].real_k[grid_index(p)];
						if (cube_data.vals[i] > K_ISO) cube_data.val |= 1 << i;
						cube_data.points[i] = p;

					}


					int edge = EDGE_TABLE[cube_data.val];
				
					for (int i = 0; i < 12; i++)
					{
						cube_data.verts[i] = cube_data.points[i];
						if (edge & (1 << i))
						{
							int a = EDGE_TO_VERTICES[i][0];
							int b = EDGE_TO_VERTICES[i][1];
							cube_data.verts[i] = interpolate_edge(
								cube_data.points[a], 
								cube_data.vals[a], cube_data.points[b], 
								cube_data.vals[b], K_ISO);
						}

						cube_data.verts[i] = vec3_add(cube_data.verts[i], offset);
						
					}
					

					int tcnt = 0;
					int next_vertex = TRIANGLE_TABLE[cube_data.val][tcnt++];
					while (next_vertex != -1)
					{
						D3DXCOLOR lc = color;
						vec3 p0 = cube_data.verts[next_vertex];
						
						next_vertex = TRIANGLE_TABLE[cube_data.val][tcnt++];
						

						vec3 p1 = cube_data.verts[next_vertex];
						next_vertex = TRIANGLE_TABLE[cube_data.val][tcnt++];

						vec3 p2 = cube_data.verts[next_vertex];
						
						vec3 n = calc_normal(p0, p1, p2);
						n = normalize(n);
						vec3 cd = centroid(p0, p1, p2); 
						float d = n.x * forward.x + n.y * forward.y + n.z * forward.z;
						if (d > 0)
						{
							n = vec3(-n.x, -n.y, -n.z);
						}

						v[kernel_vcnt++] = Vertex(p0, n, lc);
						v[kernel_vcnt++] = Vertex(p1, n, lc);
						v[kernel_vcnt++] = Vertex(p2, n, lc);
						next_vertex = TRIANGLE_TABLE[cube_data.val][tcnt++];
					}
				
				}
				z += CELL_SIZE;	
				
			}
			y += CELL_SIZE;
			
		}
		x += CELL_SIZE;
	}


}



volatile bool sim_done = true;
DWORD WINAPI sim_thread_loop(LPVOID lpParam)
{

	Timer conway, sobel, f, mu, reintegration;
	float cdiff = 0, sdiff = 0, fdiff = 0, mdiff = 0, rdiff = 0;
	while (true)
	{
		TD * t = (TD *)lpParam;
		//DWORD dwWaitResult = WaitForMultipleObjects(THREAD_EVENT_MAX, t->tevents, false, INFINITE);
		WaitForSingleObject(sim_work_event, INFINITE);


		//LONG gen = sim_thread_gen;
		if (sim_worker_gen[t->id] == sim_thread_gen || sim_done) continue;
		t->is_done = false;
		sim_worker_gen[t->id] = sim_thread_gen;
		
		
		if (!running) break;
		switch(current_sim_event)
		{
			case SIM_EVENT_FFT:
				{
					
					t->is_done = false;
					

					//start_timer(&conway);	
					start_timer(&conway);
					fft_step(t);
					
				}
				t->is_done = true;
				break;
			case SIM_EVENT_FFT_MT:
				{
					t->is_done = false;

					
					mt_fft_step(t);
					
					
					if (t->id == 0 && GRADIENT_MODE == GRADIENT_CENTRAL_DIFF && SOBEL_BORDER == BORDER_NONE)
					{
						memset(U_S, 0, sizeof(float) * GRID_PADDED * MAX_CHANNELS);
						memset(SUM_S, 0, sizeof(float) * GRID_PADDED * MAX_CHANNELS);
					}
					
				}
				t->is_done = true;
				break;
			case SIM_EVENT_BLUR:
				{
					t->is_done = false;
					if (GRADIENT_MODE == GRADIENT_CENTRAL_DIFF)
					{
						//if (t->id == 0)
						//{
						
							h_sum_blur(t);
						//sum_blur()
						//}

					}
					cdiff = end_timer(&conway);
					t->is_done = true;
				}
				break;
			case SIM_EVENT_TSB:
				{
					t->is_done = false;

					start_timer(&sobel);
					if (GRADIENT_MODE == GRADIENT_SOBEL)
					{
						switch(SOBEL_BORDER)
						{
							case BORDER_TORUS:
								transpose_sobel_borders(t);
								break;
							case BORDER_WALL:
								wall_sobel_borders(t);
								break;
							case BORDER_NONE:
								zero_sobel_borders(t);
								break;
						}
					} 
					else 
					{
						switch(SOBEL_BORDER)
						{
							case BORDER_TORUS:
								transpose_central_diff_borders(t);
								break;
							case BORDER_WALL:
								wall_central_diff_borders(t);
								break;
							case BORDER_NONE:
								//zero_central_diff_borders(t);
								break;
						}

					}
					
				}
				t->is_done = true;
				break;

			case SIM_EVENT_SOBEL_ONE:
				t->is_done = false;
				
				sobel_convolve(t);
				
				t->is_done = true;
				break;
			case SIM_EVENT_SOBEL_TWO:
				t->is_done = false;

				sobel_convolve_h(t);
				sdiff = end_timer(&sobel);
				
				t->is_done = true;
				break;
			case SIM_EVENT_F:
				{
					t->is_done = false;
					

						start_timer(&f);
					
						compute_f(t);
						fdiff = end_timer(&f);

					
				}
				t->is_done = true;
				break;
			case SIM_EVENT_MU:
				{
					t->is_done = false;

						
						
						//compute_velocity_mu_wall(t);
						start_timer(&mu);
						//compute_mu(t);
						//if (t->id == 0)
						//{
					if (MU_VEL)
					{
						if (BORDER == BORDER_WALL)
						{
							compute_velocity_mu_wall(t);
						}
						else
						{
							compute_velocity_mu(t);
						}
					}
					else
					{
						if (BORDER == BORDER_WALL)
						{
							compute_f_mu_wall(t);
						}
						else
						{
							compute_f_mu(t);
						}

					}
							
							//update_velocity(t);
						//}
						
						
					
				}
				t->is_done = true;
				break;
			case SIM_EVENT_TB:
				{
					t->is_done = false;
					transpose_borders(t);
					mdiff = end_timer(&mu);
					

					
				}
				t->is_done = true;
				break;

			case SIM_EVENT_REINTEGRATION:
				t->is_done = false;
				start_timer(&reintegration);
				reintegration_step(t);
			
				t->is_done = true;
				break;
			case SIM_EVENT_TRM:
				{
					t->is_done = false;


					transpose_to_row_major(t);
					rdiff = end_timer(&reintegration);


					
				}
				t->is_done = true;
				break;
			case SIM_EVENT_TMB:
				{
					t->is_done = false;


					transpose_mesh_borders(t);
					
					
					char str[400];
					sprintf_s(str, 399, "CW: %.4f SB: %.4f F: %.4f MU: %.4f RI: %.4f\n", cdiff, sdiff, fdiff, mdiff, rdiff);
					OutputDebugStringA(str);
				}
				t->is_done = true;
				break;
			//case SIM_EVENT_QUIT:
			//	//return 0;
			//	t->is_done = true;
			//	break;
			
			default:
				t->is_done = true;
				break;
		}
		if (InterlockedDecrement(&sim_threads_running) == 0)
		{
			sim_done = true;
			SetEvent(sim_done_event);
			
		}
	}
	
	return 0;


}
volatile bool render_done = true;
DWORD WINAPI mesh_thread_loop(LPVOID lpParam)
{
	while (true)
	{
		TD * t = (TD *)lpParam;
		//DWORD dwWaitResult = WaitForMultipleObjects(THREAD_EVENT_MAX, t->tevents, false, INFINITE);
		WaitForSingleObject(mesh_work_event, INFINITE);


		LONG gen = mesh_thread_gen;
		if (mesh_worker_gen[t->id] == gen || render_done) continue;
		mesh_worker_gen[t->id] = gen;
		t->is_done = false;
		if (!running) break;
		switch(current_mesh_event)
		{
			case MESH_EVENT_MESH:
				//populate_verts(t);
				{
					if (VIEWS == VIEW_NORMAL) 
					{
						generate_mesh_v1(t);
					} 
					else if (VIEWS == VIEW_U)
					{
						generate_mesh_v2(t);
					}
					else if (VIEWS == VIEW_H)
					{
						generate_mesh_v3(t);
					}
				}
				t->is_done = true;
				break;
			case MESH_EVENT_RESET:
				if (t->id == 0)
				{
					reset_mesh_generation();
				}
				t->is_done = true;
				break;
			default:
				break;
		}
		
		
		
		if (InterlockedDecrement(&mesh_threads_running) == 0)
		{
			SetEvent(mesh_done_event);
			render_done = true;
		}
	}
	
	return 0;


}

IDirect3DVertexShader9* VertShader = 0;
IDirect3DPixelShader9* DiffuseShader = 0;
ID3DXConstantTable* VertConstTable = 0;
ID3DXConstantTable* DiffuseConstTable = 0;
D3DXHANDLE ViewMatrixHandle = 0;
D3DXHANDLE ViewProjMatrixHandle = 0;
D3DXHANDLE AmbientMtrlHandle = 0;
D3DXHANDLE DiffuseMtrlHandle = 0;
D3DXHANDLE LightDirHandle = 0;
D3DXHANDLE ViewPositionHandle = 0;
D3DXHANDLE uTimeHandle = 0;
D3DXHANDLE LightColorHandle = 0;

//typedef struct Light
//{
//	float x, y, z, w;
//	float r, g, b, a;
//	Light(){}
//	Light(float _x, float _y, float _z, float _w, float _r, float _g, float _b, float _a)
//	{
//		x = _x; y = _y; z = _z; w = _w; r = _r; g = _g; b = _b; a = _a;
//	}
//
//} Light;
//Light lights[4];
D3DXVECTOR4 light_dir[4];
D3DXVECTOR4 light_color[4];

bool Setup()
{

	ID3DXBuffer* shader = 0;
	ID3DXBuffer* errorBuffer = 0;
	HRESULT hr = D3DXCompileShaderFromFileA(
		"./vert.txt", // shader filename
		0,
		0,
		"Main", // entry point function name
		"vs_2_0", // shader version to compile to
		D3DXSHADER_DEBUG, // debug compile
		&shader,
		&errorBuffer,
		&VertConstTable
	);
		

	if( errorBuffer )
	{
		::MessageBoxA(0, (char*)errorBuffer->GetBufferPointer(), 0, 0);
		d3d::Release<ID3DXBuffer*>(errorBuffer);
	}
	if(FAILED(hr))
	{
		::MessageBoxA(0, "D3DXCreateEffectFromFile() - FAILED", 0, 0);
		return false;
	}

	hr = Device->CreateVertexShader((const DWORD*)shader->GetBufferPointer(), &VertShader);
	if(FAILED(hr))
	{
		::MessageBoxA(0, "D3DXCreateVSHADER - FAILED", 0, 0);
		return false;
	}

	d3d::Release<ID3DXBuffer*>(shader);

	hr = D3DXCompileShaderFromFileA(
		"./diffuse.txt", // shader filename
		0,
		0,
		"Main", // entry point function name
		"ps_2_0", // shader version to compile to
		D3DXSHADER_DEBUG, // debug compile
		&shader,
		&errorBuffer,
		&DiffuseConstTable
	);
	if( errorBuffer )
	{
		::MessageBoxA(0, (char*)errorBuffer->GetBufferPointer(), 0, 0);
		d3d::Release<ID3DXBuffer*>(errorBuffer);
	}
	if(FAILED(hr))
	{
		::MessageBoxA(0, "D3DXCreateEffectFromFifffle() - FAILED", 0, 0);
		return false;
	}

	hr = Device->CreatePixelShader((const DWORD*)shader->GetBufferPointer(), &DiffuseShader);
	if(FAILED(hr))
	{
		::MessageBoxA(0, "D3DXCreateVSHADER - FAILED", 0, 0);
		return false;
	}



	ViewMatrixHandle = DiffuseConstTable->GetConstantByName(
		0, "ViewMatrix");
	ViewProjMatrixHandle = VertConstTable->GetConstantByName(
		0, "ViewProjMatrix");
	AmbientMtrlHandle = DiffuseConstTable->GetConstantByName(
		0, "AmbientMtrl");
	DiffuseMtrlHandle = DiffuseConstTable->GetConstantByName(
		0, "DiffuseMtrl");
	//LightDirHandle = DiffuseConstTable->GetConstantByName(
	//	0, "LightDirection");
	ViewPositionHandle = DiffuseConstTable->GetConstantByName(
		0, "ViewPosition");
	uTimeHandle = DiffuseConstTable->GetConstantByName(
		0, "uTime");
	LightDirHandle = DiffuseConstTable->GetConstantByName(
		0, "light_dir");
	LightColorHandle = DiffuseConstTable->GetConstantByName(
		0, "light_color");

	light_dir[0] = D3DXVECTOR4(-0.5f, 1.0f, 0.0f, 0.0f); light_color[0] = D3DXVECTOR4(1.0f, 0.0f, 1.0f, 1.0f);
	light_dir[1] = D3DXVECTOR4(-0.5f, 0.0f, 0.5f, 0.0f); light_color[1] = D3DXVECTOR4(0.5f, 0.0f, 0.0f, 1.0f);
	light_dir[2] = D3DXVECTOR4(0.5f, -1.0f, 0.0f, 0.0f); light_color[2] = D3DXVECTOR4(0.0f, 0.5f, 0.7f, 1.0f);
	light_dir[3] = D3DXVECTOR4(-0.1f, 0.5f, 0.25f, 0.0f); light_color[3] = D3DXVECTOR4(0.0f, 0.4f, 1.0f, 1.0f);



	//D3DXVECTOR4 directionToLight(-0.57f, 0.57f, -0.57f, 0.0f);
	//DiffuseConstTable->SetVector(Device, LightDirHandle,
	//&directionToLight);

	D3DXVECTOR4 ambientMtrl(0.0f, 0.0f, 1.0f, 1.0f);
	D3DXVECTOR4 diffuseMtrl(0.0f, 0.0f, 1.0f, 1.0f);
	DiffuseConstTable->SetVector(Device,AmbientMtrlHandle,&ambientMtrl);
	DiffuseConstTable->SetVector(Device,DiffuseMtrlHandle,&diffuseMtrl);
	DiffuseConstTable->SetDefaults(Device);
	


	mesh_mean = 0.0f;

	// default channels and colors
	ch[0].asymptotic = false;
	ch[0].soft_clip = false;
	ch[0].color = D3DXCOLOR(D3DCOLOR_XRGB(255, 188, 94));
	ch[1].asymptotic = true;
	ch[1].soft_clip = false;
	ch[1].color = D3DXCOLOR(D3DCOLOR_XRGB(255, 94, 129));
	ch[2].asymptotic = true;
	ch[2].soft_clip = false;
	ch[2].color = D3DXCOLOR(D3DCOLOR_XRGB(105, 255, 94));

	init_fftw_plans();
	

	init_conway();
	Device->SetRenderState(D3DRS_LIGHTING, true);

	Device->CreateVertexBuffer(
		GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH * 15 * MAX_CHANNELS * 3 * sizeof(Vertex),
		D3DUSAGE_WRITEONLY | D3DUSAGE_DYNAMIC,
		Vertex::FVF,
		D3DPOOL_DEFAULT,
		&cell_verts[0],
		0);
	Device->CreateVertexBuffer(
		GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH * 15 * MAX_CHANNELS * 3 * sizeof(Vertex),
		D3DUSAGE_WRITEONLY | D3DUSAGE_DYNAMIC,
		Vertex::FVF,
		D3DPOOL_DEFAULT,
		&cell_verts[1],
		0);
	Device->CreateVertexBuffer(
		GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH * 6 * MAX_KERNELS * sizeof(Vertex),
		D3DUSAGE_WRITEONLY | D3DUSAGE_DYNAMIC,
		Vertex::FVF,
		D3DPOOL_DEFAULT,
		&kernel_verts,
		0);
	Device->CreateIndexBuffer(
		GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH * 15 * MAX_CHANNELS * 3 * sizeof(DWORD),
		D3DUSAGE_WRITEONLY | D3DUSAGE_DYNAMIC,
		D3DFMT_INDEX32,
		D3DPOOL_DEFAULT,
		&cvi[0],
		0);
	Device->CreateIndexBuffer(
		GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH * 15 * MAX_CHANNELS * 3 * sizeof(DWORD),
		D3DUSAGE_WRITEONLY | D3DUSAGE_DYNAMIC,
		D3DFMT_INDEX32,
		D3DPOOL_DEFAULT,
		&cvi[1],
		0);


	// default kernels
	kernels[0].radius = 4;
	kernels[0].beta_size = 3;
	kernels[0].b[0] = 1.0f;
	kernels[0].b[1] = 0.5f;
	kernels[0].b[2] = 0.5f;
	kernels[0].r_a = 1.0f;
	kernels[0].m = 0.3f;
	kernels[0].s = 0.15f;
	kernels[0].h = 1.0f;
	kernels[0].c0 = 0;
	kernels[0].c1 = 0;
	init_kernel(&kernels[kernel_count++]);
	kernels[1].radius = 4;
	kernels[1].beta_size = 3;
	kernels[1].b[0] = 0.5f;
	kernels[1].b[1] = 0.2f;
	kernels[1].b[2] = 0.3f;
	kernels[1].r_a = 0.8f;
	kernels[1].m = 0.23f;
	kernels[1].s = 0.125f;
	kernels[1].h = 1.0f;
	kernels[1].c0 = 1;
	kernels[1].c1 = 1;
	init_kernel(&kernels[kernel_count++]);
	kernels[2].radius = 5;
	kernels[2].beta_size = 2;
	kernels[2].b[0] = 0.5f;
	kernels[2].b[1] = 0.7f;
	kernels[2].b[2] = 0.0;
	kernels[2].r_a = 1.0f;
	kernels[2].m = 0.4f;
	kernels[2].s = 0.2f;
	kernels[2].h = 1.0f;
	kernels[2].c0 = 1;
	kernels[2].c1 = 2;
	init_kernel(&kernels[kernel_count++]);

	blur_kernel.radius = 2;
	blur_kernel.beta_size = 1;
	blur_kernel.b[0] = 1.0f;
	blur_kernel.r_a = 1.0f;
	blur_kernel.m = 0.3f;
	blur_kernel.s = 0.15f;
	blur_kernel.h = 1.0f;
	blur_kernel.c0 = 0;
	blur_kernel.c1 = 0;
	precompute_blur_kernel(&blur_kernel);

	//D3DMATERIAL9 mtrl;
	//::ZeroMemory(&mtrl, sizeof(mtrl));
	//mtrl.Ambient = D3DXCOLOR(0.5, 0.5, 0.5, 1.0);
	//mtrl.Diffuse = D3DXCOLOR(0.5, 0.5, 0.5, 1.0);
	//mtrl.Specular = D3DXCOLOR(1.0, 0.5, 0.5, 1.0);
	//mtrl.Emissive = D3DXCOLOR(D3DCOLOR_XRGB(25, 25, 25));
	//mtrl.Power = 30.0f;
	//Device->SetMaterial(&mtrl);
	
	//
	//::ZeroMemory(&dir, sizeof(dir));
	//dir.Type = D3DLIGHT_DIRECTIONAL;
	//dir.Diffuse = D3DXCOLOR(D3DCOLOR_XRGB(255, 100, 170)) * 0.2f;
	//dir.Specular = D3DXCOLOR(D3DCOLOR_XRGB(255, 100, 170)) * 0.4f;
	//dir.Ambient = D3DXCOLOR(D3DCOLOR_XRGB(255, 100, 170)) * 0.5f;
	//dir.Range = 5.0f;
	//dir.Falloff = 1.0f;
	//dir.Attenuation0 = 0.0f;
	//dir.Attenuation1 = 0.4f;
	//dir.Attenuation2 = 0.5f;

	//dir.Position = D3DXVECTOR3(1.0f, 2.0f, 0.0f);
	//D3DXVECTOR3 d = D3DXVECTOR3(0.0f, -1.0f, 0.2f);
	//D3DXVec3Normalize(&d, &d);
	//dir.Direction = d;
	//Device->SetLight(0, &dir);
	//Device->LightEnable(0, true);

	//::ZeroMemory(&dir, sizeof(dir));
	//dir.Type = D3DLIGHT_DIRECTIONAL;
	//dir.Diffuse = D3DXCOLOR(D3DCOLOR_XRGB(180, 150, 255)) * 0.3f;
	//dir.Specular = D3DXCOLOR(D3DCOLOR_XRGB(180, 150, 255)) * 0.4f;
	//dir.Ambient = D3DXCOLOR(D3DCOLOR_XRGB(180, 150, 255)) * 0.5f;
	//dir.Range = 5.0f;
	//dir.Falloff = 1.0f;
	//dir.Attenuation0 = 0.0f;
	//dir.Attenuation1 = 0.4f;
	//dir.Attenuation2 = 0.5f;

	//dir.Position = D3DXVECTOR3(1.0f, 2.0f, 0.0f);
	//d = D3DXVECTOR3(-0.2f, 1.0f, 0.2f);
	//D3DXVec3Normalize(&d, &d);
	//dir.Direction = d;
	//Device->SetLight(1, &dir);
	//Device->LightEnable(1, true);

	//::ZeroMemory(&dir, sizeof(dir));
	//dir.Type = D3DLIGHT_DIRECTIONAL;
	//dir.Diffuse = D3DXCOLOR(D3DCOLOR_XRGB(100, 255, 190)) * 0.1f;
	//dir.Specular = D3DXCOLOR(D3DCOLOR_XRGB(100, 255, 190)) * 0.2f;
	//dir.Ambient = D3DXCOLOR(D3DCOLOR_XRGB(100, 255, 190)) * 0.3f;
	//dir.Range = 5.0f;
	//dir.Falloff = 1.0f;
	//dir.Attenuation0 = 0.0f;
	//dir.Attenuation1 = 0.4f;
	//dir.Attenuation2 = 0.5f;

	//dir.Position = D3DXVECTOR3(1.0f, 2.0f, 0.0f);
	//d = D3DXVECTOR3(0.4f, 0.5f, -0.4f);
	//D3DXVec3Normalize(&d, &d);
	//dir.Direction = d;
	//Device->SetLight(2, &dir);
	//Device->LightEnable(2, true);

	//::ZeroMemory(&dir, sizeof(dir));
	//dir.Type = D3DLIGHT_DIRECTIONAL;
	//dir.Diffuse = D3DXCOLOR(D3DCOLOR_XRGB(100, 255, 190)) * 0.3f;
	//dir.Specular = D3DXCOLOR(D3DCOLOR_XRGB(100, 255, 190)) * 0.1f;
	//dir.Ambient = D3DXCOLOR(D3DCOLOR_XRGB(100, 255, 190)) * 0.3f;
	//dir.Range = 5.0f;
	//dir.Falloff = 1.0f;
	//dir.Attenuation0 = 0.0f;
	//dir.Attenuation1 = 0.4f;
	//dir.Attenuation2 = 0.5f;

	//dir.Position = D3DXVECTOR3(1.0f, 2.0f, 0.0f);
	//d = D3DXVECTOR3(-0.3f, -0.5f, 0.6f);
	//D3DXVec3Normalize(&d, &d);
	//dir.Direction = d;
	//Device->SetLight(3, &dir);
	//Device->LightEnable(3, true);


	//Device->SetRenderState(D3DRS_NORMALIZENORMALS, true);
	//Device->SetRenderState(D3DRS_SPECULARENABLE, true);
	//Device->SetRenderState(D3DRS_AMBIENT, 0);
	//Device->SetRenderState(D3DRS_ANTIALIASEDLINEENABLE, TRUE);



	update_camera();
	


	return true;
}
void Cleanup()
{

	for (int i = 0; i < SIM_THREADS; i++)
	{
		CloseHandle(sim_threads[i]);
	}
	for (int i = 0; i < MESH_THREADS; i++)
	{
		CloseHandle(mesh_threads[i]);
	}


	d3d::Release<IDirect3DVertexBuffer9*>(cell_verts[0]);
	d3d::Release<IDirect3DVertexBuffer9*>(cell_verts[1]);
	d3d::Release<IDirect3DVertexBuffer9*>(kernel_verts);
	d3d::Release<IDirect3DIndexBuffer9*>(cvi[0]);
	d3d::Release<IDirect3DIndexBuffer9*>(cvi[1]);

	
}
float t = 0.0f;



bool Display()
{
	if (Device)
	{

		t += dt;
		Timer timer;
		start_timer(&timer);
		
		if (redraw_config)
		{
			draw_config();
			redraw_config = false;
		}

		D3DXVECTOR4 v4 = D3DXVECTOR4(position.x, position.y, position.z, 0);
		DiffuseConstTable->SetMatrix(Device, ViewMatrixHandle, &view);
		VertConstTable->SetMatrix(Device, ViewProjMatrixHandle, &viewproj);
		DiffuseConstTable->SetVector(Device, ViewPositionHandle, &v4);
		DiffuseConstTable->SetFloat(Device, uTimeHandle, t);
		DiffuseConstTable->SetVectorArray(Device, LightDirHandle, light_dir, 4);
		DiffuseConstTable->SetVectorArray(Device, LightColorHandle, light_color, 4);

		Device->Clear(0, 0, D3DCLEAR_TARGET | D3DCLEAR_ZBUFFER, D3DXCOLOR(D3DCOLOR_XRGB(240, 240, 255)), 1.0f, 0);
		Device->BeginScene();
		Device->SetFVF(Vertex::FVF);
		
		int d = cur_buf ^ 1;
		
	
		Device->SetVertexShader(VertShader);
		Device->SetPixelShader(DiffuseShader);



		update_camera();
		Device->SetStreamSource(0, cell_verts[d], 0, sizeof(Vertex));
		Device->SetIndices(cvi[d]);

		//Device->SetRenderState(D3DRS_SHADEMODE, D3DSHADE_PHONG);
		//Device->SetRenderState(D3DRS_CULLMODE , D3DCULL_NONE);
		Device->DrawIndexedPrimitive(D3DPT_TRIANGLELIST, 0, 0, voffset[d], 0, ioffset[d] / 3);


		//Device->SetStreamSource(0, kernel_verts, 0, sizeof(Vertex));
		//Device->DrawPrimitive(D3DPT_TRIANGLELIST, 0, kernel_vcnt);
		Device->EndScene();
		Device->Present(0, 0, 0, 0);

	}
	return true;
}

bool d3d::InitD3D(HINSTANCE hInstance,
				  int width, int height,
				  bool windowed,
				  D3DDEVTYPE deviceType,
				  IDirect3DDevice9 **device)
{
	IDirect3D9 * d3d9 = Direct3DCreate9(D3D_SDK_VERSION);
	if (!d3d9) return false;


	// init window
	WNDCLASS wc;

	wc.style = CS_HREDRAW | CS_VREDRAW;
	wc.lpfnWndProc = WndProc;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hInstance = hInstance;
	wc.hIcon = ::LoadIcon(0, IDI_APPLICATION);
	wc.hCursor = ::LoadCursor(0, IDC_ARROW);
	wc.hbrBackground =
	static_cast<HBRUSH>(::GetStockObject(WHITE_BRUSH));
	wc.lpszMenuName = 0;
	wc.lpszClassName = TEXT("Lenia");
	if(!::RegisterClass(&wc))
	{
		::MessageBox(0, TEXT("RegisterClass - Failed"), 0, 0);
		return false;
	}
	HWND hwnd =  ::CreateWindow(
				TEXT("Lenia"),
				TEXT("Lenia"),
				WS_OVERLAPPEDWINDOW,
				CW_USEDEFAULT,
				CW_USEDEFAULT,
				CW_USEDEFAULT,
				CW_USEDEFAULT,
				0,
				0,
				hInstance,
				0);
	if (hwnd == 0)
	{
		return false;
	}
	::ShowWindow(hwnd, true);
	::UpdateWindow(hwnd);


	D3DCAPS9 caps;
	d3d9->GetDeviceCaps(
		D3DADAPTER_DEFAULT,
		deviceType,
		&caps);
	int vp = 0;

	// if we can transform verts on device
	if (caps.DevCaps & D3DDEVCAPS_HWTRANSFORMANDLIGHT)
	{
		vp = D3DCREATE_HARDWARE_VERTEXPROCESSING;
	}
	else 
	{
		vp = D3DCREATE_SOFTWARE_VERTEXPROCESSING;
	}

	D3DPRESENT_PARAMETERS d3dpp;
	d3dpp.BackBufferWidth = width;
	d3dpp.BackBufferHeight = height;
	d3dpp.BackBufferFormat = D3DFMT_A8R8G8B8;
	d3dpp.BackBufferCount = 1;
	d3dpp.MultiSampleType = D3DMULTISAMPLE_NONE;
	d3dpp.MultiSampleQuality = 0;
	d3dpp.SwapEffect = D3DSWAPEFFECT_DISCARD;
	d3dpp.hDeviceWindow = hwnd;
	d3dpp.Windowed = true;
	d3dpp.EnableAutoDepthStencil = true;
	d3dpp.AutoDepthStencilFormat = D3DFMT_D24S8;
	d3dpp.Flags = 0;
	d3dpp.FullScreen_RefreshRateInHz = D3DPRESENT_RATE_DEFAULT;
	d3dpp.PresentationInterval = D3DPRESENT_INTERVAL_IMMEDIATE;

	// create a device using param flags
	HRESULT hr = d3d9->CreateDevice(
		D3DADAPTER_DEFAULT,
		deviceType,
		hwnd,
		vp,
		&d3dpp,
		device);
	if (FAILED(hr)){
		::MessageBox(0, TEXT("CreateDevice() - FAILED"), 0, 0);
		return false;
	}


	return true;
}