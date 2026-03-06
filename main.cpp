#include "render.h"
#include <math.h>
#include "camera.h"
#include "reintegration.h"
#include "lenia.h"
#include "config.h"
#include <time.h>

bool redraw_config = true;
HINSTANCE hinst;
IDirect3DDevice9 * Device = 0;
IDirect3DVertexBuffer9 * cell_verts[2];
IDirect3DVertexBuffer9 * kernel_verts;
IDirect3DIndexBuffer9 * cvi[2];
int kernel_vcnt = 0;
volatile bool halted = false;
int voffset[2];
int ioffset[2];


volatile int running = 0;
int Width = 1280;
int Height = 768;
float move_speed = 1.0f;
float cam_speed = 0.1f;
float dt = 0.0f;
ViewType VIEWS = VIEW_U;
bool BLUR_CENTRAL_DIFF = false;

// hyperparams
float GTA = 0.1f;
float STA = 0.1f;
float V_BETA = 0.9f;
float V_DAMP = 0.9f;
float neighborhood = 50.0f;
float ISO_VALUE = 0.5f;
float N = 2;
float THETA_A = 1;
float timestep = 0.1f;
int channels = 3;
bool ADV_VEL = true;
bool OSC_LA = false;
float OSC_LA_AMT = -0.001f;
bool OSC_FLOW = false;
float OSC_FLOW_AMT = -0.001f;
bool OSC_DIV_FLOW = false;
float OSC_DIV_FLOW_AMT = -0.001f;
bool OSC_VEL = false;
float S_NORM = 0.0f;
bool MU_VEL = false;
Border SOBEL_BORDER = BORDER_NONE;
Border BORDER = BORDER_NONE;
AlphaMode ALPHA_MODE = ALPHA_GRID;
GrowthType GROWTH_TYPE = GROWTH_RT;
ModeType MODE_TYPE = MODE_SRT;
GradientType GRADIENT_MODE = GRADIENT_SOBEL;
Kernel blur_kernel;

Kernel kernels[MAX_KERNELS];
int kernel_count = 0;

OffsetTable off_table;

ChannelData ch[MAX_CHANNELS];
fftwf_plan reverse_fft;
fftwf_plan forward_fft;


// threads
HANDLE sim_done_event;
HANDLE mesh_done_event;
HANDLE sim_work_event;
volatile LONG sim_threads_running;
volatile LONG sim_thread_gen;
volatile LONG sim_worker_gen[SIM_THREADS];
HANDLE mesh_work_event;
volatile LONG mesh_threads_running;
volatile LONG mesh_thread_gen;
volatile LONG mesh_worker_gen[MESH_THREADS];
volatile LONG render_gen;
volatile LONG mesh_render_gen;
SIM_EVENTS current_sim_event;
MESH_EVENTS current_mesh_event;
HANDLE sim_threads[8];
HANDLE mesh_threads[8];
TD sim_td[8];
TD mesh_td[8];


int cur_buf = 0;
D3DLIGHT9 dir;
CRITICAL_SECTION cc;
CRITICAL_SECTION vc;

HFONT font;
HANDLE qu;

HWND cfg_wnd;





// better memory layout wasn't working properly. probably an indexing issue. anyway, it's cursed for now
void stitch_borders()
{
	for (int i = 0; i < MESH_THREADS; i++)
	{
		if (i > 0)
		{
			for (int axis = 0; axis < 3; axis++)
			{
				for (int e = 0; e < mesh_td[i].lower_stitches; e++)
				{
					//int i0 = LB[ETIDX(LBI[ETIDX(e, axis, i)], axis, i)];
					int i0 = mesh_td[i].lower_border[axis].edges[mesh_td[i].lower_indices[axis].edges[e]];
					if (i0 != -1)
					{
						//int i1 = UB[ETIDX(LBI[ETIDX(e, axis, i)], axis, i-1)];
						int i1 = mesh_td[i-1].upper_border[axis].edges[mesh_td[i].lower_indices[axis].edges[e]];
						if (i1 != -1)
						{
							mesh_td[i].v[i0]._nx += mesh_td[i-1].v[i1]._nx;
							mesh_td[i].v[i0]._ny += mesh_td[i-1].v[i1]._ny;
							mesh_td[i].v[i0]._nz += mesh_td[i-1].v[i1]._nz;
							mesh_td[i-1].v[i1]._nx = mesh_td[i].v[i0]._nx;
							mesh_td[i-1].v[i1]._ny = mesh_td[i].v[i0]._ny;
							mesh_td[i-1].v[i1]._nz = mesh_td[i].v[i0]._nz;

						}

					}
				}
			}
			
			}
		
		if (i + 1 < MESH_THREADS)
		{
			for (int axis = 0; axis < 3; axis++)
			{
				for (int e = 0; e < mesh_td[i].upper_stitches; e++)
				{
					//int i0 = LB[ETIDX(UBI[ETIDX(e, axis, i)], axis, i)];
					int i0 = mesh_td[i].lower_border[axis].edges[mesh_td[i].upper_indices[axis].edges[e]];
					if (i0 != -1)
					{
						//int i1 = UB[ETIDX(UBI[ETIDX(e, axis, i)], axis, i+1)];
						int i1 = mesh_td[i+1].upper_border[axis].edges[mesh_td[i].upper_indices[axis].edges[e]];
						if (i1 != -1)
						{
							mesh_td[i].v[i0]._nx += mesh_td[i+1].v[i1]._nx;
							mesh_td[i].v[i0]._ny += mesh_td[i+1].v[i1]._ny;
							mesh_td[i].v[i0]._nz += mesh_td[i+1].v[i1]._nz;
							mesh_td[i+1].v[i1]._nx = mesh_td[i].v[i0]._nx;
							mesh_td[i+1].v[i1]._ny = mesh_td[i].v[i0]._ny;
							mesh_td[i+1].v[i1]._nz = mesh_td[i].v[i0]._nz;

						}

					}
					
				}
			}

		}
	}

			
}

bool transpose_override = false;
void check_and_dispatch_threads()
{


	int d = 0;
	for (int i = 0; i < SIM_THREADS; i++)
	{
		if (sim_td[i].is_done) d++;
	}
	if (sim_done && !halted)
	{
		
		ResetEvent(sim_done_event);
		if (halted)
		{
			SetEvent(sim_done_event);
			return;
		}

		if (current_sim_event == SIM_EVENT_TMB && !transpose_override)
		{
			InterlockedIncrement(&render_gen);
		}
		transpose_override = false;
		
		current_sim_event = (SIM_EVENTS)(((int)current_sim_event + 1) % (int)SIM_EVENT_MAX);
		sim_threads_running = SIM_THREADS;
		InterlockedIncrement(&sim_thread_gen);
		SetEvent(sim_work_event);
		sim_done = false;

	}
	

	d = 0;
	for (int i = 0; i < MESH_THREADS; i++)
	{
		if (mesh_td[i].is_done) d++;
	}
	if (render_done && !halted)
	{
		
		ResetEvent(mesh_done_event);
		if (halted)
		{
			SetEvent(mesh_done_event);
			return;
		}
		MESH_EVENTS next_mesh_event = (MESH_EVENTS)(((int)current_mesh_event + 1) % (int)MESH_EVENT_MAX);
		// no new mesh ready
		if (next_mesh_event == MESH_EVENT_MESH && mesh_render_gen == render_gen) return;
		
		mesh_render_gen = render_gen;
		for (int c = 0; c < MAX_CHANNELS; c++)
		{
			int cpoffset = CPOFFSET(c);
			for (int z = 0; z < GRID_WIDTH; z++)
			for (int y = 0; y < GRID_HEIGHT; y++)
			for (int x = 0; x < GRID_WIDTH; x+=4)
			{
				int cidx = PIDX(x, y, z) + cpoffset;
				_mm_store_ps(&RBB[cidx], _mm_load_ps(&GRID[cidx]));
			}
		}
		
		
		current_mesh_event = (MESH_EVENTS)(((int)current_mesh_event + 1) % (int)MESH_EVENT_MAX);
		mesh_threads_running = MESH_THREADS;
		InterlockedIncrement(&mesh_thread_gen);
		SetEvent(mesh_work_event);
		render_done = false;
	}


}

void setup_new_thread(TD * t, int i)
{
	t->id = i;
	

	t->is_done = true;

	t->bsize = (GRID_WIDTH / MESH_THREADS) * 
		(GRID_HEIGHT) * 
		(GRID_WIDTH) * 15 * MAX_CHANNELS * 3;
	t->boffset = 0;

	t->vidx[0] = 0;
	t->vidx[1] = 0;
	t->iidx[0] = 0;
	t->iidx[1] = 0;
	t->v = (Vertex*)calloc(t->bsize, sizeof(Vertex));
	t->ib = (DWORD*)calloc(t->bsize, sizeof(DWORD));

}

int d3d::EnterMsgLoop(bool (*ptr_display)())
{

	// compile shader


	MSG msg;
	::ZeroMemory(&msg, sizeof(MSG));
	srand(time(NULL));
	static float lastTime = (float)timeGetTime();

	int i = 0;
	int nf = 0;
	
	float bounds[8][6];

	float step = (GRID_WIDTH / MESH_THREADS) * CELL_SIZE;
	float start = (-GRID_WIDTH / 2) * CELL_SIZE;
	for (int i = 0; i < MESH_THREADS; i++)
	{
		bounds[i][0] = start;
		bounds[i][1] = start + step;
		bounds[i][2] = (-GRID_HEIGHT / 2) * CELL_SIZE;
		bounds[i][3] = (GRID_HEIGHT / 2) * CELL_SIZE ;
		bounds[i][4] = (-GRID_WIDTH / 2) * CELL_SIZE;
		bounds[i][5] = (GRID_WIDTH / 2) * CELL_SIZE;

		start += step;
	}


	fftwf_complex a;
	init_pos_array();

	off_table.dd = 1.0f; off_table.sigma = 0.5; 
	off_table.X = NULL; off_table.Y = NULL; off_table.Z = NULL;
	init_gaussian_kernel();
	init_offset_table(NULL);
	init_bell_lut();

	precompute_edge_info();

	sim_work_event = CreateEventA(NULL, TRUE, FALSE, "WORK");
	mesh_work_event = CreateEventA(NULL, TRUE, FALSE, "WORK");
	sim_done_event = CreateEventA(NULL, TRUE, FALSE, "DONE");
	mesh_done_event = CreateEventA(NULL, TRUE, FALSE, "DONE");

	sim_thread_gen = 0;
	mesh_thread_gen = 0;
	mesh_render_gen = 0;
	render_gen = 0;
	current_sim_event = SIM_EVENT_FFT;
	current_mesh_event = MESH_EVENT_MESH;
	for (int i = 0; i < SIM_THREADS; i++)
	{
		sim_worker_gen[i] = 0;
		setup_new_thread(&sim_td[i], i);

		sim_threads[i] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)sim_thread_loop, (void*)&sim_td[i], 0, NULL);

	}
	for (int i = 0; i < MESH_THREADS; i++)
	{
		setup_new_thread(&mesh_td[i], i);
		mesh_worker_gen[i] = 0;

		mesh_threads[i] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)mesh_thread_loop, (void*)&mesh_td[i], 0, NULL);

	}

	Timer render;
	while (msg.message != WM_QUIT)
	{
		if (::PeekMessage(&msg, 0, 0, 0, PM_REMOVE))
		{
		
			::TranslateMessage(&msg);
			::DispatchMessage(&msg);
		
		} else {
		
			float currTime = (float)timeGetTime();
			dt = (currTime - lastTime) * 0.001f;
			

			ptr_display();
			lastTime = currTime;
		}
		check_and_dispatch_threads();
		nf++;
	}
	return msg.wParam;


}

void handle_key_input(HWND hwnd, WPARAM wParam){

	vec3 proposed = vec3(0.0f, 0.0f, 0.0f);
	switch(wParam){
		case VK_ESCAPE:
			running = false;
			::DestroyWindow(hwnd);
			break;
		default:
			if (wParam == 'a'){
				proposed.z += move_speed;
			} else if (wParam == 'd'){
				proposed.z -= move_speed;
			} else if (wParam == 'w'){
				proposed.x += move_speed;
			} else if (wParam == 's'){
				proposed.x -= move_speed;
			} else if (wParam == 'r'){
				WaitForSingleObject(sim_done_event, INFINITE);
				WaitForSingleObject(mesh_done_event, INFINITE);
				generate_random();
				current_sim_event = (SIM_EVENTS)((int)SIM_EVENT_MAX - 1);
				transpose_override = true;
			} else if (wParam == VK_LEFT) {
				yaw += cam_speed;
			} else if (wParam == VK_RIGHT) {
				yaw -= cam_speed;	
			} else if (wParam == VK_UP) {
				pitch += cam_speed;
			} else if (wParam == VK_DOWN) {
				pitch -= cam_speed;
			} else if (wParam == VK_SPACE) {
				proposed.y += move_speed;
			} else if (wParam == VK_CONTROL) {
				proposed.y -= move_speed;
			} else if (wParam == VK_TAB) {
				WaitForSingleObject(sim_done_event, INFINITE);
				WaitForSingleObject(mesh_done_event, INFINITE);
				init_conway();
				current_sim_event = (SIM_EVENTS)((int)SIM_EVENT_MAX - 1);
				transpose_override = true;
			} 
			break;
	}
	move(proposed);
}

LRESULT CALLBACK d3d::WndProc(HWND hwnd, UINT msg, WPARAM wParam,
	LPARAM lParam)
{
	switch( msg )
	{
		case WM_DESTROY:
			::PostQuitMessage(0);
			break;
		case WM_KEYDOWN:
			handle_key_input(hwnd, wParam);
			break;
		case WM_CHAR:
			handle_key_input(hwnd, wParam);
			break;
	}
	return ::DefWindowProc(hwnd, msg, wParam, lParam);
}
int WINAPI WinMain(HINSTANCE hinstance,
					HINSTANCE prevInstance,
					PSTR cmdLine,
					int showCmd)
{
	running = 0;
	hinst = hinstance;
	if(!d3d::InitD3D(hinstance, Width, Height, true, D3DDEVTYPE_HAL, &Device))
	{
		::MessageBox(0, TEXT("InitD3D() - FAILED"), 0, 0);
		return 0;
	}
	if (!config_window_init(hinstance, 700, 100, false))
	{
		::MessageBox(0, TEXT("CWINDOW - FAILED"), 0, 0);
		return 0;
	}
	if(!Setup())
	{
		::MessageBox(0, TEXT("Setup() - FAILED"), 0, 0);
		return 0;
	}
	running = 1;
	d3d::EnterMsgLoop( Display );
	running = 0;
	Cleanup();
	Device->Release();
	return 0;
}







