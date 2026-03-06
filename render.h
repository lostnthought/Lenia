#ifndef RENDER_H
#define RENDER_H

#include <d3dx9.h>
#include <windows.h>
#include <math.h>
#include <fftw3.h>
#include <xmmintrin.h>

extern D3DXCOLOR c0;
extern D3DXCOLOR c1;
extern D3DXCOLOR c2;
#define MAX_THREADS 4
#define SIM_THREADS 2
#define MESH_THREADS 2
extern float neighborhood;

extern int voffset[2]; 
extern int ioffset[2];
extern float GTA;
extern float STA;
extern float V_BETA;
extern float V_DAMP;
extern bool BLUR_CENTRAL_DIFF;

#define CELL_SIZE 3.5f
#define GRID_WIDTH 32
#define GRID_HEIGHT 32
#define GRID_SIZE (GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH)
#define HALF_WIDTH (GRID_WIDTH / 2)
#define HALF_HEIGHT (GRID_HEIGHT / 2)
#define ZOFFSET GRID_WIDTH * GRID_HEIGHT
#define MAX_KERNELS 12
#define MAX_CHANNELS 4
#define MAX_BETA 12
#define BELL_LUT_SIZE 256 * 256 * 256


extern IDirect3DDevice9 * Device;
extern IDirect3DVertexBuffer9 * cell_verts[2];
extern IDirect3DVertexBuffer9 * kernel_verts;
extern IDirect3DIndexBuffer9 * cvi[2];
extern int kernel_vcnt;

extern D3DXVECTOR3 forward;
extern D3DLIGHT9 dir;
extern int Width;
extern int Height;
extern float move_speed;
extern float yaw;
extern float pitch;
extern float dt;
extern D3DXVECTOR3 forward;
extern D3DXVECTOR3 right;
extern int cur_buf;
extern int channels;
extern CRITICAL_SECTION cc;
extern CRITICAL_SECTION vc;
extern HWND cfg_wnd;
extern float timestep;
extern float ISO_VALUE;
extern bool redraw_config;
extern HINSTANCE hinst;
extern HFONT font;

extern D3DXMATRIX viewproj;
extern D3DXMATRIX view;
extern D3DXVECTOR3 position;



extern volatile int running;

extern volatile bool halted;
extern HANDLE sim_threads[8];
extern HANDLE mesh_threads[8];
extern volatile LONG sim_threads_running;
extern volatile LONG sim_thread_gen;
extern volatile LONG sim_worker_gen[SIM_THREADS];
extern volatile LONG render_gen;
extern volatile LONG mesh_render_gen;
extern HANDLE sim_done_event;
extern HANDLE mesh_done_event;
extern HANDLE sim_work_event;
extern volatile LONG mesh_threads_running;
extern volatile LONG mesh_thread_gen;
extern volatile LONG mesh_worker_gen[MESH_THREADS];
extern HANDLE mesh_work_event;
extern bool render_lock;
extern volatile bool sim_done;
extern volatile bool render_done;



typedef enum AXIS
{
	AXIS_X,
	AXIS_Y,
	AXIS_Z,
	AXIS_MAX
} AXIS;

typedef struct ivec3 {
	
	ivec3(){}
	ivec3(int _x, int _y, int _z)
	{
		x = _x; y = _y; z = _z;
		
	}
	int x, y, z;
} ivec3;

typedef struct EdgeInfo
{
	ivec3 offset;
	int corner_a;
	int corner_b;
	AXIS axis;


} EdgeInfo;


extern float N;
extern float THETA_A;

typedef enum AlphaMode
{
	ALPHA_GRID, // A
	ALPHA_GRID_GROWTH, // A + H
	ALPHA_GA, // A moving average
	ALPHA_GA_DIFF, // abs(A - A moving average)
	ALPHA_MAX

} AlphaMode;
extern AlphaMode ALPHA_MODE;

typedef enum GrowthType
{
	// in stages, prealpha (before F), rt (during reintegration), transposition (after rt)
	GROWTH_NOGROWTH,
	GROWTH_RT,
	GROWTH_TRANSPOSITION,
	GROWTH_MAX
} GrowthType;

extern GrowthType GROWTH_TYPE;

typedef enum ModeType
{
	MODE_SRT, // SINGULAR REINTEGRATION
	MODE_RT, // REINTEGRATION
	MODE_ADV, // ADVECTION
	MODE_MAX
} ModeType;

extern ModeType MODE_TYPE;
extern bool ADV_VEL;
extern bool MU_VEL;
extern float S_NORM;

extern bool OSC_LA;
extern float OSC_LA_AMT;
extern bool OSC_FLOW;
extern float OSC_FLOW_AMT;
extern bool OSC_DIV_FLOW;
extern float OSC_DIV_FLOW_AMT;
extern bool OSC_VEL;

typedef enum ThreadEvents
{
	THREAD_EVENT_SOBEL_ONE,
	THREAD_EVENT_SOBEL_TWO,
	THREAD_EVENT_RENDER,
	THREAD_EVENT_REINTEGRATION,
	THREAD_EVENT_QUIT,
	THREAD_EVENT_MAX
} ThreadEvents;

typedef enum SIM_EVENTS
{
	SIM_EVENT_FFT,
	SIM_EVENT_FFT_MT,
	SIM_EVENT_BLUR,
	SIM_EVENT_TSB,
	SIM_EVENT_SOBEL_ONE,
	SIM_EVENT_SOBEL_TWO,
	SIM_EVENT_F,
	SIM_EVENT_MU,
	SIM_EVENT_TB,
	SIM_EVENT_REINTEGRATION,
	SIM_EVENT_TRM,
	SIM_EVENT_TMB,

	SIM_EVENT_MAX
} SIM_EVENTS;
typedef enum MESH_EVENTS
{
	MESH_EVENT_MESH,
	MESH_EVENT_RESET,
	MESH_EVENT_MAX
} MESH_EVENTS;

extern SIM_EVENTS current_sim_event;
extern MESH_EVENTS current_mesh_event;
typedef enum Border
{
	BORDER_TORUS,
	BORDER_WALL,
	BORDER_NONE,
	BORDER_MAX
}Border;
extern Border SOBEL_BORDER;
extern Border BORDER;

typedef enum GradientType
{
	GRADIENT_CENTRAL_DIFF,
	GRADIENT_SOBEL,
	GRADIENT_MAX
} GradientType;
extern GradientType GRADIENT_MODE;

typedef enum ViewType
{
	VIEW_NORMAL,
	VIEW_U,
	VIEW_H,
	VIEW_MAX
} ViewType;

extern ViewType VIEWS;

typedef struct vec3 {
	float x, y, z;
	vec3(){}
	vec3(float _x, float _y, float _z)
	{
		x = _x; y = _y; z = _z;
		
	}
} vec3;



typedef struct Timer {

	LARGE_INTEGER freq;
	LARGE_INTEGER start;
	LARGE_INTEGER end;

} Timer;

static inline void start_timer (Timer * timer)
{
	QueryPerformanceFrequency(&timer->freq);
	QueryPerformanceCounter(&timer->start);
}
static inline float end_timer (Timer * timer)
{
	QueryPerformanceCounter(&timer->end);
	return (float)((timer->end.QuadPart - timer->start.QuadPart)) / timer->freq.QuadPart;
}



extern EdgeInfo EI[12];
extern __declspec(align(16)) float gaussian_kernel[3][3][3];

extern __declspec(align(16)) float POS_ARRAY_X[GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH];
extern __declspec(align(16)) float POS_ARRAY_Y[GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH];
extern __declspec(align(16)) float POS_ARRAY_Z[GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH];



// unpadded

extern __declspec(align(16)) float GA[GRID_SIZE * MAX_CHANNELS];



extern __declspec(align(16)) float NA_x[GRID_SIZE];
extern __declspec(align(16)) float NA_y[GRID_SIZE];
extern __declspec(align(16)) float NA_z[GRID_SIZE];

extern __declspec(align(16)) float NU_x[GRID_SIZE * MAX_CHANNELS];
extern __declspec(align(16)) float NU_y[GRID_SIZE * MAX_CHANNELS];
extern __declspec(align(16)) float NU_z[GRID_SIZE * MAX_CHANNELS];

#define PADDING 4

#define GRID_PADDED (GRID_WIDTH + (PADDING * 2)) * (GRID_HEIGHT + (PADDING * 2)) * (GRID_WIDTH + (PADDING * 2))


// padding
extern __declspec(align(16)) float SUM_S[GRID_PADDED];
extern __declspec(align(16)) float U_S[GRID_PADDED * MAX_CHANNELS];

extern __declspec(align(16)) float ALPHA[GRID_PADDED * MAX_CHANNELS];

extern __declspec(align(16)) float SUM[GRID_PADDED];

extern __declspec(align(16)) float V_x[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float V_y[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float V_z[GRID_PADDED * MAX_CHANNELS];

extern __declspec(align(16)) float F_x[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float F_y[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float F_z[GRID_PADDED * MAX_CHANNELS];

extern __declspec(align(16)) float MU_x[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float MU_y[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float MU_z[GRID_PADDED * MAX_CHANNELS];

extern __declspec(align(16)) float GRID[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float BACK[GRID_PADDED * MAX_CHANNELS];

extern __declspec(align(16)) float RBB[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float RB[GRID_PADDED * MAX_CHANNELS];

extern __declspec(align(16)) float H_RB[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float H[GRID_PADDED * MAX_CHANNELS];

extern __declspec(align(16)) float U_SUM_RB[GRID_PADDED * MAX_CHANNELS];
extern __declspec(align(16)) float U_SUM[GRID_PADDED * MAX_CHANNELS];

extern __declspec(align(16)) float TT[GRID_PADDED * MAX_CHANNELS];

#define PW (GRID_WIDTH + PADDING * 2)
#define PH (GRID_HEIGHT + PADDING * 2)
#define PZOFFSET PW * PH

#define IDX(x, y, z) ((x) + (y * GRID_WIDTH) + (z * ZOFFSET))
#define PIDX(x, y, z) ((x + PADDING) + ((y+PADDING) * PW) + ((z+PADDING) * PZOFFSET))
#define COFFSET(c) (c * GRID_SIZE)
#define CPOFFSET(c) (c * GRID_PADDED)
#define CIDX(i, c) ((c * GRID_SIZE) + (i))
#define CPIDX(i, c) ((c * GRID_PADDED) + (i))
#define TTIOFFSET SIM_THREADS * MAX_CHANNELS
#define TTIDX(t, i, c) (t * GRID_PADDED * MAX_CHANNELS + c * GRID_PADDED + i)


#define B8(b7, b6, b5, b4, b3, b2, b1, b0)\
	((b7<<7)|(b6<<6)|(b5<<5)|(b4<<4)|(b3<<3)|(b2<<2)|(b1<<1)|(b0<<0))\

typedef struct Kernel{
	fftwf_complex grid[GRID_SIZE];
	fftwf_complex kg[GRID_SIZE];
	float real_k[GRID_SIZE];
	int beta_size;
	float b[MAX_BETA];
	float r_a;
	float m;
	float s;
	float h;
	int radius;
	int c0; // input channel
	int c1; // output channel
	
} Kernel;


typedef struct CubeData{

	int val;
	float vals[8];
	vec3 points[8];
	vec3 verts[12];

} CubeData;


typedef struct OffsetTable {

	float sigma;
	float dd;
	int row_size;
	float ma;
	int * X;
	int * Y;
	int * Z;

	int offsets[GRID_WIDTH * GRID_SIZE];

} OffsetTable;

extern OffsetTable off_table;



namespace d3d
{
	bool InitD3D(
		HINSTANCE hInstance,
		int width, int height,
		bool windowed,

		D3DDEVTYPE deviceType,
		IDirect3DDevice9** device);

	int EnterMsgLoop(
		bool (*ptr_display)());

	LRESULT CALLBACK WndProc(
		HWND hwnd,
		UINT msg,
		WPARAM wParam,
		LPARAM lParam);

	template<class T> void Release(T t)
	{
		if (t)
		{
			t->Release();
			t = 0;
		}
	}

	template<class T> void Delete(T t)
	{
		if (t)
		{
			delete t;
			t = 0;
		}
	}

	const D3DXCOLOR GREEN( D3DCOLOR_XRGB( 0, 255, 0) );
	const D3DXCOLOR BLUE( D3DCOLOR_XRGB( 0, 0, 255) );
	const D3DXCOLOR YELLOW( D3DCOLOR_XRGB(255, 255, 0) );
	const D3DXCOLOR CYAN( D3DCOLOR_XRGB( 0, 255, 255) );
	const D3DXCOLOR MAGENTA( D3DCOLOR_XRGB(255, 0, 255) );

}

bool Setup();
void Cleanup();
bool Display();



typedef struct Vertex
{
	Vertex(){}
	Vertex(vec3 a)
	{
		_x = a.x; _y = a.y; _z = a.z;
		_nx = 0.0f; _ny = 1.0f; _nz = 0.0f;
		color = D3DCOLOR_XRGB(255, 255, 255);
	}
	Vertex(float x, float y, float z)
	{
		_x = x; _y = y; _z = z;
		_nx = 0.0f; _ny = 1.0f; _nz = 0.0f;
	}
	Vertex(vec3 a, vec3 n, D3DXCOLOR c)
	{
		_x = a.x; _y = a.y; _z = a.z;
		_nx = n.x; _ny = n.y; _nz = n.z;
		color = c;
	}
	Vertex(vec3 a, D3DXCOLOR c)
	{
		_x = a.x; _y = a.y; _z = a.z;
		_nx = 0; _ny = 0; _nz = 0;
		color = c;
	}
	Vertex(float x, float y, float z, D3DXCOLOR c)
	{
		_x = x; _y = y; _z = z;
		_nx = 0; _ny = 0; _nz = 0;
		color = c;
	}
	Vertex(float x, float y, float z, float nx, float ny, float nz)
	{
		_x = x; _y = y; _z = z;
		_nx = 0.0f; _ny = 1.0f; _nz = 0.0f;
	}
	float _x, _y, _z;
	float _nx, _ny, _nz;
	D3DCOLOR color;
	static const DWORD FVF = D3DFVF_XYZ | D3DFVF_NORMAL | D3DFVF_DIFFUSE;
} Vertex;

typedef struct ChannelData
{
	bool asymptotic;
	bool soft_clip;
	D3DXCOLOR color;

}ChannelData;


typedef struct EdgeArr
{
	int edges[GRID_WIDTH + 1][GRID_HEIGHT + 1][GRID_WIDTH + 1];
} EdgeArr;

typedef struct BorderArr
{
	int edges[(GRID_WIDTH + 1) * (GRID_HEIGHT + 1) * (GRID_WIDTH + 1)];
} BorderArr;

typedef struct EA
{
	// axis
	EdgeArr arr[3];
} EA;

#define EDGE_PADDING 2
#define EW GRID_WIDTH + EDGE_PADDING
#define EH GRID_HEIGHT + EDGE_PADDING
#define EDGE_GRID_SIZE (EW) * (EH) * (EW)

#define EIDX(x, y, z) (x + y * GRID_WIDTH + z * GRID_WIDTH * GRID_HEIGHT)
#define ETIDX(i, a, t) (t * EDGE_GRID_SIZE * AXIS_MAX + a * EDGE_GRID_SIZE + i)
#define EDGE_ARR_SIZE (MESH_THREADS * EDGE_GRID_SIZE * AXIS_MAX)

extern int V1_EDGES[EDGE_ARR_SIZE];
extern int V2_EDGES[EDGE_ARR_SIZE];
extern int V3_EDGES[EDGE_ARR_SIZE];

#define BIDX(x, y, z) (x + y * (GRID_WIDTH + 1) + z * (GRID_WIDTH + 1) * (GRID_HEIGHT + 1))

extern int LB[EDGE_ARR_SIZE];
extern int LBI[EDGE_ARR_SIZE];
extern int UB[EDGE_ARR_SIZE];
extern int UBI[EDGE_ARR_SIZE];

vec3 normalize(vec3 a);
typedef struct TD
{
	float bounds[6];
	int vidx[2];
	int iidx[2];
	Vertex * v;
	DWORD * ib;
	bool is_done;
	int id;
	HANDLE tevents[THREAD_EVENT_MAX];
	int bsize;
	int boffset;

	// view
	EA edges[3];


	BorderArr lower_border[3];
	int lower_stitches;
	BorderArr lower_indices[3];
	BorderArr upper_border[3];
	int upper_stitches;
	BorderArr upper_indices[3];

} TD;

struct Button;

typedef enum DataTypes
{
	DT_FLOAT,
	DT_INT,
	DT_BOOL,
	DT_D3DXCOLOR,
	DT_ENUM,
	DT_NULL
} DataTypes;

typedef enum ActionDataTypes
{

	ADT_KERNEL,
	ADT_CHANNEL,


} ActionDataTypes;


typedef struct Button
{
	int cid;
	void (*fn)(Button*);
	HWND box_id;
	void * data;
	DataTypes dtype;
	bool redraw;

	void (*action_fn)(Button *);
	void * action_data;
	ActionDataTypes adt;
	union {
		int ex_int;
		float ex_float;
	};

	

} Buttons;

typedef struct TextBox
{
	int cid;

} TextBox;


static inline int grid_index(vec3 p){

	return (int)(p.x / CELL_SIZE + p.y / CELL_SIZE * GRID_WIDTH + p.z / CELL_SIZE * GRID_HEIGHT * GRID_WIDTH);
}

static inline bool bounds_check(int idx)
{
	if (idx < 0 || idx >= GRID_WIDTH * GRID_HEIGHT * GRID_WIDTH) return false;
	return true;
}

static inline vec3 vec3_sub(vec3 a, vec3 b)
{
	return vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline vec3 vec3_add_s(vec3 a, float b)
{
	return vec3(a.x + b, a.y + b, a.z + b);
}
static inline vec3 vec3_add(vec3 a, vec3 b)
{
	return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline vec3 vec3_mul(vec3 a, vec3 b)
{
	return vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

static inline vec3 vec3_scale(vec3 a, float b)
{
	return vec3(a.x * b, a.y * b, a.z * b);
}

static inline int emod(int a, int b)
{

	int r = a % b;
	return (r < 0) ? r + b : r;

	//if (a < 0)
	//{
	//	return b + a;
	//}
	//else if (a > b)
	//{
	//	return a - b;
	//}
	//return a;
}

static inline float emodf(float a, float b)
{

	float r = fmodf(a, b);
	return (r < 0) ? r + b : r;
}

#define MIN(a, b) \
	(a < b) ? (a) : (b) \

#define MAX(a, b) \
	(a > b) ? (a) : (b) \


static inline float bell(float x, float m, float s)
{
	return exp(-pow((x-m)/s, 2)/2);
}

static inline float soft_clip(float i, float min, float max)
{
	return 1 / (1 + (float)exp(-4 * (i - 0.5)));
}

static inline float clamp(float a, float b, float c){
	if (a < b) a = b;
	if (a > c) a = c;
	return a;
}

static inline vec3 world_to_cell(vec3 a)
{

	a.x += GRID_WIDTH * CELL_SIZE / 2;
	a.y += GRID_WIDTH * CELL_SIZE / 2;
	a.z += GRID_HEIGHT * CELL_SIZE / 2;
	return a;
}
static inline vec3 interpolate_edge(vec3 a, float va, vec3 b, float vb, float iso)
{
	float t = (iso - va) / (vb - va);
	vec3 pos = vec3_add(a, vec3_scale(vec3_sub(b, a), t));
	return pos;
}
static inline int normal_index(vec3 a)
{

	return (int)(((a.z * GRID_WIDTH * GRID_HEIGHT) + (a.y * GRID_WIDTH) + (a.x)) * 1000);
}

static inline __m128 abs_ps(__m128 x)
{
	static const __m128 sign_mask = _mm_set1_ps(-0.f);
	return _mm_andnot_ps(sign_mask, x);
}

void rebuild_kernel_mesh(Vertex * v);
void populate_verts(TD * td);
void precompute_edge_info();
void init_mesh_registers();

extern Kernel kernels[MAX_KERNELS];
extern Kernel blur_kernel;
extern int kernel_count;
extern ChannelData ch[MAX_CHANNELS];
extern fftwf_plan reverse_fft;
extern fftwf_plan forward_fft;
extern TD mesh_td[8];
extern TD sim_td[8];

DWORD WINAPI sim_thread_loop(LPVOID lpParam);
DWORD WINAPI mesh_thread_loop(LPVOID lpParam);
void stitch_borders();

#endif