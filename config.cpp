#include "render.h"
#include "windows.h"
#include "lenia.h"
#include <time.h>

#define MAX_BUTTONS 400
#define MAX_TEXT_BOX 400
#define MAX_COMPONENTS 800
Button buttons[MAX_BUTTONS];
TextBox boxes[MAX_TEXT_BOX];
HWND components[MAX_COMPONENTS];
int component_cnt = 0;
int box_cnt = 0;
int button_cnt = 0;
int u_cid = 0;
float scale = 0.1f;

const char * BOOL_VAL[] =
{
	"FALSE", "TRUE"
};
const char * BORDER_TYPES[] =
{
	"BORDER_TORUS",
	"BORDER_WALL",
	"BORDER_NONE"
};

const char * VIEW_TYPES[] =
{
	"VIEW_NORMAL",
	"VIEW_U",
	"VIEW_H"
};

const char * GRADIENT_TYPES[] =
{
	"GRADIENT_CENTRAL_DIFF",
	"GRADIENT_SOBEL"
};


const char * ALPHA_TYPES[] =
{
	"ALPHA_GRID",
	"ALPHA_GRID_GROWTH",
	"ALPHA_GA",
	"ALPHA_GA_DIFF"
};

const char * GROWTH_TYPES[] =
{
	"GROWTH_NOGROWTH",
	"GROWTH_RT",
	"GROWTH_TRANSPOSITION",
};

const char * MODE_TYPES[] =
{
	"MODE_SRT",
	"MODE_RT",
	"MODE_ADV",
};


void destroy_kernel(Kernel * k){

	//if (k->grid) fftwf_free(k->grid);
	//k->grid = NULL;
	//if (k->kg) fftwf_free(k->kg);
	//k->kg = NULL;
}


void save_config(Button * button)
{
	//char path[128];
	//GetWindowTextA(button->box_id, path, 122); 

	//const char * bi = ".bin";

	//strncat(path, bi, 127);
	//FILE * f = NULL;
	//fopen_s(&f, path, "wb");

	//fwrite(&ISO_VALUE, sizeof(float), 1, f);
	//fwrite(&timestep, sizeof(float), 1, f);
	//fwrite(&channels, sizeof(int), 1, f);
	//fwrite(&GTA, sizeof(float), 1, f);
	//fwrite(&STA, sizeof(float), 1, f);
	//fwrite(&V_BETA, sizeof(float), 1, f);
	//fwrite(&V_DAMP, sizeof(float), 1, f);
	//fwrite(&neighborhood, sizeof(float), 1, f);
	//fwrite(&N, sizeof(float), 1, f);
	//fwrite(&THETA_A, sizeof(float), 1, f);
	//fwrite(&ADV_VEL, sizeof(bool), 1, f);

	//fwrite(&OSC_LA, sizeof(bool), 1, f);
	//fwrite(&OSC_LA_AMT, sizeof(float), 1, f);

	//fwrite(&OSC_FLOW, sizeof(bool), 1, f);
	//fwrite(&OSC_FLOW_AMT, sizeof(float), 1, f);

	//fwrite(&OSC_DIV_FLOW, sizeof(bool), 1, f);
	//fwrite(&OSC_DIV_FLOW_AMT, sizeof(float), 1, f);

	//fwrite(&OSC_VEL, sizeof(bool), 1, f);
	//fwrite(&ADV_VEL, sizeof(bool), 1, f);
	//fwrite(&MU_VEL, sizeof(bool), 1, f);

	//fwrite(&BORDER, sizeof(int), 1, f);
	//fwrite(&SOBEL_BORDER, sizeof(int), 1, f);
	//fwrite(&ALPHA_MODE, sizeof(int), 1, f);
	//fwrite(&GROWTH_TYPE, sizeof(int), 1, f);
	//fwrite(&MODE_TYPE, sizeof(int), 1, f);

	//for (int i = 0; i < GRID_SIZE; i++)
	//{
	//	fwrite(cells[i].grid, sizeof(fftwf_complex), MAX_CHANNELS, f);
	//	fwrite(cells[i].GA, sizeof(float), MAX_CHANNELS, f);
	//	fwrite(&cells[i].SUM, sizeof(float), 1, f);
	//	fwrite(cells[i].V_x, sizeof(float), MAX_CHANNELS, f);
	//	fwrite(cells[i].V_y, sizeof(float), MAX_CHANNELS, f);
	//	fwrite(cells[i].V_z, sizeof(float), MAX_CHANNELS, f);
	//}


	//for (int i = 0; i < channels; i++)
	//{
	//	fwrite(&ch[i].asymptotic, sizeof(bool), 1, f);
	//	fwrite(&ch[i].soft_clip, sizeof(bool), 1, f);
	//	fwrite(&ch[i].color, sizeof(D3DXCOLOR), 1, f);
	//}

	//fwrite(&kernel_count, sizeof(int), 1, f);
	//for (int i = 0; i < kernel_count; i++)
	//{
	//	fwrite(&kernels[i].radius, sizeof(int), 1, f);
	//	fwrite(&kernels[i].r_a, sizeof(float), 1, f);
	//	fwrite(&kernels[i].m, sizeof(float), 1, f);
	//	fwrite(&kernels[i].s, sizeof(float), 1, f);
	//	fwrite(&kernels[i].c0, sizeof(int), 1, f);
	//	fwrite(&kernels[i].c1, sizeof(int), 1, f);
	//	fwrite(&kernels[i].h, sizeof(float), 1, f);
	//	fwrite(&kernels[i].beta_size, sizeof(int), 1, f);
	//	fwrite(kernels[i].b, sizeof(float), kernels[i].beta_size, f);
	//}
	//fclose(f);

}

void load_config(Button * button)
{
	//WaitForSingleObject(sim_done_event, INFINITE);
	//WaitForSingleObject(mesh_done_event, INFINITE);
	//
	////WaitForMultipleObjects(MAX_THREADS, threads, true, INFINITE);

	////for (int i = 0; i < MAX_THREADS; i++)
	////{
	////	CloseHandle(threads[i]);
	////	cell_verts[i][cur_buf]->Unlock();
	////}

	////halted = true;
	////

	//char path[128];
	//char name[64];
	//GetWindowTextA(button->box_id, name, 122); 
	//const char * bi = ".bin";
	//path[0] = '.';
	//path[1] = '/';
	//path[2] = '\0';
	//strncat(path, name, 127);

	//strncat(path, bi, 127);

	//OutputDebugStringA(path);

	//FILE * f = NULL;
	//errno_t d = fopen_s(&f, path, "rb");
	//if (!f) return;
	//OutputDebugStringA(path);

	//fread(&ISO_VALUE, sizeof(float), 1, f);
	//fread(&timestep, sizeof(float), 1, f);
	//fread(&channels, sizeof(int), 1, f);
	//fread(&GTA, sizeof(float), 1, f);
	//fread(&STA, sizeof(float), 1, f);
	//fread(&V_BETA, sizeof(float), 1, f);
	//fread(&V_DAMP, sizeof(float), 1, f);
	//fread(&neighborhood, sizeof(float), 1, f);
	//fread(&N, sizeof(float), 1, f);
	//fread(&THETA_A, sizeof(float), 1, f);
	//fread(&ADV_VEL, sizeof(bool), 1, f);

	//fread(&OSC_LA, sizeof(bool), 1, f);
	//fread(&OSC_LA_AMT, sizeof(float), 1, f);

	//fread(&OSC_FLOW, sizeof(bool), 1, f);
	//fread(&OSC_FLOW_AMT, sizeof(float), 1, f);

	//fread(&OSC_DIV_FLOW, sizeof(bool), 1, f);
	//fread(&OSC_DIV_FLOW_AMT, sizeof(float), 1, f);

	//fread(&OSC_VEL, sizeof(bool), 1, f);
	//fread(&ADV_VEL, sizeof(bool), 1, f);
	//fread(&MU_VEL, sizeof(bool), 1, f);

	//fread(&BORDER, sizeof(int), 1, f);
	//fread(&SOBEL_BORDER, sizeof(int), 1, f);
	//fread(&ALPHA_MODE, sizeof(int), 1, f);
	//fread(&GROWTH_TYPE, sizeof(int), 1, f);
	//fread(&MODE_TYPE, sizeof(int), 1, f);

	//for (int i = 0; i < GRID_SIZE; i++)
	//{
	//	fread(cells[i].grid, sizeof(fftwf_complex), MAX_CHANNELS, f);
	//	fread(cells[i].GA, sizeof(float), MAX_CHANNELS, f);
	//	fread(&cells[i].SUM, sizeof(float), 1, f);
	//	fread(cells[i].V_x, sizeof(float), MAX_CHANNELS, f);
	//	fread(cells[i].V_y, sizeof(float), MAX_CHANNELS, f);
	//	fread(cells[i].V_z, sizeof(float), MAX_CHANNELS, f);
	//}

	//for (int i = 0; i < GRID_SIZE; i++)
	//{
	//	for (int c = 0; c < MAX_CHANNELS; c++)
	//	{
	//		cell_fft[c].val[i][0] = cells[i].grid[c]; 
	//		cell_fft[c].val[i][1] = 0;
	//	}
	//}

	//for (int i = 0; i < channels; i++)
	//{
	//	fread(&ch[i].asymptotic, sizeof(bool), 1, f);
	//	fread(&ch[i].soft_clip, sizeof(bool), 1, f);
	//	fread(&ch[i].color, sizeof(D3DXCOLOR), 1, f);
	//}

	//fwrite(&kernel_count, sizeof(int), 1, f);
	//for (int i = 0; i < kernel_count; i++)
	//{
	//	fread(&kernels[i].radius, sizeof(int), 1, f);
	//	fread(&kernels[i].r_a, sizeof(float), 1, f);
	//	fread(&kernels[i].m, sizeof(float), 1, f);
	//	fread(&kernels[i].s, sizeof(float), 1, f);
	//	fread(&kernels[i].c0, sizeof(int), 1, f);
	//	fread(&kernels[i].c1, sizeof(int), 1, f);
	//	fread(&kernels[i].h, sizeof(float), 1, f);
	//	fread(&kernels[i].beta_size, sizeof(int), 1, f);
	//	fread(kernels[i].b, sizeof(float), kernels[i].beta_size, f);
	//	init_kernel(&kernels[i]);
	//}

	//redraw_config = true;

	//fclose(f);
}


HWND make_edit(int width, int height, int x_offset, int y_offset, HINSTANCE hinst, HWND hwnd, int ci)
{
	HWND hwndEdit = CreateWindowA( 
		"EDIT",
		"",
		WS_TABSTOP | WS_VISIBLE | WS_CHILD | ES_LEFT, 
		x_offset,
		y_offset,
		width,
		height,
		hwnd,
		(HMENU)ci,
		(HINSTANCE)GetWindowLong(hwnd, GWL_HINSTANCE), 
    NULL);

	components[component_cnt++] = hwndEdit;
	SendMessage(hwndEdit, WM_SETFONT, (WPARAM)font, TRUE);

	return hwndEdit;

}


HWND make_button(int width, int height, int x_offset, int y_offset, LPSTR text, HINSTANCE hinst, HWND hwnd, int ci, HWND ti, void (*fn)(Button*), void (*action_fn)(Button *))
{
	HWND hwndButton = CreateWindowA( 
		"BUTTON",
		text,
		WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON, 
		x_offset,
		y_offset,
		width,
		height,
		hwnd,
		(HMENU)ci,
		(HINSTANCE)GetWindowLong(hwnd, GWL_HINSTANCE), 
    NULL);

	components[component_cnt++] = hwndButton;

	buttons[button_cnt].box_id = ti;
	buttons[button_cnt].cid = ci;
	buttons[button_cnt].fn = fn;
	buttons[button_cnt].action_fn = action_fn;

	button_cnt++;

	SendMessage(hwndButton, WM_SETFONT, (WPARAM)font, TRUE);
	return hwndButton;
}

HWND make_text_box(int width, int height, int x_offset, int y_offset, DWORD style, LPCSTR text, HINSTANCE hinst, HWND hwnd, int ci)
{

	HWND tbox = CreateWindowA( 
		"STATIC",
		text,
		WS_TABSTOP | WS_VISIBLE | WS_CHILD | style, 
		x_offset,
		y_offset,
		width,
		height,
		hwnd,
		HMENU(ci),
		(HINSTANCE)GetWindowLong(hwnd, GWL_HINSTANCE), 
    NULL);
	boxes[box_cnt++].cid = ci;
	components[component_cnt++] = tbox;
	SendMessage(tbox, WM_SETFONT, (WPARAM)font, TRUE);

	return tbox;
}


void increase_scale(Button * button)
{
	scale += 0.005f;
	char ts[24];
	sprintf_s(ts, 23, "%f", scale);
	SetWindowTextA(button->box_id, ts);
}
void decrease_scale(Button * button)
{
	scale -= 0.005f;
	char ts[24];
	sprintf_s(ts, 23, "%f", scale);
	SetWindowTextA(button->box_id, ts);
}

void color_dialog_init(Button * button){


	CHOOSECOLOR cc;
	::ZeroMemory(&cc, sizeof(CHOOSECOLOR));
	COLORREF cust[16];
	cc.hInstance = cfg_wnd;
	cc.Flags = CC_ANYCOLOR;
	cc.hwndOwner = cfg_wnd;
	cc.lpCustColors = cust;
	cc.lStructSize = sizeof(CHOOSECOLOR);
	bool res = ChooseColor(&cc);

	if (res)
	{
		(*(D3DXCOLOR*)button->data) = D3DXCOLOR(
			GetRValue(cc.rgbResult), GetGValue(cc.rgbResult), GetBValue(cc.rgbResult), 0.5);
	}
}


void remove_channel(Button * button)
{

	for (int i = button->ex_int; i < MAX_CHANNELS - 1; i++)
	{
		ch[i] = ch[i+1];
	}
	channels--;

}


void remove_kernel(Button * button)
{
	destroy_kernel(&kernels[button->ex_int]);

	for (int i = button->ex_int; i < MAX_KERNELS - 1; i++)
	{
		kernels[i] = kernels[i+1];
	}
	kernel_count--;

}

void recompute_kernel(Button * button)
{
	init_kernel(&kernels[button->ex_int]);

}



void increment(Button * button)
{

	switch(button->dtype)
	{
		case DT_INT:
			{
				int * d = (int*)button->data;
				(*d) += 1;
				char ts[24];
				sprintf_s(ts, 23, "%d", *d);
				SetWindowTextA(button->box_id, ts);
			}
			break;
		case DT_FLOAT:
			{
				float * d = (float*)button->data;
				(*d) += scale;
				char ts[24];
				sprintf_s(ts, 23, "%.3f", *d);
				SetWindowTextA(button->box_id, ts);
			}
			
			break;
		case DT_BOOL:
		{
	
			bool * d = (bool*)button->data;
			(*d) = !(*d);
			char ts[24];
			sprintf_s(ts, 23, "%s", BOOL_VAL[*d]);
			SetWindowTextA(button->box_id, ts);

		}
		break;
		case DT_ENUM:
			{
		
				*(int*)button->data = (*(int*)button->data + 1) % button->ex_int;
				
				char ts[24];
				sprintf_s(ts, 23, "%s", ((const char (**)[])button->action_data)[*(int*)button->data]);
				SetWindowTextA(button->box_id, ts);

			}
			break;
		default:
			break;
	}
	if (button->redraw) redraw_config = true;

}
void decrement(Button * button)
{
	switch(button->dtype)
	{
		case DT_INT:
			{
				int * d = (int*)button->data;
				(*d) -= 1;
				char ts[24];
				sprintf_s(ts, 23, "%d", *d);
				SetWindowTextA(button->box_id, ts);
			}			
			break;
		case DT_FLOAT:
			{
				
				float * d = (float*)button->data;
				(*d) -= scale;
				char ts[24];
				sprintf_s(ts, 23, "%.3f", *d);
				SetWindowTextA(button->box_id, ts);
			}
			break;
		case DT_BOOL:
			{
		
				bool * d = (bool*)button->data;
				(*d) = !(*d);
				char ts[24];
				sprintf_s(ts, 23, "%s", BOOL_VAL[*d]);
				SetWindowTextA(button->box_id, ts);

			}
			break;
		case DT_ENUM:
			{
		
				*(int*)button->data = (*(int*)button->data + 1) % button->ex_int;

				char ts[24];
				sprintf_s(ts, 23, "%s", (*(const char**)button->action_data)[*(int*)button->data]);
				SetWindowTextA(button->box_id, ts);

			}
			break;
		default:
			break;
	}

	if (button->redraw) redraw_config = true;

}



//
// kernel data
// r
// r_a
// m
// s
// beta_size
// beta values
// c0
// c1


bool draw_config()
{


	UpdateWindow(cfg_wnd);
	HDC hdc = GetDC(cfg_wnd);

	if (hdc)
	{
		RECT c;
		GetClientRect(cfg_wnd, &c);
		HBRUSH h = CreateSolidBrush(RGB(255, 255, 255));
		
		if (h)
		{
			FillRect(hdc, &c, h);
			DeleteObject(h);
		}
		ReleaseDC(cfg_wnd, hdc);
	}

	for (int i = 0; i < component_cnt; i++)
	{
		PostMessage(components[i], WM_CLOSE, 0, 0);
	}

	memset(buttons, 0, sizeof(Button) * MAX_BUTTONS);
	memset(boxes, 0, sizeof(TextBox) * MAX_TEXT_BOX);
	memset(components, 0, sizeof(HWND) * MAX_COMPONENTS);

	box_cnt = 0;
	component_cnt = 0;
	button_cnt = 0;
	u_cid = 0;
	int offset = 3;
	int row = offset;
	int border = 5;

	char ts[24];

	int height = 15;


	// scale
	// timestep
	// field size
	// iso value
	// thread count
	int llabel = 50;
	int bt = height;
	int mlabel = 50;
	int off = border;
	SendMessage(cfg_wnd, WM_SETREDRAW, false, 0);
	sprintf_s(ts, 23, "Scale:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off = border + llabel + border + bt + border;
	sprintf_s(ts, 23, "%.3f", scale);
	HWND t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = border + llabel + border;
	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrease_scale, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&scale;

	off = border + llabel + border + bt + border + mlabel + border;
	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increase_scale, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&scale;

	
	off += bt + border;
	int placeholder = off;
	sprintf_s(ts, 23, "Timestep:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + llabel + border + bt + border;

	sprintf_s(ts, 23, "%.3f", timestep);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder +  llabel + border;

	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&timestep;
	
	off = placeholder + llabel + border + bt + border + mlabel + border;

	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&timestep;

	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "Iso Value:");
	make_text_box(llabel, height, border, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);
	
	off = border + llabel + border + bt + border;
	sprintf_s(ts, 23, "%.3f", ISO_VALUE);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = border + llabel + border;
	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&ISO_VALUE;
	
	off = border + llabel + border + bt + border + mlabel + border;
	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&ISO_VALUE;

	off += bt + border;
	placeholder = off;

	sprintf_s(ts, 23, "Channels:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + border + llabel + border + bt + border;
	sprintf_s(ts, 23, "%d", channels);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + border + llabel + border;
	make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_INT;
	buttons[button_cnt-1].data = (void*)&channels;
	buttons[button_cnt-1].redraw = true;

	off = placeholder + border + llabel + border + bt + border + mlabel + border;
	make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_INT;
	buttons[button_cnt-1].data = (void*)&channels;
	buttons[button_cnt-1].redraw = true;

	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "Mode:");
	make_text_box(40, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += 40 + border;
	sprintf_s(ts, 23, "%s", MODE_TYPES[MODE_TYPE]);

	t0 = make_button(75, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_ENUM;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&MODE_TYPE;
	buttons[button_cnt-1].ex_int = MODE_MAX;
	buttons[button_cnt-1].action_data = (void*)&MODE_TYPES;

	off += 75 + border;

	sprintf_s(ts, 23, "Alpha:");
	make_text_box(35, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += 30 + border;
	sprintf_s(ts, 23, "%s", ALPHA_TYPES[ALPHA_MODE]);

	t0 = make_button(135, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_ENUM;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&ALPHA_MODE;
	buttons[button_cnt-1].ex_int = ALPHA_MAX;
	buttons[button_cnt-1].action_data = (void*)&ALPHA_TYPES;

	row += height + offset;

	off = border;

	sprintf_s(ts, 23, "Growth:");
	make_text_box(40, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += 40 + border;
	sprintf_s(ts, 23, "%s", GROWTH_TYPES[GROWTH_TYPE]);

	t0 = make_button(150, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_ENUM;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&GROWTH_TYPE;
	buttons[button_cnt-1].ex_int = GROWTH_MAX;
	buttons[button_cnt-1].action_data = (void*)&GROWTH_TYPES;

	off += 150 + border;

	placeholder = off;
	sprintf_s(ts, 23, "S. Norm:");
	make_text_box(50, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + 50 + border + bt + border;

	sprintf_s(ts, 23, "%.3f", S_NORM);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder +  50 + border;

	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&S_NORM;
	
	off = placeholder + 50 + border + bt + border + mlabel + border;

	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&S_NORM;

	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "Theta A:");
	make_text_box(llabel, height, border, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);
	
	off = border + llabel + border + bt + border;
	sprintf_s(ts, 23, "%.3f", THETA_A);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = border + llabel + border;
	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&THETA_A;
	
	off = border + llabel + border + bt + border + mlabel + border;
	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&THETA_A;

	off += bt + border;
	placeholder = off;
	sprintf_s(ts, 23, "N:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + llabel + border + bt + border;

	sprintf_s(ts, 23, "%.3f", N);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder +  llabel + border;

	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&N;
	
	off = placeholder + llabel + border + bt + border + mlabel + border;

	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&N;

	row += height + offset;


	int odf = 80;
	off = border;
	sprintf_s(ts, 23, "Osc LA:");
	make_text_box(odf, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += odf + border;
	sprintf_s(ts, 23, "%s", BOOL_VAL[OSC_LA]);

	t0 = make_button(mlabel, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_BOOL;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&OSC_LA;

	off += mlabel + border;
	placeholder = off;
	off = placeholder + bt + border;

	sprintf_s(ts, 23, "%.3f", OSC_LA_AMT);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + border;

	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&OSC_LA_AMT;
	
	off = placeholder + bt + border + mlabel + border;

	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&OSC_LA_AMT;

	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "Osc Flow:");
	make_text_box(odf, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += odf + border;
	sprintf_s(ts, 23, "%s", BOOL_VAL[OSC_FLOW]);

	t0 = make_button(mlabel, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_BOOL;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&OSC_FLOW;

	off += mlabel + border;
	placeholder = off;
	off = placeholder + bt + border;

	sprintf_s(ts, 23, "%.3f", OSC_FLOW_AMT);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + border;

	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&OSC_FLOW_AMT;
	
	off = placeholder + bt + border + mlabel + border;

	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&OSC_FLOW_AMT;

	row += height + offset;

	off = border;
	
	sprintf_s(ts, 23, "Osc Div Flow:");
	make_text_box(odf, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += odf + border;
	sprintf_s(ts, 23, "%s", BOOL_VAL[OSC_DIV_FLOW]);

	t0 = make_button(mlabel, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_BOOL;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&OSC_DIV_FLOW;

	off += mlabel + border;
	placeholder = off;
	off = placeholder + bt + border;

	sprintf_s(ts, 23, "%.3f", OSC_DIV_FLOW_AMT);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + border;

	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&OSC_DIV_FLOW_AMT;
	
	off = placeholder + bt + border + mlabel + border;

	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&OSC_DIV_FLOW_AMT;

	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "Adv. Vel:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += llabel + border;
	sprintf_s(ts, 23, "%s", BOOL_VAL[ADV_VEL]);

	t0 = make_button(50, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_BOOL;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&ADV_VEL;


	off += 50 + border;
	placeholder = off;

	sprintf_s(ts, 23, "Mu Vel:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += llabel + border;
	sprintf_s(ts, 23, "%s", BOOL_VAL[MU_VEL]);

	t0 = make_button(50, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_BOOL;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&MU_VEL;
	off += 50 + border;

	sprintf_s(ts, 23, "Osc Vel:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += llabel + border;
	sprintf_s(ts, 23, "%s", BOOL_VAL[OSC_VEL]);

	t0 = make_button(50, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_BOOL;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&OSC_VEL;

	row += height + offset;


	off = border;
	sprintf_s(ts, 23, "Sigma:");
	make_text_box(llabel, height, border, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);
	
	off = border + llabel + border + bt + border;
	sprintf_s(ts, 23, "%.3f", off_table.sigma);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = border + llabel + border;
	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, init_offset_table);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&off_table.sigma;
	
	off = border + llabel + border + bt + border + mlabel + border;
	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, init_offset_table);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&off_table.sigma;

	off += bt + border;
	placeholder = off;
	sprintf_s(ts, 23, "DD:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + llabel + border + bt + border;

	sprintf_s(ts, 23, "%.3f", off_table.dd);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder +  llabel + border;

	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, init_offset_table);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&off_table.dd;
	
	off = placeholder + llabel + border + bt + border + mlabel + border;

	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, init_offset_table);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&off_table.dd;


	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "GTA:");
	make_text_box(llabel, height, border, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);
	
	off = border + llabel + border + bt + border;
	sprintf_s(ts, 23, "%.3f", GTA);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = border + llabel + border;
	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&GTA;
	
	off = border + llabel + border + bt + border + mlabel + border;
	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&GTA;

	off += bt + border;
	placeholder = off;
	sprintf_s(ts, 23, "STA:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + llabel + border + bt + border;

	sprintf_s(ts, 23, "%.3f", STA);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder +  llabel + border;

	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&STA;
	
	off = placeholder + llabel + border + bt + border + mlabel + border;

	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&STA;
	
	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "VB:");
	make_text_box(llabel, height, border, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);
	
	off = border + llabel + border + bt + border;
	sprintf_s(ts, 23, "%.3f", V_BETA);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = border + llabel + border;
	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&V_BETA;
	
	off = border + llabel + border + bt + border + mlabel + border;
	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&V_BETA;

	off += bt + border;
	placeholder = off;
	sprintf_s(ts, 23, "V Damp:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder + llabel + border + bt + border;

	sprintf_s(ts, 23, "%.3f", V_DAMP);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = placeholder +  llabel + border;

	make_button(bt, height, off, row, "<", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&V_DAMP;
	
	off = placeholder + llabel + border + bt + border + mlabel + border;

	make_button(bt, height, off, row, ">", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_FLOAT;
	buttons[button_cnt-1].data = (void*)&V_DAMP;
	
	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "Border:");
	make_text_box(40, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += 40 + border;
	sprintf_s(ts, 23, "%s", BORDER_TYPES[BORDER]);

	t0 = make_button(100, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_ENUM;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&BORDER;
	buttons[button_cnt-1].ex_int = BORDER_MAX;
	buttons[button_cnt-1].action_data = (void*)&BORDER_TYPES;

	off += 100 + border;

	sprintf_s(ts, 23, "Sobel Border:");
	make_text_box(40, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += 40 + border;
	sprintf_s(ts, 23, "%s", BORDER_TYPES[SOBEL_BORDER]);

	t0 = make_button(100, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_ENUM;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&SOBEL_BORDER;
	buttons[button_cnt-1].ex_int = BORDER_MAX;
	buttons[button_cnt-1].action_data = (void*)&BORDER_TYPES;

	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "Views:");
	make_text_box(40, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += 40 + border;
	sprintf_s(ts, 23, "%s", VIEW_TYPES[VIEWS]);

	t0 = make_button(100, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_ENUM;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&VIEWS;
	buttons[button_cnt-1].ex_int = VIEW_MAX;
	buttons[button_cnt-1].action_data = (void*)&VIEW_TYPES;

	off += 100 + border;

	sprintf_s(ts, 23, "Blur CD:");
	make_text_box(40, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += 40 + border;
	sprintf_s(ts, 23, "%s", BOOL_VAL[BLUR_CENTRAL_DIFF]);

	t0 = make_button(100, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_BOOL;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&BLUR_CENTRAL_DIFF;
	//buttons[button_cnt-1].ex_int = GRADIENT_MAX;
	buttons[button_cnt-1].action_data = (void*)&BOOL_VAL;

	row += height + offset;

	off = border;
	sprintf_s(ts, 23, "Gradient:");
	make_text_box(40, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

	off += 40 + border;
	sprintf_s(ts, 23, "%s", GRADIENT_TYPES[GRADIENT_MODE]);

	t0 = make_button(160, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, increment, NULL);
	buttons[button_cnt-1].dtype = DT_ENUM;
	buttons[button_cnt-1].box_id = t0;
	buttons[button_cnt-1].data = (void*)&GRADIENT_MODE;
	buttons[button_cnt-1].ex_int = GRADIENT_MAX;
	buttons[button_cnt-1].action_data = (void*)&GRADIENT_TYPES;

	row += height + offset;

	t0 = make_edit(80, height, border, row, hinst, cfg_wnd, u_cid++);

	make_button(mlabel, 20, border + 80 + border, row, "Save", hinst, cfg_wnd, u_cid++, t0, save_config, NULL);
	buttons[button_cnt-1].data = NULL;
	row += height + offset;

	t0 = make_edit(80, height, border, row, hinst, cfg_wnd, u_cid++);

	make_button(mlabel, 20, border + 80 + border, row, "Load", hinst, cfg_wnd, u_cid++, t0, load_config, NULL);
	buttons[button_cnt-1].data = NULL;
	row += height + offset;





	// channel data
	// asymptotic 
	// soft clip
	// color
	for (int i = 0; i < channels; i++)
	{

		off = border;
		llabel = 60;
		mlabel = 40;
		
		sprintf_s(ts, 23, "Channel %d:", i);
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off = border + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, remove_channel);
		buttons[button_cnt-1].dtype = DT_NULL;
		buttons[button_cnt-1].data = NULL;
		buttons[button_cnt-1].redraw = true;
		buttons[button_cnt-1].ex_int = i;
		
		row += height + offset;

		sprintf_s(ts, 23, "Asymptotic:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off += llabel + border;
		sprintf_s(ts, 23, "%s", BOOL_VAL[ch[i].asymptotic]);

		t0 = make_button(mlabel, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
		buttons[button_cnt-1].dtype = DT_BOOL;
		buttons[button_cnt-1].box_id = t0;
		buttons[button_cnt-1].data = (void*)&ch[i].asymptotic;

		off += mlabel + border;
		llabel = 50;
		sprintf_s(ts, 23, "Soft Clip:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off += llabel + border;
		sprintf_s(ts, 23, "%s", BOOL_VAL[ch[i].soft_clip]);

		t0 = make_button(mlabel, 20, off, row - 3, ts, hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
		buttons[button_cnt-1].dtype = DT_BOOL;
		buttons[button_cnt-1].box_id = t0;
		buttons[button_cnt-1].data = (void*)&ch[i].soft_clip;

		off += mlabel + border;

		make_button(mlabel, 20, off, row - 3, "Color", hinst, cfg_wnd, u_cid++, t0, color_dialog_init, NULL);
		buttons[button_cnt-1].dtype = DT_D3DXCOLOR;
		buttons[button_cnt-1].data = (void*)&ch[i].color;

		row += height + offset;
	}

	// kernel data
	// r
	// r_a
	// m
	// s
	// beta_size
	// beta values
	// c0
	// c1

	off = border;
	sprintf_s(ts, 23, "Kernels:");
	make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);
	
	off = border + llabel + border + bt + border;
	sprintf_s(ts, 23, "%d", kernel_count);
	t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

	off = border + llabel + border + bt + border + mlabel + border;
	make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
	buttons[button_cnt-1].dtype = DT_INT;
	buttons[button_cnt-1].data = (void*)&kernel_count;
	buttons[button_cnt-1].ex_int = kernel_count;
	buttons[button_cnt-1].redraw = true;


	row += height + offset;

	for (int i = 0; i < kernel_count; i++)
	{

		off = border;
		
		sprintf_s(ts, 23, "Kernel %d:", i);
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off = border + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, remove_kernel);
		buttons[button_cnt-1].dtype = DT_NULL;
		buttons[button_cnt-1].data = NULL;
		buttons[button_cnt-1].redraw = true;
		buttons[button_cnt-1].ex_int = i;
		

		row += height + offset;
		off = border;
		sprintf_s(ts, 23, "Radius:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off = border + llabel + border + bt + border;
		sprintf_s(ts, 23, "%d", kernels[i].radius);
		t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

		off = border + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_INT;
		buttons[button_cnt-1].data = (void*)&kernels[i].radius;
		buttons[button_cnt-1].ex_int = i;

		off = border + llabel + border + bt + border + mlabel + border;
		make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_INT;
		buttons[button_cnt-1].data = (void*)&kernels[i].radius;
		buttons[button_cnt-1].ex_int = i;

		off += bt + border;
		placeholder = off;
		sprintf_s(ts, 23, "R_A:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off = placeholder + llabel + border + bt + border;
		sprintf_s(ts, 23, "%.3f", kernels[i].r_a);
		t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

		off = placeholder + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_FLOAT;
		buttons[button_cnt-1].data = (void*)&kernels[i].r_a;
		buttons[button_cnt-1].ex_int = i;

		off = placeholder + llabel + border + bt + border + mlabel + border;
		make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_FLOAT;
		buttons[button_cnt-1].data = (void*)&kernels[i].r_a;
		buttons[button_cnt-1].ex_int = i;
		
		row += height + offset;

		off = border;
		sprintf_s(ts, 23, "m:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);


		off = border + llabel + border + bt + border;
		sprintf_s(ts, 23, "%.3f", kernels[i].m);
		t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

		off = border + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_FLOAT;
		buttons[button_cnt-1].data = (void*)&kernels[i].m;
		buttons[button_cnt-1].ex_int = i;

		off = border + llabel + border + mlabel + border + bt + border;
		make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_FLOAT;
		buttons[button_cnt-1].data = (void*)&kernels[i].m;
		buttons[button_cnt-1].ex_int = i;
		

		off += bt + border;
		placeholder = off;
		sprintf_s(ts, 23, "s:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off = placeholder + llabel + border + bt + border;

		sprintf_s(ts, 23, "%.3f", kernels[i].s);
		t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

		off = placeholder + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_FLOAT;
		buttons[button_cnt-1].data = (void*)&kernels[i].s;
		buttons[button_cnt-1].ex_int = i;

		off = placeholder + llabel + border + bt + border + mlabel + border;
		make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_FLOAT;
		buttons[button_cnt-1].data = (void*)&kernels[i].s;
		buttons[button_cnt-1].ex_int = i;
		
		row += height + offset;

		off = border;
		sprintf_s(ts, 23, "c0:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);


		off = border + llabel + border + bt + border;
		sprintf_s(ts, 23, "%d", kernels[i].c0);
		t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

		off = border + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_INT;
		buttons[button_cnt-1].data = (void*)&kernels[i].c0;
		buttons[button_cnt-1].ex_int = i;

		off = border + llabel + border + mlabel + border + bt + border;
		make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_INT;
		buttons[button_cnt-1].data = (void*)&kernels[i].c0;
		buttons[button_cnt-1].ex_int = i;
		

		off += bt + border;
		placeholder = off;
		sprintf_s(ts, 23, "c1:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off = placeholder + llabel + border + bt + border;

		sprintf_s(ts, 23, "%d", kernels[i].c1);
		t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

		off = placeholder + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_INT;
		buttons[button_cnt-1].data = (void*)&kernels[i].c1;
		buttons[button_cnt-1].ex_int = i;

		off = placeholder + llabel + border + bt + border + mlabel + border;
		make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_INT;
		buttons[button_cnt-1].data = (void*)&kernels[i].c1;
		buttons[button_cnt-1].ex_int = i;
		
		row += height + offset;

		off = border;
		sprintf_s(ts, 23, "Beta Size:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off = border + llabel + border + bt + border;
		sprintf_s(ts, 23, "%d", kernels[i].beta_size);
		t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

		off = border + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_INT;
		buttons[button_cnt-1].data = (void*)&kernels[i].beta_size;
		buttons[button_cnt-1].redraw = true;
		buttons[button_cnt-1].ex_int = i;

		off = border + llabel + border + mlabel + border + bt + border;
		make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
		buttons[button_cnt-1].dtype = DT_INT;
		buttons[button_cnt-1].data = (void*)&kernels[i].beta_size;
		buttons[button_cnt-1].redraw = true;
		buttons[button_cnt-1].ex_int = i;

		off += bt + border;
		placeholder = off;
		sprintf_s(ts, 23, "H:");
		make_text_box(llabel, height, off, row, SS_LEFT, ts, hinst, cfg_wnd, u_cid++);

		off = placeholder + llabel + border + bt + border;
		sprintf_s(ts, 23, "%d", kernels[i].h);
		t0 = make_text_box(mlabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

		off = placeholder + llabel + border;
		make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, NULL);
		buttons[button_cnt-1].dtype = DT_FLOAT;
		buttons[button_cnt-1].data = (void*)&kernels[i].h;
		buttons[button_cnt-1].redraw = true;

		off = placeholder + llabel + border + mlabel + border + bt + border;
		make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, NULL);
		buttons[button_cnt-1].dtype = DT_FLOAT;
		buttons[button_cnt-1].data = (void*)&kernels[i].h;
		buttons[button_cnt-1].redraw = true;



		row += height + offset;

		for (int b = 0; b < kernels[i].beta_size; b+=3)
		{

			off = border + bt + border;
			sprintf_s(ts, 23, "%.3f", kernels[i].b[b]);
			t0 = make_text_box(llabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

			off = border;
			make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
			buttons[button_cnt-1].dtype = DT_FLOAT;
			buttons[button_cnt-1].data = (void*)&kernels[i].c0;
			buttons[button_cnt-1].ex_int = i;

			off = border + llabel + bt + border;
			make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
			buttons[button_cnt-1].dtype = DT_FLOAT;
			buttons[button_cnt-1].data = (void*)&kernels[i].c0;
			buttons[button_cnt-1].ex_int = i;
			

			if (b+1 >= kernels[i].beta_size) {
				row += height + offset;
				break;
			}

			off += bt + border;
			placeholder = off;

			off = placeholder + bt + border;

			sprintf_s(ts, 23, "%.3f", kernels[i].b[b+1]);
			t0 = make_text_box(llabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

			off = placeholder;
			make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
			buttons[button_cnt-1].dtype = DT_FLOAT;
			buttons[button_cnt-1].data = (void*)&kernels[i].b[b+1];
			buttons[button_cnt-1].ex_int = i;

			off = placeholder + bt + border + llabel + border;
			make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
			buttons[button_cnt-1].dtype = DT_FLOAT;
			buttons[button_cnt-1].data = (void*)&kernels[i].b[b+1];
			buttons[button_cnt-1].ex_int = i;
			
			if (b+2 >= kernels[i].beta_size) {
				row += height + offset;
				break;
			}

			off += bt + border;
			placeholder = off;

			off = placeholder + bt + border;

			sprintf_s(ts, 23, "%.3f", kernels[i].b[b+2]);
			t0 = make_text_box(llabel, height, off, row, SS_CENTER, ts, hinst, cfg_wnd, u_cid++);

			off = placeholder;
			make_button(bt, height, off, row, "-", hinst, cfg_wnd, u_cid++, t0, decrement, recompute_kernel);
			buttons[button_cnt-1].dtype = DT_FLOAT;
			buttons[button_cnt-1].data = (void*)&kernels[i].b[b+2];
			buttons[button_cnt-1].ex_int = i;

			off = placeholder + bt + border + llabel + border;
			make_button(bt, height, off, row, "+", hinst, cfg_wnd, u_cid++, t0, increment, recompute_kernel);
			buttons[button_cnt-1].dtype = DT_FLOAT;
			buttons[button_cnt-1].data = (void*)&kernels[i].b[b+2];
			buttons[button_cnt-1].ex_int = i;



			row += height + offset;

		}
		
	}

	
	

	SendMessage(cfg_wnd, WM_SETREDRAW, TRUE, 0);
	InvalidateRect(cfg_wnd, NULL, TRUE);

	//ScrollWindow(

	return true;
}

void dispatch_control(WPARAM wParam)
{
	halted = true;
	WaitForSingleObject(sim_done_event, INFINITE);
	// dont even ask
	while (!sim_done)
	{
		int a = 0;
	}
	//WaitForSingleObject(mesh_done_event, INFINITE);
	
	for (int i = 0; i < button_cnt; i++)
	{
		if (wParam == buttons[i].cid)
		{
			
			buttons[i].fn(&buttons[i]);
			if (buttons[i].action_fn)
			{
				buttons[i].action_fn(&buttons[i]);
			}
			break;
		}
	}
	//if (!sim_done)
	//{
	//	int i = 1;
	//}
	halted = false;

}
int g_scrollY = 0;

LRESULT CALLBACK cfg_wnd_proc(HWND hwnd, UINT msg, WPARAM wParam,
	LPARAM lParam)
{
	switch( msg )
	{
		case WM_DESTROY:
			::PostQuitMessage(0);
			break;
		case WM_COMMAND:
			{
				dispatch_control(wParam);
			}
			break;
		case WM_VSCROLL:
			{
				// khayal said to put a ? here. explain it to me.
	
				WORD action = LOWORD(wParam);
				int pos = -1;
				if (action == SB_THUMBPOSITION || action == SB_THUMBTRACK)
				{
					pos = HIWORD(wParam);
				}

				if (pos == -1) break;
				SCROLLINFO si = {0};
				si.cbSize = sizeof(SCROLLINFO);
				si.fMask = SIF_POS;
				si.nPos = pos;
				si.nTrackPos = 0;
				
				SetScrollInfo(cfg_wnd, SB_VERT, &si, true);
				GetScrollInfo(cfg_wnd, SB_VERT, &si);

				pos = si.nPos;
				POINT pt;
				pt.x = 0;
				pt.y = pos - g_scrollY;
				
				HDC hdc = GetDC(cfg_wnd);
				LPtoDP(hdc, &pt, 1);
				ReleaseDC(cfg_wnd, hdc);
				ScrollWindow(cfg_wnd, 0, -pt.y, NULL, NULL);
				g_scrollY = pos;
				return 0;






			}
			break;
	}
	return ::DefWindowProc(hwnd, msg, wParam, lParam);
}

bool config_window_init(HINSTANCE hinst,
				  int width, int height,
				  bool windowed)
{





	// init window
	WNDCLASS wc;
	wc.style = CS_HREDRAW | CS_VREDRAW;
	wc.lpfnWndProc = cfg_wnd_proc;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hInstance = hinst;
	wc.hIcon = ::LoadIcon(0, IDI_APPLICATION);
	wc.hCursor = ::LoadCursor(0, IDC_ARROW);
	wc.hbrBackground =
	static_cast<HBRUSH>(::GetStockObject(WHITE_BRUSH));
	wc.lpszMenuName = 0;
	wc.lpszClassName = TEXT("Config");
	if(!::RegisterClass(&wc))
	{
		::MessageBox(0, TEXT("RegisterClass - Failed"), 0, 0);
		return false;
	}
	cfg_wnd =  ::CreateWindow(
				TEXT("Config"),
				TEXT("Settings"),
				WS_OVERLAPPEDWINDOW | WS_VSCROLL,
				width,
				height,
				CW_USEDEFAULT,
				CW_USEDEFAULT,
				0,
				0,
				hinst,
				0);
	if (cfg_wnd == 0)
	{
		return false;
	}


	HDC hdc = GetDC(NULL);
	long lfHeight;
	lfHeight = -MulDiv(6, GetDeviceCaps(hdc, LOGPIXELSY), 72);
	ReleaseDC(NULL, hdc);
	font = CreateFont(
	   lfHeight,                        // nHeight
	   0,                         // nWidth
	   0,                         // nEscapement
	   0,                         // nOrientation
	   0,                 // nWeight
	   FALSE,                     // bItalic
	   FALSE,                     // bUnderline
	   0,                         // cStrikeOut
	   ANSI_CHARSET,              // nCharSet
	   OUT_DEFAULT_PRECIS,        // nOutPrecision
	   CLIP_DEFAULT_PRECIS,       // nClipPrecision
	   DEFAULT_QUALITY,           // nQuality
	   DEFAULT_PITCH | FF_SWISS,  // nPitchAndFamily
	   TEXT("Times New Roman"));                 // lpszFacename

	SendMessage(cfg_wnd, WM_SETFONT, (WPARAM)font, TRUE);



	SetScrollRange(cfg_wnd, SB_VERT, 0, 500, false);
	
	::ShowWindow(cfg_wnd, true);
	::UpdateWindow(cfg_wnd);
	

	//make_button(100, 100, 10, 10, hinst, cfg_wnd);
	//make_text_box(100, 100, 10, 10, hinst, cfg_wnd);

	return true;
}
