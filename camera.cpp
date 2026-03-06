#include "camera.h"
#include "render.h"
#define _USE_MATH_DEFINES
#include <math.h>

D3DXVECTOR3 position(2.0f, 1.0f, -5.0f);
D3DXVECTOR3 target(0.0f, 0.0f, 0.0f);
D3DXVECTOR3 up(0.0f, 1.0f, 0.0f);
D3DXVECTOR3 forward(0.0f, 0.0f, 0.0f);
D3DXVECTOR3 right(0.0f, 0.0f, 0.0f);
float yaw = 0.0f;
float pitch = 0.0f;

D3DXMATRIX viewproj;
D3DXMATRIX view;

void move(vec3 v){

	float x = v.x * forward.x + v.z * right.x;
	float z = v.x * forward.z + v.z * right.z;

	position.x += x;
	position.y += v.y;
	position.z += z;
}

void update_camera(){

	forward = D3DXVECTOR3(
		cos(yaw) * cos(pitch), 
		sin(pitch), 
		sin(yaw) * cos(pitch));
	right = D3DXVECTOR3(-sin(yaw), 0.0f, cos(yaw));

	target = forward + position;
 

	D3DXMatrixLookAtLH(&view, &position, &target, &up);
	Device->SetTransform(D3DTS_VIEW, &view);
	D3DXMATRIX proj;
	D3DXMatrixPerspectiveFovLH(
		&proj,
		D3DX_PI * 0.5f,
		(float)Width / (float)Height,
		1.0f,
		1000.0f);
	Device->SetTransform(D3DTS_PROJECTION, &proj);
	viewproj = view * proj;
}

