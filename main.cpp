//--------------------------------------------------------------------------------------
// File: cpp
//
// Empty starting point for new Direct3D 9 and/or Direct3D 11 applications
//
// Copyright (c) Microsoft Corporation. All rights reserved.
//--------------------------------------------------------------------------------------
#include "DXUT.h"
#include "Atmosphere.h"
#include "Model/Model.h"
#include "FrameBuffer.h"

#ifdef _DEBUG 
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <malloc.h>    // 解决 malloc.h 与 crtdbg.h 顺序导致的 Debug Assertion Failed, "Corrupted pointer passed to _freea" 。
#include <crtdbg.h>
#define new new( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#endif  // _DEBUG

Atmosphere*							m_pAtmosphere;
Model*								m_pModel;
FrameBuffer*						m_pFrameBuffer;
FrameBuffer*						m_pShadowMapFrameBuffer;
int									screen_width = 1920;
int									screen_height = 1080;
int									m_shadowMapDim = 1024;
D3D11_VIEWPORT						g_viewport;
D3DXMATRIX							g_World;
D3DXMATRIX							g_View;
D3DXMATRIX							g_Projection;
float								m_EyeHeight = 0.f;	// Unit:km 
float								m_ModelScaling = 1;
D3DXVECTOR3							g_Eye(10 * m_ModelScaling, m_EyeHeight, 20 * m_ModelScaling);
D3DXVECTOR3							g_At(0.0f, m_EyeHeight, 20 * m_ModelScaling);
D3DXVECTOR3							g_CamDir(-10 * m_ModelScaling, 0, 0);
DirectionalLight					g_DirectionalLight;
CFirstPersonCamera					m_Camera;
ID3D11Device*						g_pd3dDevice = nullptr;
ID3D11DeviceContext*				g_pd3dImmediateContext = nullptr;
ID3D11RenderTargetView*				g_pRenderTargetView = nullptr;
ID3D11DepthStencilView*				g_pDepthStencilView = nullptr;

void UpdateMatrix()
{
	D3DXMatrixIdentity(&g_World);
	g_View = *(m_Camera.GetViewMatrix());
	g_Projection = *(m_Camera.GetProjMatrix());
	g_Eye = *(m_Camera.GetEyePt());
	g_At = *(m_Camera.GetLookAtPt());
	g_CamDir = g_At - g_Eye;
}

//--------------------------------------------------------------------------------------
// Reject any D3D11 devices that aren't acceptable by returning false
//--------------------------------------------------------------------------------------
bool CALLBACK IsD3D11DeviceAcceptable( const CD3D11EnumAdapterInfo *AdapterInfo, UINT Output, const CD3D11EnumDeviceInfo *DeviceInfo,
                                       DXGI_FORMAT BackBufferFormat, bool bWindowed, void* pUserContext )
{
    return true;
}


//--------------------------------------------------------------------------------------
// Called right before creating a D3D9 or D3D11 device, allowing the app to modify the device settings as needed
//--------------------------------------------------------------------------------------
bool CALLBACK ModifyDeviceSettings( DXUTDeviceSettings* pDeviceSettings, void* pUserContext )
{
    return true;
}


//--------------------------------------------------------------------------------------
// Create any D3D11 resources that aren't dependant on the back buffer
//--------------------------------------------------------------------------------------
HRESULT CALLBACK OnD3D11CreateDevice( ID3D11Device* pd3dDevice, const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc,
                                      void* pUserContext )
{
	HRESULT hr = S_OK;

	auto pContext = DXUTGetD3D11DeviceContext();

	V_RETURN(RenderStates::Initialize(pd3dDevice));

	g_DirectionalLight.Direction = D3DXVECTOR3(-10.0f, -3.0f, -1.0f);
	g_DirectionalLight.Ambient = D3DXVECTOR4(0.05f, 0.05f, 0.05f, 1.0f);
	g_DirectionalLight.Diffuse = D3DXVECTOR4(0.70f, 0.70f, 0.70f, 1.0f);
	g_DirectionalLight.Specular = D3DXVECTOR4(0.70f, 0.70f, 0.70f, 0.70f);

	D3DXMatrixIdentity(&g_World);
	m_Camera.SetViewParams(&g_Eye, &g_At);
	g_View = *(m_Camera.GetViewMatrix());

	m_pModel = new Model();
	m_pModel->Initialize();
	m_pModel->SetModelHeight(m_EyeHeight);

	m_pFrameBuffer = new FrameBuffer();
	m_pFrameBuffer->Initialize();

	m_pShadowMapFrameBuffer = new FrameBuffer();
	m_pShadowMapFrameBuffer->Initialize();

	m_pAtmosphere = new Atmosphere();
	m_pAtmosphere->Initialize();

	V_RETURN(m_pFrameBuffer->OnD3D11CreateDevice(pd3dDevice, pContext));
	V_RETURN(m_pShadowMapFrameBuffer->OnD3D11CreateDevice(pd3dDevice, pContext));
	V_RETURN(m_pModel->OnD3D11CreateDevice(pd3dDevice, pContext));
	V_RETURN(m_pAtmosphere->OnD3D11CreateDevice(pd3dDevice,pContext));

	m_pShadowMapFrameBuffer->Resize(m_shadowMapDim, m_shadowMapDim);

    return S_OK;
}


//--------------------------------------------------------------------------------------
// Create any D3D11 resources that depend on the back buffer
//--------------------------------------------------------------------------------------
HRESULT CALLBACK OnD3D11ResizedSwapChain( ID3D11Device* pd3dDevice, IDXGISwapChain* pSwapChain,
                                          const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext )
{
	HRESULT hr;
	screen_width = pBackBufferSurfaceDesc->Width;
	screen_height = pBackBufferSurfaceDesc->Height;
	g_viewport.Width = static_cast<float>(screen_width);
	g_viewport.Height = static_cast<float>(screen_height);
	float fAspect = g_viewport.Width / g_viewport.Height;
	float fNear = 0.1, fFar = 1000;

	g_pRenderTargetView = DXUTGetD3D11RenderTargetView();
	g_pDepthStencilView = DXUTGetD3D11DepthStencilView();

	m_Camera.SetProjParams(D3DX_PI * 0.25f, fAspect, fNear, fFar);
	g_Projection = *(m_Camera.GetProjMatrix());

	m_pAtmosphere->Resize(screen_width, screen_height, D3DX_PI * 25.0f, fAspect, fNear, fFar);
	m_pModel->Resize(pBackBufferSurfaceDesc);
	m_pFrameBuffer->Resize(screen_width, screen_height);
    return S_OK;
}


//--------------------------------------------------------------------------------------
// Handle updates to the scene.  This is called regardless of which D3D API is used
//--------------------------------------------------------------------------------------
void CALLBACK OnFrameMove( double fTime, float fElapsedTime, void* pUserContext )
{
	m_pAtmosphere->OnFrameMove(fTime,fElapsedTime);
	m_Camera.FrameMove(fElapsedTime);
}


//--------------------------------------------------------------------------------------
// Render the scene using the D3D11 device
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D11FrameRender( ID3D11Device* pd3dDevice, ID3D11DeviceContext* pd3dImmediateContext,
                                  double fTime, float fElapsedTime, void* pUserContext )
{
	float ClearColor[4] = { 0.176f, 0.196f, 0.667f, 0.0f };

	ID3D11RenderTargetView* pRTV = DXUTGetD3D11RenderTargetView();
	ID3D11DepthStencilView* pDSV = DXUTGetD3D11DepthStencilView();
	pd3dImmediateContext->ClearRenderTargetView(pRTV, ClearColor);
	pd3dImmediateContext->ClearDepthStencilView(pDSV, D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 1.0, 0);

	m_pAtmosphere->PreCompute(pd3dDevice, pd3dImmediateContext,pRTV);
	{
		m_pFrameBuffer->Activate();
		m_pFrameBuffer->ActivateDepth(true);

		pd3dImmediateContext->RSSetState(RenderStates::CullClockWiseRS);
		UpdateMatrix();

		D3DXMATRIX viewProj = g_View*g_Projection;
		m_pModel->SetViewProj(viewProj);
		m_pModel->SetViewPos(g_Eye);
		m_pModel->SetModelHeight(g_Eye.y);
		m_pModel->SetLight(&g_DirectionalLight);
		m_pModel->Render(pd3dDevice, pd3dImmediateContext);

		m_pFrameBuffer->DeactivateDepth();
	}
	{
		m_pShadowMapFrameBuffer->Activate();
		m_pShadowMapFrameBuffer->ActivateDepth(true);
		
		//pd3dImmediateContext->RSSetState(RenderStates::NoCullRS);
		m_pModel->RenderShadowMap(pd3dDevice, pd3dImmediateContext, m_pAtmosphere->GetSunDir(),m_shadowMapDim);
		pd3dImmediateContext->RSSetState(RenderStates::CullCounterClockWiseRS);

		m_pShadowMapFrameBuffer->DeactivateDepth();
	}
	m_pAtmosphere->SetCamParam(g_Eye, g_CamDir, g_View, g_Projection);
	m_pAtmosphere->Render(pd3dDevice, pd3dImmediateContext, pRTV, 
		m_pFrameBuffer->GetDepthSRV(), m_pShadowMapFrameBuffer->GetDepthSRV(),
		m_shadowMapDim);
}


//--------------------------------------------------------------------------------------
// Release D3D11 resources created in OnD3D11ResizedSwapChain 
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D11ReleasingSwapChain( void* pUserContext )
{
}


//--------------------------------------------------------------------------------------
// Release D3D11 resources created in OnD3D11CreateDevice 
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D11DestroyDevice(void* pUserContext)
{
	RenderStates::Release();
	if (m_pAtmosphere)
	{
		m_pAtmosphere->Release();
		delete m_pAtmosphere;
		m_pAtmosphere = nullptr;
	}
	if(m_pFrameBuffer)
	{
		m_pFrameBuffer->Release();
		delete m_pFrameBuffer;
		m_pFrameBuffer = nullptr;
	}
	if (m_pShadowMapFrameBuffer)
	{
		m_pShadowMapFrameBuffer->Release();
		delete m_pShadowMapFrameBuffer;
		m_pShadowMapFrameBuffer = nullptr;
	}
	if(m_pModel)
	{
		m_pModel->Release();
		delete m_pModel;
		m_pModel = nullptr;
	}
}


//--------------------------------------------------------------------------------------
// Handle messages to the application
//--------------------------------------------------------------------------------------
LRESULT CALLBACK MsgProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam,
                          bool* pbNoFurtherProcessing, void* pUserContext )
{
	if (m_pAtmosphere)
	{
		m_pAtmosphere->MsgProc(hWnd, uMsg, wParam, lParam);
		g_DirectionalLight.Direction = m_pAtmosphere->GetSunDir();
	}
	m_Camera.HandleMessages(hWnd, uMsg, wParam, lParam);

    return 0;
}


//--------------------------------------------------------------------------------------
// Call if device was removed.  Return true to find a new device, false to quit
//--------------------------------------------------------------------------------------
bool CALLBACK OnDeviceRemoved( void* pUserContext )
{
    return true;
}


//--------------------------------------------------------------------------------------
// Initialize everything and go into a render loop
//--------------------------------------------------------------------------------------
int WINAPI wWinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance, LPWSTR lpCmdLine, int nCmdShow )
{
    // Enable run-time memory check for debug builds.
#if defined(DEBUG) | defined(_DEBUG)
    _CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	//_CrtSetBreakAlloc(267);
#endif

    // DXUT will create and use the best device (either D3D9 or D3D11) 
    // that is available on the system depending on which D3D callbacks are set below

    // Set general DXUT callbacks
    DXUTSetCallbackFrameMove( OnFrameMove );
    DXUTSetCallbackMsgProc( MsgProc );
    DXUTSetCallbackDeviceChanging( ModifyDeviceSettings );
    DXUTSetCallbackDeviceRemoved( OnDeviceRemoved );

    // Set the D3D11 DXUT callbacks. Remove these sets if the app doesn't need to support D3D11
    DXUTSetCallbackD3D11DeviceAcceptable( IsD3D11DeviceAcceptable );
    DXUTSetCallbackD3D11DeviceCreated( OnD3D11CreateDevice );
    DXUTSetCallbackD3D11SwapChainResized( OnD3D11ResizedSwapChain );
    DXUTSetCallbackD3D11FrameRender( OnD3D11FrameRender );
    DXUTSetCallbackD3D11SwapChainReleasing( OnD3D11ReleasingSwapChain );
    DXUTSetCallbackD3D11DeviceDestroyed( OnD3D11DestroyDevice );

    // Perform any application-level initialization here

    DXUTInit( true, true, nullptr ); // Parse the command line, show msgboxes on error, no extra command line params
    DXUTSetCursorSettings( true, true ); // Show the cursor and clip it when in full screen
    DXUTCreateWindow( L"Atmosphere" );

    // Only require 10-level hardware
    DXUTCreateDevice( D3D_FEATURE_LEVEL_10_0, true, screen_width, screen_height );
    DXUTMainLoop(); // Enter into the DXUT ren  der loop

    // Perform any application-level cleanup here

    return DXUTGetExitCode();
}


