#include "DXUT.h"
#include "DXUTcamera.h"
#include "SDKmisc.h"
#include "Atmosphere.h"
#include "RenderStates.h"

#ifdef _DEBUG 
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <malloc.h>    // 解决 malloc.h 与 crtdbg.h 顺序导致的 Debug Assertion Failed, "Corrupted pointer passed to _freea" 。
#include <crtdbg.h>
#define new new( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#endif  // _DEBUG


CFirstPersonCamera								mCamera;
ID3D11Device*									gd3dDevice;
ID3D11DeviceContext*							gd3dImmediateContext;
Atmosphere*										mAtmosphere;
float	fCameraHeight = 6357.0f;
D3DXVECTOR3 v3CameraPos = D3DXVECTOR3(0,0, fCameraHeight);
D3DXVECTOR3 LookAt = D3DXVECTOR3(1, 0, fCameraHeight);
int	screen_width = 1024, screen_height = 768;
float fAspect;
float fNear = 0.1f, fFar = 100.0f;
float fFov = D3DX_PI * 0.25f;

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
	gd3dDevice = pd3dDevice;
	gd3dImmediateContext = DXUTGetD3D11DeviceContext();

	RenderStates::Initialize(gd3dDevice);

	mAtmosphere = new Atmosphere();
	mAtmosphere->Initialize(gd3dDevice,gd3dImmediateContext);
	mAtmosphere->OnD3D11CreateDevice();

	mCamera.SetViewParams(&v3CameraPos, &LookAt);

    return hr;
}


//--------------------------------------------------------------------------------------
// Create any D3D11 resources that depend on the back buffer
//--------------------------------------------------------------------------------------
HRESULT CALLBACK OnD3D11ResizedSwapChain( ID3D11Device* pd3dDevice, IDXGISwapChain* pSwapChain,
                                          const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext )
{
	screen_width = pBackBufferSurfaceDesc->Width;
	screen_height = pBackBufferSurfaceDesc->Height;
	fAspect = static_cast<float>(screen_width) / static_cast<float>(screen_height);
	mCamera.SetProjParams(fFov, fAspect, fNear, fFar);


    return S_OK;
}


//--------------------------------------------------------------------------------------
// Handle updates to the scene.  This is called regardless of which D3D API is used
//--------------------------------------------------------------------------------------
void CALLBACK OnFrameMove( double fTime, float fElapsedTime, void* pUserContext )
{
	mCamera.FrameMove(fElapsedTime);
}


//--------------------------------------------------------------------------------------
// Render the scene using the D3D11 device
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D11FrameRender( ID3D11Device* pd3dDevice, ID3D11DeviceContext* pd3dImmediateContext,
                                  double fTime, float fElapsedTime, void* pUserContext )
{
    // Clear render target and the depth stencil 
    float ClearColor[4] = { 0.176f, 0.196f, 0.667f, 0.0f };

    ID3D11RenderTargetView* pRTV = DXUTGetD3D11RenderTargetView();
    ID3D11DepthStencilView* pDSV = DXUTGetD3D11DepthStencilView();
    pd3dImmediateContext->ClearRenderTargetView( pRTV, ClearColor );
	pd3dImmediateContext->ClearDepthStencilView(pDSV, D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 1.0, 0);

	v3CameraPos = *mCamera.GetEyePt();
	D3DXMATRIX viewPro;
	D3DXMatrixMultiply(&viewPro, mCamera.GetViewMatrix(), mCamera.GetProjMatrix());
	
	D3DXVECTOR3 lookAhead = *mCamera.GetWorldAhead();
	D3DXVECTOR3 lookUp = *mCamera.GetWorldUp();
	D3DXVECTOR3 lookRight = *mCamera.GetWorldRight();
	
	lookAhead *= (fNear + fFar) / 2;
	float height_half = tanf(fFov / 2) * (fNear + fFar) / 2;
	float width_half = height_half * fAspect;
	
	std::vector<Vertex> v=
	{
		{ v3CameraPos + lookAhead - height_half * lookUp - width_half * lookRight,{0.f,0.f,0.f},{ 0.f,0.f,0.f } ,{ 0.f,0.f }},
		{ v3CameraPos + lookAhead + height_half * lookUp - width_half * lookRight,{ 0.f,0.f,0.f },{ 0.f,0.f,0.f } ,{ 0.f,0.f }},
		{ v3CameraPos + lookAhead - height_half * lookUp + width_half * lookRight,{ 0.f,0.f,0.f },{ 0.f,0.f,0.f } ,{ 0.f,0.f }},
		{ v3CameraPos + lookAhead + height_half * lookUp + width_half * lookRight,{ 0.f,0.f,0.f },{ 0.f,0.f,0.f } ,{ 0.f,0.f}}
	};

	std::vector<GeometryGenerator::uint32> index =
	{
		0,1,2,
		2,1,3
	};

	GeometryGenerator::MeshData mesh;
	mesh.Vertices = v;
	mesh.Indices32 = index;

	mAtmosphere->Render(viewPro, v3CameraPos,mesh);
	
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
void CALLBACK OnD3D11DestroyDevice( void* pUserContext )
{
	gd3dDevice = nullptr;
	gd3dImmediateContext = nullptr;	
	if(mAtmosphere)
	{
		mAtmosphere->Release();
		delete mAtmosphere;
	}
	
	RenderStates::Release();
}


//--------------------------------------------------------------------------------------
// Handle messages to the application
//--------------------------------------------------------------------------------------
LRESULT CALLBACK MsgProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam,
                          bool* pbNoFurtherProcessing, void* pUserContext )
{
	mCamera.HandleMessages(hWnd,uMsg,wParam,lParam);
    return 0;
}


//--------------------------------------------------------------------------------------
// Handle key presses
//--------------------------------------------------------------------------------------
void CALLBACK OnKeyboard( UINT nChar, bool bKeyDown, bool bAltDown, void* pUserContext )
{
}


//--------------------------------------------------------------------------------------
// Handle mouse button presses
//--------------------------------------------------------------------------------------
void CALLBACK OnMouse( bool bLeftButtonDown, bool bRightButtonDown, bool bMiddleButtonDown,
                       bool bSideButton1Down, bool bSideButton2Down, int nMouseWheelDelta,
                       int xPos, int yPos, void* pUserContext )
{
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
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);
    _CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

    // DXUT will create and use the best device (either D3D9 or D3D11) 
    // that is available on the system depending on which D3D callbacks are set below

    // Set general DXUT callbacks
    DXUTSetCallbackFrameMove( OnFrameMove );
    DXUTSetCallbackKeyboard( OnKeyboard );
    DXUTSetCallbackMouse( OnMouse );
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

    DXUTInit( true, true, NULL ); // Parse the command line, show msgboxes on error, no extra command line params
    DXUTSetCursorSettings( true, true ); // Show the cursor and clip it when in full screen
    DXUTCreateWindow( L"Atmosphere" );

    // Only require 10-level hardware
    DXUTCreateDevice( D3D_FEATURE_LEVEL_10_0, true, screen_width, screen_height );
    DXUTMainLoop(); // Enter into the DXUT ren  der loop

#if defined(DEBUG) | defined(_DEBUG)
	_CrtDumpMemoryLeaks();
#endif
    // Perform any application-level cleanup here
    return DXUTGetExitCode();
}
