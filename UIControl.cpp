#include "DXUT.h"
#include "UIControl.h"


UIControl::UIControl()
{
}


UIControl::~UIControl()
{
}


HRESULT UIControl::OnD3D11CreateDevice(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	V_RETURN(m_DialogResourceManager.OnD3D11CreateDevice(pDevice, pContext));
	V_RETURN(m_D3DSettingsDlg.OnD3D11CreateDevice(pDevice));

}


void UIControl::Initialize(PCALLBACKDXUTGUIEVENT OnGUIEvent)
{
	m_D3DSettingsDlg.Init(&m_DialogResourceManager);
	m_SwitchUI.Init(&m_DialogResourceManager);
	m_AtmosphereParaUI.Init(&m_DialogResourceManager);

	m_SwitchUI.SetCallback(OnGUIEvent);
}


void UIControl::Release()
{
	m_DialogResourceManager.OnD3D11DestroyDevice();
	m_D3DSettingsDlg.OnD3D11DestroyDevice();

}


void UIControl::Render(float fElapsedTime)
{
	DXUT_BeginPerfEvent(DXUT_PERFEVENTCOLOR, L"AtmosphereParams");
	m_AtmosphereParaUI.OnRender(fElapsedTime);
	m_SwitchUI.OnRender(fElapsedTime);
	DXUT_EndPerfEvent();
}


LRESULT UIControl::MsgProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam,
					bool* pbNoFurtherProcessing, void* pUserContext)
{
	*pbNoFurtherProcessing = m_DialogResourceManager.MsgProc(hWnd, uMsg, wParam, lParam);
	if (*pbNoFurtherProcessing)
		return 0;
	if (m_D3DSettingsDlg.IsActive())
	{
		m_DialogResourceManager.MsgProc(hWnd, uMsg, wParam, lParam);
		return 0;
	}
	*pbNoFurtherProcessing = m_SwitchUI.MsgProc(hWnd, uMsg, wParam, lParam);
	if (*pbNoFurtherProcessing)
		return 0;
	*pbNoFurtherProcessing = m_AtmosphereParaUI.MsgProc(hWnd, uMsg, wParam, lParam);
	if (*pbNoFurtherProcessing)
		return 0;
	return 0;
}


HRESULT UIControl::OnD3D11ResizedSwapChain(ID3D11Device* pd3dDevice, IDXGISwapChain* pSwapChain,
	const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext)
{
	HRESULT hr;
	V_RETURN(m_DialogResourceManager.OnD3D11ResizedSwapChain(pd3dDevice, pBackBufferSurfaceDesc));
	V_RETURN(m_D3DSettingsDlg.OnD3D11ResizedSwapChain(pd3dDevice, pBackBufferSurfaceDesc));

	m_SwitchUI.SetLocation(0, 0);
	m_SwitchUI.SetSize(200, 300);
	
	m_AtmosphereParaUI.SetLocation(pBackBufferSurfaceDesc->Width - 200, 0);
	m_AtmosphereParaUI.SetSize(200, 300);

	return S_OK;
}

void CALLBACK UIControl::OnGUIEvent(UINT nEvent, int nControlID, CDXUTControl* pControl, void* pUserContext)
{
	WCHAR sz[100];
	WCHAR m_wstr[100];

	swprintf_s(m_wstr, 100, L"");
	switch (nControlID)
	{

	}
}