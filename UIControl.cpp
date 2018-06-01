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


void UIControl::Initialize()
{

}


void UIControl::Release()
{
	m_DialogResourceManager.OnD3D11DestroyDevice();
	m_D3DSettingsDlg.OnD3D11DestroyDevice();

}

LRESULT UIControl::MsgProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam,
					bool* pbNoFurtherProcessing, void* pUserContext)
{

}