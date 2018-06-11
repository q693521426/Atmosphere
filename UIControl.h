#pragma once
#ifndef UICONTROL_H
#define UICONTROL_H

#include "DXUT.h"
#include "DXUTgui.h"
#include "DXUTsettingsDlg.h"
#include "SDKmisc.h"

class UIControl
{
public:
	UIControl();
	~UIControl();

	bool g_IsRenderModel = true;
	bool g_UseEpipolarLine = true;

	HRESULT OnD3D11CreateDevice(ID3D11Device*, ID3D11DeviceContext*);
	void Initialize(PCALLBACKDXUTGUIEVENT);
	void Release();

	void Render(float fElapsedTime);

	HRESULT OnD3D11ResizedSwapChain(ID3D11Device* pd3dDevice, IDXGISwapChain* pSwapChain,
		const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext);

	LRESULT MsgProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam,
			bool* pbNoFurtherProcessing, void* pUserContext);

private:
	CDXUTDialogResourceManager	m_DialogResourceManager;
	CD3DSettingsDlg				m_D3DSettingsDlg;
	CDXUTDialog					m_AtmosphereParaUI;
	CDXUTDialog					m_SwitchUI;

	void CALLBACK OnGUIEvent(UINT nEvent, int nControlID, CDXUTControl* pControl, void* pUserContext);
};

#endif