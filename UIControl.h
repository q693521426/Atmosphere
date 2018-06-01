#pragma once
#ifndef UICONTROL_H
#define UICONTROL_H

#include "DXUT.h"
#include "DXUTgui.h"
#include "DXUTsettingsDlg.h"

class UIControl
{
public:
	UIControl();
	~UIControl();

	HRESULT OnD3D11CreateDevice(ID3D11Device*, ID3D11DeviceContext*);
	void Initialize();
	void Release();

	LRESULT MsgProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam,
			bool* pbNoFurtherProcessing, void* pUserContext);
private:
	CDXUTDialogResourceManager	m_DialogResourceManager;
	CD3DSettingsDlg				m_D3DSettingsDlg;
	CDXUTDialog					m_AtmosphereParaUI;
};

#endif