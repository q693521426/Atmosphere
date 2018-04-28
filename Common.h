#pragma once

#ifndef COMMON_H
#define COMMON_H

#include "DXUT.h"
#include "SDKmisc.h"
#include <algorithm>
#include <string>
#include <vector>
#include <atlcomcli.h>
#include <unordered_map>
#include <sstream>
#include <set>

#define D3D_COMPILE_STANDARD_FILE_INCLUDE ((ID3DInclude*)(UINT_PTR)1)

struct MiscDynamicParams
{
	D3DXVECTOR2 f2WQ;
	float scatter_order;
	UINT uiMinMaxLevelMax;

	UINT uiSrcMinMaxOffsetX;
	UINT uiSrcMinMaxOffsetY;
	UINT uiDstMinMaxOffsetX;
	UINT uiDstMinMaxOffsetY;

	float fEnableLightShaft;
	float padding[3];
};

struct CameraParams
{
	D3DXVECTOR3 f3CameraPos;
	float fNearZ;

	D3DXVECTOR3 f3CameraDir;
	float fFarZ;

	D3DXMATRIX View;
	D3DXMATRIX Proj;
	D3DXMATRIX ViewProj;
	D3DXMATRIX InvViewProj;
};

struct LightParams
{
	D3DXVECTOR3 f3LightDir;
	float padding;

	D3DXVECTOR4 f4LightScreenPos;

	D3DXMATRIX View;
	D3DXMATRIX Proj;
	D3DXMATRIX ViewProj;
	D3DXMATRIX InvViewProj;
};

class GameObject
{
public:
	GameObject(){}
	~GameObject(){}
	HRESULT OnD3D11CreateDevice(ID3D11Device*, ID3D11DeviceContext*, 
				std::wstring EffectStr, const std::vector<std::string>& TechStr,
				const std::vector<std::string>& VarStr,
				const std::vector<std::string>& ShaderResourceVarStr);
	void Release();
protected:
	CComPtr<ID3DX11Effect>		pEffect;

	std::unordered_map<std::string, CComPtr<ID3DX11EffectTechnique>>				TechMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectMatrixVariable>>			MatrixVarMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectVectorVariable>>			VectorVarMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectScalarVariable>>			ScalarVarMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectVariable>>					VarMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectShaderResourceVariable>>	ShaderResourceVarMap;

	MiscDynamicParams misc;
private:

};

struct Vertex
{
	Vertex() {}
	Vertex(
		const D3DXVECTOR3& p,
		const D3DXVECTOR3& n,
		const D3DXVECTOR3& t,
		const D3DXVECTOR2& uv) :
		Position(p),
		Normal(n),
		TangentU(t),
		TexC(uv) {}
	Vertex(
		float px, float py, float pz,
		float nx, float ny, float nz,
		float tx, float ty, float tz,
		float u, float v) :
		Position(px, py, pz),
		Normal(nx, ny, nz),
		TangentU(tx, ty, tz),
		TexC(u, v) {}

	D3DXVECTOR3 Position;
	D3DXVECTOR3 Normal;
	D3DXVECTOR3 TangentU;
	D3DXVECTOR2 TexC;
};

inline HRESULT CompileShaderFromFile(LPCTSTR strFilePath,
	LPCSTR strFunctionName,
	const D3D_SHADER_MACRO* pDefines,
	LPCSTR profile,
	ID3DBlob **ppBlobOut)
{
	DWORD dwShaderFlags = D3D10_SHADER_ENABLE_STRICTNESS;
#if defined( DEBUG ) || defined( _DEBUG )
	dwShaderFlags |= D3D10_SHADER_DEBUG;
#else
#endif
	HRESULT hr;
	do
	{
		CComPtr<ID3DBlob> errors;
		hr = D3DX11CompileFromFile(strFilePath, pDefines, NULL, strFunctionName, profile, dwShaderFlags, 0, NULL, ppBlobOut, &errors, NULL);
		if (errors)
		{
			OutputDebugStringA((char*)errors->GetBufferPointer());
			if (FAILED(hr) &&
				IDRETRY != MessageBoxA(NULL, (char*)errors->GetBufferPointer(), "FX Error", MB_ICONERROR | MB_ABORTRETRYIGNORE))
			{
				break;
			}
		}
	} while (FAILED(hr));
	return hr;
}

inline HRESULT CompileEffectFromFile(ID3D11Device* pd3dDevice, ID3DX11Effect** pEffect, wchar_t* FileName)
{
	HRESULT hr = S_OK;
	DWORD dwShaderFlags = D3DCOMPILE_ENABLE_STRICTNESS;
#ifdef _DEBUG
	dwShaderFlags |= D3DCOMPILE_DEBUG;
	dwShaderFlags |= D3DCOMPILE_SKIP_OPTIMIZATION;
#endif

	ID3DBlob* pEffectBuffer = nullptr;

	WCHAR szShaderPath[MAX_PATH];
	V_RETURN(DXUTFindDXSDKMediaFileCch(szShaderPath, MAX_PATH, FileName));

	ID3DBlob* pErrorBlob = nullptr;
	hr = D3DX11CompileFromFile(szShaderPath, nullptr, nullptr, nullptr, "fx_5_0", dwShaderFlags, 0, nullptr, &pEffectBuffer, &pErrorBlob, nullptr);
	if (pErrorBlob)
	{
		OutputDebugStringA(reinterpret_cast<const char*>(pErrorBlob->GetBufferPointer()));
		pErrorBlob->Release();
		pErrorBlob = nullptr;
	}
	if (FAILED(hr))
		return hr;

	hr = D3DX11CreateEffectFromMemory(pEffectBuffer->GetBufferPointer(), pEffectBuffer->GetBufferSize(), 0, pd3dDevice, pEffect);
	//hr = D3DX11CompileEffectFromFile( szShaderPath, nullptr, D3D_COMPILE_STANDARD_FILE_INCLUDE, dwShaderFlags, 0, pd3dDevice, pEffect, &pErrorBlob );

	if (pErrorBlob)
	{
		OutputDebugStringA(reinterpret_cast<const char*>(pErrorBlob->GetBufferPointer()));
		pErrorBlob->Release();
	}

	if (FAILED(hr))
		return hr;

	SAFE_RELEASE(pEffectBuffer);

	return hr;
}

//std::set<std::string> m_DefinitionsSet;
//
//template<class T>
//void AddMarco(std::vector<D3D_SHADER_MACRO>& Marcos, LPCSTR Name,T Val)
//{
//	std::string Definition = std::to_string(Val);
//	D3D_SHADER_MACRO Marco =
//	{
//		Name,
//		m_DefinitionsSet.insert(Definition).first->c_str()
//	};
//	Marcos.push_back(Marco);
//}
//
//inline
//void FinishMarco(std::vector<D3D_SHADER_MACRO>& Marcos)
//{
//	D3D_SHADER_MACRO Marco =
//	{
//		nullptr,
//		nullptr
//	};
//	Marcos.push_back(Marco);
//	m_DefinitionsSet.clear();
//}

inline void UnbindResources(ID3D11DeviceContext *pDeviceCtx)
{
	ID3D11ShaderResourceView *pDummySRVs[12] = { nullptr };
	ID3D11UnorderedAccessView *pDummyUAVs[8] = { nullptr };
	pDeviceCtx->PSSetShaderResources(0, _countof(pDummySRVs), pDummySRVs);
	pDeviceCtx->VSSetShaderResources(0, _countof(pDummySRVs), pDummySRVs);
	pDeviceCtx->GSSetShaderResources(0, _countof(pDummySRVs), pDummySRVs);
	pDeviceCtx->CSSetShaderResources(0, _countof(pDummySRVs), pDummySRVs);
	pDeviceCtx->CSSetUnorderedAccessViews(0, _countof(pDummyUAVs), pDummyUAVs, nullptr);
}

inline D3DXMATRIX InverseTranspose(D3DXMATRIX* m)
{
	D3DXMATRIX out;
	float det = D3DXMatrixDeterminant(m);
	D3DXMatrixInverse(&out, &det, m);
	D3DXMatrixTranspose(&out, &out);
	return out;
}

inline HRESULT CreateTexture2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext,
								UINT width,UINT height, DXGI_FORMAT format, 
								std::vector<ID3D11Texture2D**> ppTex2Ds,
								std::vector<ID3D11ShaderResourceView**> ppSRVs,
								std::vector<ID3D11RenderTargetView**> ppRTVs)
{
	HRESULT hr = S_OK;

	D3D11_TEXTURE2D_DESC PreCompute2DTexDesc;
	ZeroMemory(&PreCompute2DTexDesc, sizeof(PreCompute2DTexDesc));
	PreCompute2DTexDesc.Width = width;
	PreCompute2DTexDesc.Height = height;
	PreCompute2DTexDesc.MipLevels = 1;
	PreCompute2DTexDesc.ArraySize = 1;
	PreCompute2DTexDesc.Format = format;
	PreCompute2DTexDesc.SampleDesc.Count = 1;
	PreCompute2DTexDesc.SampleDesc.Quality = 0;
	PreCompute2DTexDesc.Usage = D3D11_USAGE_DEFAULT;
	PreCompute2DTexDesc.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
	PreCompute2DTexDesc.CPUAccessFlags = 0;
	PreCompute2DTexDesc.MiscFlags = 0;

	D3D11_RENDER_TARGET_VIEW_DESC PreComputeRTVDesc;
	PreComputeRTVDesc.Format = PreCompute2DTexDesc.Format;
	PreComputeRTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
	PreComputeRTVDesc.Texture2D.MipSlice = 0;

	D3D11_SHADER_RESOURCE_VIEW_DESC PreComputeSRVDesc;
	PreComputeSRVDesc.Format = PreCompute2DTexDesc.Format;
	PreComputeSRVDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
	PreComputeSRVDesc.Texture2D.MostDetailedMip = 0;
	PreComputeSRVDesc.Texture2D.MipLevels = 1;

	for(int i = 0 ; i<ppTex2Ds.size(); ++i)
	{
		SAFE_RELEASE((*ppTex2Ds[i]));
		SAFE_RELEASE((*ppRTVs[i]));
		SAFE_RELEASE((*ppSRVs[i]));
		V_RETURN(pDevice->CreateTexture2D(&PreCompute2DTexDesc, nullptr, ppTex2Ds[i]));
		V_RETURN(pDevice->CreateRenderTargetView(*ppTex2Ds[i], &PreComputeRTVDesc, ppRTVs[i]));
		V_RETURN(pDevice->CreateShaderResourceView(*ppTex2Ds[i], &PreComputeSRVDesc, ppSRVs[i]));
	}

	return hr;
}

inline HRESULT CreateTexture3D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext,
	UINT width, UINT height, UINT depth, DXGI_FORMAT format, 
	std::vector<ID3D11Texture3D**> ppTex3Ds,std::vector<ID3D11ShaderResourceView**> ppSRVs)
{
	HRESULT hr = S_OK;

	D3D11_TEXTURE3D_DESC PreCompute3DTexDesc;
	ZeroMemory(&PreCompute3DTexDesc, sizeof(PreCompute3DTexDesc));
	PreCompute3DTexDesc.Width = width;
	PreCompute3DTexDesc.Height = height;
	PreCompute3DTexDesc.Depth = depth;
	PreCompute3DTexDesc.MipLevels = 1;
	PreCompute3DTexDesc.Format = format;
	PreCompute3DTexDesc.Usage = D3D11_USAGE_DEFAULT;
	PreCompute3DTexDesc.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
	PreCompute3DTexDesc.CPUAccessFlags = 0;
	PreCompute3DTexDesc.MiscFlags = 0;

	for(int i=0;i<ppTex3Ds.size();++i)
	{
		SAFE_RELEASE((*ppTex3Ds[i]));
		SAFE_RELEASE((*ppSRVs[i]));
		V_RETURN(pDevice->CreateTexture3D(&PreCompute3DTexDesc, nullptr, ppTex3Ds[i]));
		V_RETURN(pDevice->CreateShaderResourceView(*ppTex3Ds[i], nullptr, ppSRVs[i]));
	}

	return hr;
}

inline void RenderQuad(ID3D11DeviceContext* pContext,
	ID3DX11EffectTechnique* activeTech,
	int iWidth = 0, int iHeight = 0,
	int iTopLeftX = 0, int iTopLeftY = 0,
	int iNumInstances = 1)
{
	D3D11_VIEWPORT ViewPort;
	ViewPort.TopLeftX = static_cast<float>(iTopLeftX);
	ViewPort.TopLeftY = static_cast<float>(iTopLeftY);
	ViewPort.Width = static_cast<float>(iWidth);
	ViewPort.Height = static_cast<float>(iHeight);
	ViewPort.MinDepth = 0;
	ViewPort.MaxDepth = 1;
	// Set the viewport
	pContext->RSSetViewports(1, &ViewPort);

	D3DX11_TECHNIQUE_DESC techDesc;
	activeTech->GetDesc(&techDesc);
	for (UINT p = 0; p<techDesc.Passes; ++p)
	{
		CComPtr<ID3DX11EffectPass> pass = activeTech->GetPassByIndex(p);

		UINT offset[1] = { 0 };
		UINT stride[1] = { 0 };
		ID3D11Buffer *ppBuffers[1] = { 0 };
		pContext->IASetVertexBuffers(0, 1, ppBuffers, stride, offset);
		pContext->IASetInputLayout(nullptr);
		pContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);
		pass->Apply(0, pContext);
		if (iNumInstances == 1)
		{
			// Draw 4 vertices (two triangles )
			pContext->Draw(4, 0);
		}
		else
		{
			// Draw 4 vertices (two triangles ) x number of instances
			pContext->DrawInstanced(4, iNumInstances, 0, 0);
		}
	}
}

template<typename T1,typename T2>
void MapRelease(std::unordered_map<T1,CComPtr<T2>>& m)
{
	std::for_each(m.begin(), m.end(), [](std::pair<const T1, CComPtr<T2>>& v) {v.second.Release(); });
	m.clear();
}

template<typename T1>
void VectorRelease(std::vector<T1>& v)
{
	v.clear();
	v.shrink_to_fit();
	std::vector<T1>().swap(v);
}

#endif