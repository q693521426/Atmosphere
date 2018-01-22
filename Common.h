#pragma once

#ifndef COMMON_H
#define COMMON_H

#include "DXUT.h"
#include "SDKmisc.h"
#include <algorithm>

#define D3D_COMPILE_STANDARD_FILE_INCLUDE ((ID3DInclude*)(UINT_PTR)1)

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

inline D3DXMATRIX InverseTranspose(D3DXMATRIX* m)
{
	D3DXMATRIX out;
	float det = D3DXMatrixDeterminant(m);
	D3DXMatrixInverse(&out, &det, m);
	D3DXMatrixTranspose(&out, &out);
	return out;
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