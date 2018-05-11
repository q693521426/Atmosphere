#pragma once
#ifndef CLOUD_H
#define CLOUD_H

#include "DXUT.h"
#include "Common.h"

struct CloudTypeLayer
{
	D3DXVECTOR2 f2LayerHeightScale;
	D3DXVECTOR2 f2LayerDensityPoint;
};

struct CloudParams
{
	CloudTypeLayer mCloudTypeLayer[3];

	D3DXVECTOR2 f2CloudLayerHeightScale;
	float fTransition;
	float fUpperDensity;
};

class Cloud : public GameObject
{
public:
	Cloud();
	~Cloud();

	void Initialize();
	void Release();

	HRESULT OnD3D11CreateDevice(ID3D11Device*, ID3D11DeviceContext*);
	void Resize(int, int, float, float, float, float);

	HRESULT PreCompute(ID3D11Device*, ID3D11DeviceContext*, ID3D11RenderTargetView*);
	HRESULT PreComputePerlinWorleyTex3D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT PreComputeWorleyTex3D(ID3D11Device*, ID3D11DeviceContext*);
	//HRESULT PreComputePerlinWorleyTex_Test(ID3D11Device*, ID3D11DeviceContext*);

	void Render(ID3D11Device*, ID3D11DeviceContext*, ID3D11RenderTargetView*, ID3D11ShaderResourceView*, ID3D11ShaderResourceView*);

	void SetLightParams(LightParams* light);
	void SetCamParams(CameraParams* cam);
	void SetAtmosphereParams(AtmosphereParams* atmosphere);
private:
	const int PERLIN_WORLEY_TEXTURE_DIM = 128;
	const int WORLEY_TEXTURE_DIM = 32;
	const int CURL_TEXTURE_DIM = 128;

	float fLowerLayer = 1500;
	float fUpperLayer = 8000;

	bool IsPreComputed = true;
	CComPtr<ID3D11Texture3D>				pPerlinWorleyTex3D;
	CComPtr<ID3D11ShaderResourceView>		pPerlinWorleySRV;

	CComPtr<ID3D11Texture3D>				pWorleyTex3D;
	CComPtr<ID3D11ShaderResourceView>		pWorleySRV;

	CComPtr<ID3D11Texture3D>				pCurlTex3D;
	CComPtr<ID3D11ShaderResourceView>		pCurlSRV;

	CloudParams cloudParams;
	AtmosphereParams* pAtmosphereParams;
	CameraParams* pCameraParams;
	LightParams* pLightParams;

	std::vector<std::string> TechStr
	{
		"ComputePerlinWorleyNoiseTex3DTech",
		"ComputeWorleyNoiseTex3DTech",
		"DrawCloudTech"
	};

	std::vector<std::string> VarStr
	{
		"misc",
		"camera",
		"light",
		"atmosphere",
		"misc",
		"cloud"
	};

	std::vector<std::string> ShaderResourceVarStr
	{
		"g_tex3DPerlinWorleyNoise",
		"g_tex3DWorleyNoise",
		"g_tex2DSpaceDepth",
		"g_tex2DColorBuffer"
	};
};

#endif