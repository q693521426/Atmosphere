#include "DXUT.h"
#include "Cloud.h"
#include <sstream>

Cloud::Cloud()
{
}


Cloud::~Cloud()
{
}


void Cloud::Initialize()
{
	
}


void Cloud::Release()
{
	pPerlinWorleyTex3D.Release();
	pPerlinWorleySRV.Release();
	pWorleyTex3D.Release();
	pWorleySRV.Release();
	pCurlTex3D.Release();
	pCurlSRV.Release();

	GameObject::Release();
}


HRESULT Cloud::OnD3D11CreateDevice(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	V_RETURN(GameObject::OnD3D11CreateDevice(pDevice, pContext, L"Cloud.fx", TechStr, VarStr, ShaderResourceVarStr));

	return hr;
}


HRESULT Cloud::PreCompute(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11RenderTargetView* pRTV)
{
	HRESULT hr = S_OK;

	V_RETURN(PreComputePerlinWorleyTex3D(pDevice, pContext));
	V_RETURN(PreComputeWorleyTex3D(pDevice, pContext));
	
	return hr;
}


HRESULT Cloud::PreComputePerlinWorleyTex3D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	pPerlinWorleyTex3D.Release();
	pPerlinWorleySRV.Release();
	V_RETURN(CreateTexture3D(pDevice, pContext, PERLIN_WORLEY_TEXTURE_DIM, PERLIN_WORLEY_TEXTURE_DIM, PERLIN_WORLEY_TEXTURE_DIM, format,
	{ &pPerlinWorleyTex3D.p }, { &pPerlinWorleySRV.p }));

	std::vector<CComPtr<ID3D11RenderTargetView>> pNoiseRTVs(PERLIN_WORLEY_TEXTURE_DIM);
	MiscDynamicParams misc;

	for (UINT depthSlice = 0; depthSlice < PERLIN_WORLEY_TEXTURE_DIM; ++depthSlice)
	{
		D3D11_RENDER_TARGET_VIEW_DESC CurrSliceRTVDesc;
		CurrSliceRTVDesc.Format = format;
		CurrSliceRTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
		CurrSliceRTVDesc.Texture3D.MipSlice = 0;
		CurrSliceRTVDesc.Texture3D.FirstWSlice = depthSlice;
		CurrSliceRTVDesc.Texture3D.WSize = 1;

		V_RETURN(pDevice->CreateRenderTargetView(pPerlinWorleyTex3D, &CurrSliceRTVDesc, &pNoiseRTVs[depthSlice]));

		ID3DX11EffectTechnique* activeTech = TechMap["ComputePerlinWorleyNoiseTex3DTech"];

		misc.f2WQ.x = (depthSlice + 0.5) / PERLIN_WORLEY_TEXTURE_DIM;
		assert(0 < misc.f2WQ.x && misc.f2WQ.x < 1);

		VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));

		ID3D11RenderTargetView* pRTVs[] =
		{
			pNoiseRTVs[depthSlice].p
		};

		UINT size = ARRAYSIZE(pRTVs);
		for (UINT i = 0; i<size; ++i)
		{
			float zero[] = { 0.f, 0.f,0.f,0.f };
			pContext->ClearRenderTargetView(pRTVs[i], zero);
		}
		pContext->OMSetRenderTargets(size, pRTVs, nullptr);

		RenderQuad(pContext, activeTech, PERLIN_WORLEY_TEXTURE_DIM, PERLIN_WORLEY_TEXTURE_DIM);
	}

#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pPerlinWorleyTex3D, D3DX11_IFF_DDS, L"Texture/PerlinWorley.dds"));
#endif
}


HRESULT Cloud::PreComputeWorleyTex3D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	pWorleyTex3D.Release();
	pWorleySRV.Release();
	V_RETURN(CreateTexture3D(pDevice, pContext, WORLEY_TEXTURE_DIM, WORLEY_TEXTURE_DIM, WORLEY_TEXTURE_DIM, format,
	{ &pWorleyTex3D.p }, { &pWorleySRV.p }));

	std::vector<CComPtr<ID3D11RenderTargetView>> pNoiseRTVs(WORLEY_TEXTURE_DIM);
	MiscDynamicParams misc;

	for (UINT depthSlice = 0; depthSlice < WORLEY_TEXTURE_DIM; ++depthSlice)
	{
		D3D11_RENDER_TARGET_VIEW_DESC CurrSliceRTVDesc;
		CurrSliceRTVDesc.Format = format;
		CurrSliceRTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
		CurrSliceRTVDesc.Texture3D.MipSlice = 0;
		CurrSliceRTVDesc.Texture3D.FirstWSlice = depthSlice;
		CurrSliceRTVDesc.Texture3D.WSize = 1;

		V_RETURN(pDevice->CreateRenderTargetView(pWorleyTex3D, &CurrSliceRTVDesc, &pNoiseRTVs[depthSlice]));

		ID3DX11EffectTechnique* activeTech = TechMap["ComputeWorleyNoiseTex3DTech"];

		misc.f2WQ.x = (depthSlice + 0.5) / WORLEY_TEXTURE_DIM;
		assert(0 < misc.f2WQ.x && misc.f2WQ.x < 1);

		VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));

		ID3D11RenderTargetView* pRTVs[] =
		{
			pNoiseRTVs[depthSlice].p
		};

		UINT size = ARRAYSIZE(pRTVs);
		for (UINT i = 0; i<size; ++i)
		{
			float zero[] = { 0.f, 0.f,0.f,0.f };
			pContext->ClearRenderTargetView(pRTVs[i], zero);
		}
		pContext->OMSetRenderTargets(size, pRTVs, nullptr);

		RenderQuad(pContext, activeTech, WORLEY_TEXTURE_DIM, WORLEY_TEXTURE_DIM);
	}

#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pWorleyTex3D, D3DX11_IFF_DDS, L"Texture/Worley.dds"));
#endif
}


void Cloud::Render(ID3D11Device*, ID3D11DeviceContext*, ID3D11RenderTargetView*)
{
	
}

