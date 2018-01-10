#include "DXUT.h"
#include "Atmosphere.h"
#include "Common.h"

Atmosphere::Atmosphere()
{
}


Atmosphere::~Atmosphere()
{
}


void Atmosphere::Initialize()
{
	AtmosphereParams = AtmosphereParameters
	{
		
	};
}


void Atmosphere::Release()
{
	MapRelease(AtmosphereTechMap);
	MapRelease(MatrixVarMap);
	MapRelease(VectorVarMap);
	MapRelease(ScalarVarMap);
	MapRelease(VarMap);
	MapRelease(ShaderResourceVarMap);
	pAtmosphereEffect.Release();
	pTransmittanceTex2D.Release();
	pTransmittanceSRV.Release();
	pTransmittanceRTV.Release();
}


HRESULT Atmosphere::OnD3D11CreateDevice(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	V_RETURN(CompileEffectFromFile(pDevice, &pAtmosphereEffect, L"Atmosphere.fx"));

	for (int i = 0; i < TechStr.size(); i++)
	{
		CComPtr<ID3DX11EffectTechnique> pTech = pAtmosphereEffect->GetTechniqueByName(TechStr[i].c_str());
		AtmosphereTechMap.emplace(TechStr[i], pTech);
	}

	for (int i = 0; i < MatrixVarStr.size(); i++)
	{
		CComPtr<ID3DX11EffectMatrixVariable> pMatrixVar = pAtmosphereEffect->GetVariableByName(MatrixVarStr[i].c_str())->AsMatrix();
		MatrixVarMap.emplace(MatrixVarStr[i], pMatrixVar);
	}

	for (int i = 0; i < VectorVarStr.size(); i++)
	{
		CComPtr<ID3DX11EffectVectorVariable> pVectorVar = pAtmosphereEffect->GetVariableByName(VectorVarStr[i].c_str())->AsVector();
		VectorVarMap.emplace(VectorVarStr[i], pVectorVar);
	}

	for (int i = 0; i < ScalarVarStr.size(); i++)
	{
		CComPtr<ID3DX11EffectScalarVariable> pScalarVar = pAtmosphereEffect->GetVariableByName(ScalarVarStr[i].c_str())->AsScalar();
		ScalarVarMap.emplace(ScalarVarStr[i], pScalarVar);
	}

	for (int i = 0; i < VarStr.size(); i++)
	{
		CComPtr<ID3DX11EffectVariable> pVar = pAtmosphereEffect->GetVariableByName(VarStr[i].c_str());
		VarMap.emplace(VarStr[i], pVar);
	}

	for (int i = 0; i < ShaderResourceVarStr.size(); i++)
	{
		CComPtr<ID3DX11EffectShaderResourceVariable> pShaderResourceVar = pAtmosphereEffect->GetVariableByName(ShaderResourceVarStr[i].c_str())->AsShaderResource();
		ShaderResourceVarMap.emplace(ShaderResourceVarStr[i], pShaderResourceVar);
	}
	return hr;
}


void Atmosphere::Render(ID3D11Device* pDevice, ID3D11DeviceContext* pContext,ID3D11RenderTargetView* pRTV)
{
	//if(!IsPreComputed)
	//{
	//	PreComputeTransmittanceTex2D(pDevice, pContext,pRTV);
	//	PreComputeSingleSctrTex3D(pDevice, pContext);
	//	IsPreComputed = true;
	//}
	PreComputeTransmittanceTex2D(pDevice, pContext, pRTV);
}

HRESULT Atmosphere::PreComputeTransmittanceTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11RenderTargetView* pRTV = nullptr)
{
	HRESULT hr = S_OK;

	D3D11_TEXTURE2D_DESC PreCompute2DTexDesc;
	ZeroMemory(&PreCompute2DTexDesc, sizeof(PreCompute2DTexDesc));
	PreCompute2DTexDesc.Width = TRANSMITTANCE_TEXTURE_WIDTH;
	PreCompute2DTexDesc.Height = TRANSMITTANCE_TEXTURE_HEIGHT;
	PreCompute2DTexDesc.MipLevels = 1;
	PreCompute2DTexDesc.ArraySize = 1;
	PreCompute2DTexDesc.Format = DXGI_FORMAT_R16G16B16A16_FLOAT;
	PreCompute2DTexDesc.SampleDesc.Count = 1;
	PreCompute2DTexDesc.SampleDesc.Quality = 0;
	PreCompute2DTexDesc.Usage = D3D11_USAGE_DEFAULT;
	PreCompute2DTexDesc.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
	PreCompute2DTexDesc.CPUAccessFlags = 0;
	PreCompute2DTexDesc.MiscFlags = 0;
	pTransmittanceTex2D.Release();
	V_RETURN(pDevice->CreateTexture2D(&PreCompute2DTexDesc, nullptr, &pTransmittanceTex2D));

	D3D11_RENDER_TARGET_VIEW_DESC PreComputeRTVDesc;
	PreComputeRTVDesc.Format = PreCompute2DTexDesc.Format;
	PreComputeRTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
	PreComputeRTVDesc.Texture2D.MipSlice = 0;
	pTransmittanceRTV.Release();
	V_RETURN(pDevice->CreateRenderTargetView(pTransmittanceTex2D, &PreComputeRTVDesc, &pTransmittanceRTV));

	D3D11_SHADER_RESOURCE_VIEW_DESC PreComputeSRVDesc;
	PreComputeSRVDesc.Format = PreCompute2DTexDesc.Format;
	PreComputeSRVDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
	PreComputeSRVDesc.Texture2D.MostDetailedMip = 0;
	PreComputeSRVDesc.Texture2D.MipLevels = 1;
	pTransmittanceSRV.Release();
	V_RETURN(pDevice->CreateShaderResourceView(pTransmittanceTex2D, &PreComputeSRVDesc, &pTransmittanceSRV));

	ID3DX11EffectTechnique* activeTech = AtmosphereTechMap["ComputeTransmittanceTex2DTech"];

	RenderQuad(pContext, activeTech, &pTransmittanceRTV.p, TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
	
	//V_RETURN(D3DX11SaveTextureToFile(pContext, pTransmittanceTex2D, D3DX11_IFF_DDS, L"transmittance.dds"));
	return hr;
}

HRESULT Atmosphere::PreComputeSingleSctrTex3D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;
	return hr;
}