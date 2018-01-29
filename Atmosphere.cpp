#include "DXUT.h"
#include "Atmosphere.h"
#include "Common.h"
#include <minwinbase.h>
#include <sstream>

Atmosphere::Atmosphere()
{
}


Atmosphere::~Atmosphere()
{
}

void Atmosphere::Initialize()
{
	constexpr double kPi = 3.1415926;
	constexpr double kSunAngularRadius = 0.00935 / 2.0;
	constexpr double kSunSolidAngle = kPi * kSunAngularRadius * kSunAngularRadius;

	std::vector<double> lambda =
	{
		680.0,//R
		550.0,//G
		440.0 //B
	};

	constexpr int lambdaMin = 360;//nm
	constexpr int lambdaMax = 830;//nm

	constexpr double kSolarIrradiance[48] = {
		1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
		1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
		1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
		1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
		1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
		1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
	};

	constexpr double kOzoneCrossSection[48] = {
		1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
		8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26, 9.016e-26,
		1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25, 4.266e-25,
		4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25, 3.74e-25, 3.215e-25,
		2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25, 1.209e-25, 9.423e-26, 7.455e-26,
		6.566e-26, 5.105e-26, 4.15e-26, 4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26,
		2.534e-26, 1.624e-26, 1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
	};
	// From https://en.wikipedia.org/wiki/Dobson_unit, in molecules.m^-2.
	constexpr double kDobsonUnit = 2.687e20;
	// Maximum number density of ozone molecules, in m^-3 (computed so at to get
	// 300 Dobson units of ozone - for this we divide 300 DU by the integral of
	// the ozone density profile defined below, which is equal to 15km).
	constexpr double kMaxOzoneNumberDensity = 300.0 * kDobsonUnit / 15000.0;

	constexpr double kRayleigh = 1.24062e-6;
	constexpr double kRayleighScaleHeight = 8000.0;
	constexpr double kMieScaleHeight = 1200.0;
	constexpr double kMieAngstromAlpha = 0.0;
	constexpr double kMieAngstromBeta = 5.328e-3;
	constexpr double kMieSingleScatteringAlbedo = 0.9;
	constexpr double kGroundAlbedo = 0.1;

	//std::vector<float> solar_irradiance;
	//std::vector<float> rayleigh_scatter;
	//std::vector<float> mie_scatter;
	//std::vector<float> mie_extinction;
	//std::vector<float> absorption_extinction;	//ozone

	auto lerp = [](double x,double y,double s)
	{
		return static_cast<float>(x*(1 - s) + y*s);
	};

	ZeroMemory(&atmosphereParams, sizeof(AtmosphereParameters));
	for(int i=0;i<lambda.size();++i)
	{
		double l = lambda[i];
		int index = (l - lambdaMin) / 10;
		int lambda_x = index * 10 + lambdaMin;
		int lambda_y = lambda_x == lambdaMax ?
							lambdaMax : lambda_x + 10;
		double lambda_x_um = static_cast<double>(lambda_x) * 1e-3;
		double lambda_y_um = static_cast<double>(lambda_y) * 1e-3;

		double s = (l - static_cast<double>(lambda_x)) / 10;
		
		double mie = lerp(kMieAngstromBeta / kMieScaleHeight * pow(lambda_x_um, -kMieAngstromAlpha),
					kMieAngstromBeta / kMieScaleHeight * pow(lambda_y_um, -kMieAngstromAlpha), s)*1000;

		atmosphereParams.solar_irradiance[i]=
			lerp(kSolarIrradiance[index], kSolarIrradiance[index + 1], s);
		atmosphereParams.rayleigh_scattering[i] =
			lerp(kRayleigh * pow(lambda_x_um, -4), kRayleigh *pow(lambda_y_um, -4), s)*1000;
		atmosphereParams.mie_scattering[i] = mie;
		atmosphereParams.mie_extinction[i] = mie*kMieSingleScatteringAlbedo;
		atmosphereParams.absorption_extinction[i] =
			lerp(kMaxOzoneNumberDensity * kOzoneCrossSection[index],
				kMaxOzoneNumberDensity * kOzoneCrossSection[index+1],s) * 1000;
	}
	
	atmosphereParams.mu_s_min = cos(102.0*D3DX_PI/180.0);
	atmosphereParams.sun_angular_radius = 0.2678 * D3DX_PI / 180.0;
	atmosphereParams.bottom_radius = 6360.f;
	atmosphereParams.top_radius = 6420.f;
	atmosphereParams.mie_g = 0.8f;
	atmosphereParams.ground_albedo = 0.1f;
	atmosphereParams.ozone_width = 25.0f;
	atmosphereParams.rayleigh_density = DensityProfileLayer
	{
		1.f, -1.0 / kMieScaleHeight * 1000.0, 0.f, 0.f
	};
	atmosphereParams.mie_density = DensityProfileLayer
	{
		1.f, -1.0 / kRayleighScaleHeight * 1000.0, 0.f, 0.f
	};
	atmosphereParams.ozone_density[0] = DensityProfileLayer
	{
		0.0, 0.0, 1.0 / 15.0, -2.0 / 3.0
	};
	atmosphereParams.ozone_density[1] = DensityProfileLayer
	{
		0.0, 0.0, -1.0 / 15.0, 8.0 / 3.0
	};
	
	scatter_order_num = 3;
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

	pSingleScatterTex3D.Release();
	pSingleScatterSRV.Release();

	pSingleScatterRayleighTex3D.Release();
	pSingleScatterRayleighSRV.Release();

	pSingleScatterMieTex3D.Release();
	pSingleScatterMieSRV.Release();

	pDirectIrradianceTex2D.Release();
	pDirectIrradianceSRV.Release();

	pIndirectIrradianceTex2D.Release();
	pIndirectIrradianceSRV.Release();

	pMultiScatterTex3D.Release();
	pMultiScatterSRV.Release();
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


void Atmosphere::Render(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11RenderTargetView* pRTV)
{
	if(!IsPreComputed)
	{
		PreComputeTransmittanceTex2D(pDevice, pContext);

		PreComputeDirectIrradianceTex2D(pDevice, pContext);

		//PreComputeSingleSctrTex3D_Test(pDevice, pContext);
		PreComputeSingleSctrTex3D(pDevice, pContext);

		for (int scatter_order = 2; scatter_order <= scatter_order_num; ++scatter_order)
		{
			PreComputeInDirectIrradianceTex2D(pDevice, pContext, scatter_order - 1);
			PreComputeMultiSctrTex3D_Test(pDevice, pContext, scatter_order);
			if(scatter_order!=scatter_order_num)
				PreComputeMultiSctrTex3D(pDevice, pContext, scatter_order);
		}
		IsPreComputed = true;
	}
	//PreComputeTransmittanceTex2D(pDevice, pContext);

	//PreComputeDirectIrradianceTex2D(pDevice, pContext);

	////PreComputeSingleSctrTex3D_Test(pDevice, pContext);
	//PreComputeSingleSctrTex3D(pDevice, pContext);

	//for(int scatter_order = 2;scatter_order<=scatter_order_num;++scatter_order)
	//{
	//	PreComputeInDirectIrradianceTex2D(pDevice, pContext, scatter_order-1);
	//	PreComputeMultiSctrTex3D_Test(pDevice, pContext, scatter_order - 1);
	//}
}

HRESULT Atmosphere::PreComputeTransmittanceTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	pTransmittanceTex2D.Release();
	pTransmittanceSRV.Release();
	CComPtr<ID3D11RenderTargetView>	pTransmittanceRTV;
	V_RETURN(CreateTexture2D(pDevice, pContext, TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT, format,
							&pTransmittanceTex2D.p, &pTransmittanceSRV.p, &pTransmittanceRTV.p));

	ID3DX11EffectTechnique* activeTech = AtmosphereTechMap["ComputeTransmittanceTex2DTech"];

	VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

	ID3D11RenderTargetView* pRTVs[] =
	{
		pTransmittanceRTV.p,
	};
	pContext->OMSetRenderTargets(ARRAYSIZE(pRTVs), pRTVs, nullptr);

	RenderQuad(pContext, activeTech, TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
	//pContext->OMSetRenderTargets(0, nullptr, nullptr);
	V_RETURN(D3DX11SaveTextureToFile(pContext, pTransmittanceTex2D, D3DX11_IFF_DDS, L"Texture/Transmittance.dds"));
	return hr;
}

HRESULT Atmosphere::PreComputeDirectIrradianceTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	pDirectIrradianceTex2D.Release();
	pDirectIrradianceSRV.Release();
	CComPtr<ID3D11RenderTargetView>	pDirectIrradianceRTV;
	V_RETURN(CreateTexture2D(pDevice, pContext, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, format,
							&pDirectIrradianceTex2D.p, &pDirectIrradianceSRV.p, &pDirectIrradianceRTV.p));

	ID3DX11EffectTechnique* activeTech = AtmosphereTechMap["ComputeDirectIrradiance2DTech"];

	ID3D11RenderTargetView* pRTVs[] =
	{
		pDirectIrradianceRTV.p,
	};
	pContext->OMSetRenderTargets(ARRAYSIZE(pRTVs), pRTVs, nullptr);

	ShaderResourceVarMap["g_tex2DTransmittanceLUT"]->SetResource(pTransmittanceSRV);

	RenderQuad(pContext, activeTech, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);

	V_RETURN(D3DX11SaveTextureToFile(pContext, pDirectIrradianceTex2D, D3DX11_IFF_DDS, L"Texture/DirectIrradiance.dds"));
	return hr;
}

HRESULT Atmosphere::PreComputeSingleSctrTex3D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	pSingleScatterTex3D.Release();
	pSingleScatterRayleighTex3D.Release();
	pSingleScatterMieTex3D.Release();
	pSingleScatterSRV.Release();
	pSingleScatterRayleighSRV.Release();
	pSingleScatterMieSRV.Release();
	V_RETURN(CreateTexture3D(pDevice, pContext, 
							SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH,format,
							{ &pSingleScatterTex3D.p ,&pSingleScatterRayleighTex3D.p ,&pSingleScatterMieTex3D.p },
							{ &pSingleScatterSRV,&pSingleScatterRayleighSRV.p,&pSingleScatterMieSRV.p }));
	
	ShaderResourceVarMap["g_tex2DTransmittanceLUT"]->SetResource(pTransmittanceSRV);

	std::vector<CComPtr<ID3D11RenderTargetView>> pSingleScaterRTVs(SCATTERING_TEXTURE_DEPTH);
	std::vector<CComPtr<ID3D11RenderTargetView>> pSingleScaterRayleighRTVs(SCATTERING_TEXTURE_DEPTH);
	std::vector<CComPtr<ID3D11RenderTargetView>> pSingleScaterMieRTVs(SCATTERING_TEXTURE_DEPTH);
	for(UINT depthSlice = 0;depthSlice<SCATTERING_TEXTURE_DEPTH;++depthSlice)
	{
		D3D11_RENDER_TARGET_VIEW_DESC CurrSliceRTVDesc;
		CurrSliceRTVDesc.Format = format;
		CurrSliceRTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
		CurrSliceRTVDesc.Texture3D.MipSlice = 0;
		CurrSliceRTVDesc.Texture3D.FirstWSlice = depthSlice;
		CurrSliceRTVDesc.Texture3D.WSize = 1;

		V_RETURN(pDevice->CreateRenderTargetView(pSingleScatterTex3D, &CurrSliceRTVDesc, &pSingleScaterRTVs[depthSlice]));
		V_RETURN(pDevice->CreateRenderTargetView(pSingleScatterRayleighTex3D, &CurrSliceRTVDesc, &pSingleScaterRayleighRTVs[depthSlice]));
		V_RETURN(pDevice->CreateRenderTargetView(pSingleScatterMieTex3D, &CurrSliceRTVDesc, &pSingleScaterMieRTVs[depthSlice]));

		ID3DX11EffectTechnique* activeTech = AtmosphereTechMap["ComputeSingleScatterTex3DTech"];

		VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

		MiscDynamicParams misc;
		UINT uiW = depthSlice % SCATTERING_TEXTURE_MU_S_SIZE;
		UINT uiQ = depthSlice / SCATTERING_TEXTURE_MU_S_SIZE;
		misc.f2WQ.x = ((float)uiW + 0.5f) / SCATTERING_TEXTURE_MU_S_SIZE;
		assert(0 < misc.f2WQ.x && misc.f2WQ.x < 1);
		misc.f2WQ.y = ((float)uiQ + 0.5f) / SCATTERING_TEXTURE_NU_SIZE;
		assert(0 < misc.f2WQ.y && misc.f2WQ.y < 1);
		VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));
		
		ID3D11RenderTargetView* pRTVs[] =
		{
			pSingleScaterRTVs[depthSlice].p,
			pSingleScaterRayleighRTVs[depthSlice].p,
			pSingleScaterMieRTVs[depthSlice].p
		};

		UINT size = ARRAYSIZE(pRTVs);
		for(UINT i = 0;i<size;++i)
		{
			float zero[] = { 0.f, 0.f,0.f,0.f };
			pContext->ClearRenderTargetView(pRTVs[i], zero);
		}
		pContext->OMSetRenderTargets(size, pRTVs, nullptr);

		RenderQuad(pContext, activeTech, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
	}
	V_RETURN(D3DX11SaveTextureToFile(pContext, pSingleScatterTex3D, D3DX11_IFF_DDS, L"Texture/SingleScatter3D.dds"));

	return hr;
}

HRESULT Atmosphere::PreComputeInDirectIrradianceTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, int scatter_order)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	CComPtr<ID3D11RenderTargetView>	pIndirectIrradianceRTV;
	V_RETURN(CreateTexture2D(pDevice, pContext, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT,
		format,&pIndirectIrradianceTex2D.p, &pIndirectIrradianceSRV.p, &pIndirectIrradianceRTV.p));

	ID3DX11EffectTechnique* activeTech = AtmosphereTechMap["ComputeIndirectIrradiance2DTech"];

	ID3D11RenderTargetView* pRTVs[] =
	{
		pIndirectIrradianceRTV.p,
	};
	pContext->OMSetRenderTargets(ARRAYSIZE(pRTVs), pRTVs, nullptr);

	if(1==scatter_order)
	{
		ShaderResourceVarMap["g_tex3DSingleScatteringLUT"]->SetResource(pSingleScatterSRV);
	}
	else
	{
		ShaderResourceVarMap["g_tex3DMultiScatteringLUT"]->SetResource(pMultiScatterSRV);
	}
	MiscDynamicParams misc;
	misc.scatter_order = scatter_order;
	VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));

	RenderQuad(pContext, activeTech, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);

	std::wstringstream ss;
	ss << "Texture/Irradiance/IndirectIrradiance_" << scatter_order << ".dds";
	V_RETURN(D3DX11SaveTextureToFile(pContext, pIndirectIrradianceTex2D, D3DX11_IFF_DDS, ss.str().c_str()));
	return hr;
}

HRESULT Atmosphere::PreComputeMultiSctrTex3D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext,int scatter_order)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	if(1==scatter_order-1)
	{
		ShaderResourceVarMap["g_tex3DSingleScatteringLUT"]->SetResource(pSingleScatterSRV);
	}
	else
	{
		ShaderResourceVarMap["g_tex3DMultiScatteringLUT"]->SetResource(pMultiScatterSRV);
	}
	ShaderResourceVarMap["g_tex2DTransmittanceLUT"]->SetResource(pTransmittanceSRV);

	CComPtr<ID3D11Texture3D> pTex3D;
	CComPtr<ID3D11ShaderResourceView> pSRV;
	V_RETURN(CreateTexture3D(pDevice, pContext, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH,
							format,{ &pTex3D.p },{ &pSRV.p }));

	std::vector<CComPtr<ID3D11RenderTargetView>> pMultiScaterRTVs(SCATTERING_TEXTURE_DEPTH);
	for (UINT depthSlice = 0; depthSlice<SCATTERING_TEXTURE_DEPTH; ++depthSlice)
	{
		D3D11_RENDER_TARGET_VIEW_DESC CurrSliceRTVDesc;
		CurrSliceRTVDesc.Format = format;
		CurrSliceRTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
		CurrSliceRTVDesc.Texture3D.MipSlice = 0;
		CurrSliceRTVDesc.Texture3D.FirstWSlice = depthSlice;
		CurrSliceRTVDesc.Texture3D.WSize = 1;

		V_RETURN(pDevice->CreateRenderTargetView(pTex3D, &CurrSliceRTVDesc, &pMultiScaterRTVs[depthSlice]));
		
		ID3DX11EffectTechnique* activeTech = AtmosphereTechMap["ComputeMultiScatterTex3DTech"];

		VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

		MiscDynamicParams misc;
		misc.scatter_order = scatter_order;
		UINT uiW = depthSlice % SCATTERING_TEXTURE_MU_S_SIZE;
		UINT uiQ = depthSlice / SCATTERING_TEXTURE_MU_S_SIZE;
		misc.f2WQ.x = ((float)uiW + 0.5f) / SCATTERING_TEXTURE_MU_S_SIZE;
		assert(0 < misc.f2WQ.x && misc.f2WQ.x < 1);
		misc.f2WQ.y = ((float)uiQ + 0.5f) / SCATTERING_TEXTURE_NU_SIZE;
		assert(0 < misc.f2WQ.y && misc.f2WQ.y < 1);
		VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));

		ID3D11RenderTargetView* pRTVs[] =
		{
			pMultiScaterRTVs[depthSlice].p,
		};

		UINT size = ARRAYSIZE(pRTVs);
		for (UINT i = 0; i<size; ++i)
		{
			float zero[] = { 0.f, 0.f,0.f,0.f };
			pContext->ClearRenderTargetView(pRTVs[i], zero);
		}
		pContext->OMSetRenderTargets(size, pRTVs, nullptr);

		RenderQuad(pContext, activeTech, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
		//pContext->OMSetRenderTargets(0, nullptr, nullptr);
	}

	std::wstringstream ss;
	ss << "Texture/MultiScatter/MultiScatter_" << scatter_order << ".dds";
	V_RETURN(D3DX11SaveTextureToFile(pContext, pTex3D, D3DX11_IFF_DDS, ss.str().c_str()));

	pMultiScatterTex3D.Release();
	pMultiScatterSRV.Release();

	pMultiScatterTex3D = pTex3D;
	pMultiScatterSRV = pSRV;

	return hr;
}

HRESULT Atmosphere::PreComputeSingleSctrTex3D_Test(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;
	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	//D3D11_TEXTURE3D_DESC PreCompute3DTexDesc;
	D3D11_TEXTURE2D_DESC PreCompute2DTexDesc;
	ZeroMemory(&PreCompute2DTexDesc, sizeof(PreCompute2DTexDesc));
	PreCompute2DTexDesc.Width = SCATTERING_TEXTURE_WIDTH;
	PreCompute2DTexDesc.Height = SCATTERING_TEXTURE_HEIGHT;
	//PreCompute3DTexDesc.Depth = SCATTERING_TEXTURE_DEPTH;
	PreCompute2DTexDesc.MipLevels = 1;
	PreCompute2DTexDesc.ArraySize = 1;//DEBUG 2D
	PreCompute2DTexDesc.SampleDesc.Count = 1;//DEBUG 2D
	PreCompute2DTexDesc.SampleDesc.Quality = 0;//DEBUG 2D
	PreCompute2DTexDesc.Format = format;
	PreCompute2DTexDesc.Usage = D3D11_USAGE_DEFAULT;
	PreCompute2DTexDesc.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
	PreCompute2DTexDesc.CPUAccessFlags = 0;
	PreCompute2DTexDesc.MiscFlags = 0;

	CComPtr<ID3D11Texture2D> pSingleScatterTex2D;
	CComPtr<ID3D11Texture2D> pSingleScatterRayleighTex2D;
	CComPtr<ID3D11Texture2D> pSingleScatterMieTex2D;

	V_RETURN(pDevice->CreateTexture2D(&PreCompute2DTexDesc, nullptr, &pSingleScatterTex2D));
	V_RETURN(pDevice->CreateTexture2D(&PreCompute2DTexDesc, nullptr, &pSingleScatterRayleighTex2D));
	V_RETURN(pDevice->CreateTexture2D(&PreCompute2DTexDesc, nullptr, &pSingleScatterMieTex2D));

	pSingleScatterSRV.Release();
	pSingleScatterRayleighSRV.Release();
	pSingleScatterMieSRV.Release();
	V_RETURN(pDevice->CreateShaderResourceView(pSingleScatterTex2D, nullptr, &pSingleScatterSRV));
	V_RETURN(pDevice->CreateShaderResourceView(pSingleScatterTex2D, nullptr, &pSingleScatterRayleighSRV));
	V_RETURN(pDevice->CreateShaderResourceView(pSingleScatterTex2D, nullptr, &pSingleScatterMieSRV));

	ShaderResourceVarMap["g_tex2DTransmittanceLUT"]->SetResource(pTransmittanceSRV);

	std::vector<CComPtr<ID3D11RenderTargetView>> pSingleScatterRTVs(SCATTERING_TEXTURE_DEPTH);
	std::vector<CComPtr<ID3D11RenderTargetView>> pSingleScatterRayleighRTVs(SCATTERING_TEXTURE_DEPTH);
	std::vector<CComPtr<ID3D11RenderTargetView>> pSingleScatterMieRTVs(SCATTERING_TEXTURE_DEPTH);
	for (UINT depthSlice = 0; depthSlice < SCATTERING_TEXTURE_DEPTH; ++depthSlice)
	{
		V_RETURN(pDevice->CreateRenderTargetView(pSingleScatterTex2D, nullptr, &pSingleScatterRTVs[depthSlice]));
		V_RETURN(pDevice->CreateRenderTargetView(pSingleScatterRayleighTex2D, nullptr, &pSingleScatterRayleighRTVs[depthSlice]));
		V_RETURN(pDevice->CreateRenderTargetView(pSingleScatterMieTex2D, nullptr, &pSingleScatterMieRTVs[depthSlice]));

		ID3DX11EffectTechnique* activeTech = AtmosphereTechMap["ComputeSingleScatterTex3DTech"];

		VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

		MiscDynamicParams misc;
		UINT uiW = depthSlice % SCATTERING_TEXTURE_MU_S_SIZE;
		UINT uiQ = depthSlice / SCATTERING_TEXTURE_MU_S_SIZE;
		misc.f2WQ.x = ((float)uiW + 0.5f) / SCATTERING_TEXTURE_MU_S_SIZE;
		assert(0 < misc.f2WQ.x && misc.f2WQ.x < 1);
		misc.f2WQ.y = ((float)uiQ + 0.5f) / SCATTERING_TEXTURE_NU_SIZE;
		assert(0 < misc.f2WQ.y && misc.f2WQ.y < 1);
		VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));

		ID3D11RenderTargetView* pRTVs[] =
		{
			pSingleScatterRTVs[depthSlice].p,
			pSingleScatterRayleighRTVs[depthSlice].p,
			pSingleScatterMieRTVs[depthSlice].p
		};

		UINT size = ARRAYSIZE(pRTVs);
		for (UINT i = 0; i < size; ++i)
		{
			float zero[] = { 0.f, 0.f,0.f,0.f };
			pContext->ClearRenderTargetView(pRTVs[i], zero);
		}
		pContext->OMSetRenderTargets(size, pRTVs, nullptr);

		RenderQuad(pContext, activeTech, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
		//pContext->OMSetRenderTargets(0, nullptr, nullptr);

		std::wstringstream ss;
		ss << L"Texture/SingleScatter/SingleScatter_" << depthSlice<<".dds";

		V_RETURN(D3DX11SaveTextureToFile(pContext, pSingleScatterTex2D, D3DX11_IFF_DDS, ss.str().c_str()));
	}
	return hr;
}

HRESULT Atmosphere::PreComputeMultiSctrTex3D_Test(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, int scatter_order)
{
	HRESULT hr = S_OK;
	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	if (1 == scatter_order - 1)
	{
		ShaderResourceVarMap["g_tex3DSingleScatteringLUT"]->SetResource(pSingleScatterSRV);
	}
	else
	{
		ShaderResourceVarMap["g_tex3DMultiScatteringLUT"]->SetResource(pMultiScatterSRV);
	}
	ShaderResourceVarMap["g_tex2DTransmittanceLUT"]->SetResource(pTransmittanceSRV);

	D3D11_TEXTURE2D_DESC PreCompute2DTexDesc;
	ZeroMemory(&PreCompute2DTexDesc, sizeof(PreCompute2DTexDesc));
	PreCompute2DTexDesc.Width = SCATTERING_TEXTURE_WIDTH;
	PreCompute2DTexDesc.Height = SCATTERING_TEXTURE_HEIGHT;
	//PreCompute3DTexDesc.Depth = SCATTERING_TEXTURE_DEPTH;
	PreCompute2DTexDesc.MipLevels = 1;
	PreCompute2DTexDesc.ArraySize = 1;//DEBUG 2D
	PreCompute2DTexDesc.SampleDesc.Count = 1;//DEBUG 2D
	PreCompute2DTexDesc.SampleDesc.Quality = 0;//DEBUG 2D
	PreCompute2DTexDesc.Format = format;
	PreCompute2DTexDesc.Usage = D3D11_USAGE_DEFAULT;
	PreCompute2DTexDesc.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
	PreCompute2DTexDesc.CPUAccessFlags = 0;
	PreCompute2DTexDesc.MiscFlags = 0;

	CComPtr<ID3D11Texture2D> pMultiScatterTex2D;

	V_RETURN(pDevice->CreateTexture2D(&PreCompute2DTexDesc, nullptr, &pMultiScatterTex2D));

	std::vector<CComPtr<ID3D11RenderTargetView>> pMultiScatterRTVs(SCATTERING_TEXTURE_DEPTH);
	std::vector<CComPtr<ID3D11RenderTargetView>> pMultiScatterSRVs(SCATTERING_TEXTURE_DEPTH);
	for (UINT depthSlice = 0; depthSlice<SCATTERING_TEXTURE_DEPTH; ++depthSlice)
	{
		V_RETURN(pDevice->CreateRenderTargetView(pMultiScatterTex2D, nullptr, &pMultiScatterSRVs[depthSlice]));
		V_RETURN(pDevice->CreateRenderTargetView(pMultiScatterTex2D, nullptr, &pMultiScatterRTVs[depthSlice]));

		ID3DX11EffectTechnique* activeTech = AtmosphereTechMap["ComputeMultiScatterTex3DTech"];

		VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

		MiscDynamicParams misc;
		misc.scatter_order = scatter_order;
		UINT uiW = depthSlice % SCATTERING_TEXTURE_MU_S_SIZE;
		UINT uiQ = depthSlice / SCATTERING_TEXTURE_MU_S_SIZE;
		misc.f2WQ.x = ((float)uiW + 0.5f) / SCATTERING_TEXTURE_MU_S_SIZE;
		assert(0 < misc.f2WQ.x && misc.f2WQ.x < 1);
		misc.f2WQ.y = ((float)uiQ + 0.5f) / SCATTERING_TEXTURE_NU_SIZE;
		assert(0 < misc.f2WQ.y && misc.f2WQ.y < 1);
		VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));

		ID3D11RenderTargetView* pRTVs[] =
		{
			pMultiScatterRTVs[depthSlice].p,
		};

		UINT size = ARRAYSIZE(pRTVs);
		for (UINT i = 0; i<size; ++i)
		{
			float zero[] = { 0.f, 0.f,0.f,0.f };
			pContext->ClearRenderTargetView(pRTVs[i], zero);
		}
		pContext->OMSetRenderTargets(size, pRTVs, nullptr);

		RenderQuad(pContext, activeTech, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
		//pContext->OMSetRenderTargets(0, nullptr, nullptr);
		std::wstringstream ss;
		ss << "Texture/MultiScatter/scatter_order_" << scatter_order<<"/"<< "MultiScatter_" <<depthSlice << ".dds";
		V_RETURN(D3DX11SaveTextureToFile(pContext, pMultiScatterTex2D, D3DX11_IFF_DDS, ss.str().c_str()));
	}

	return hr;
}