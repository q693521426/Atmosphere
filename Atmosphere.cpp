#include "DXUT.h"
#include "Atmosphere.h"
#include <minwinbase.h>
#include <sstream>
#include <fstream>

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
			kMieAngstromBeta / kMieScaleHeight * pow(lambda_y_um, -kMieAngstromAlpha), s) * 1000;

		atmosphereParams.solar_irradiance[i]=
			lerp(kSolarIrradiance[index], kSolarIrradiance[index + 1], s) ;
		atmosphereParams.rayleigh_scattering[i] =
			lerp(kRayleigh * pow(lambda_x_um, -4), kRayleigh *pow(lambda_y_um, -4), s) * 1000;
		atmosphereParams.mie_scattering[i] = mie;
		atmosphereParams.mie_extinction[i] = mie*kMieSingleScatteringAlbedo;
		atmosphereParams.absorption_extinction[i] =
			lerp(kMaxOzoneNumberDensity * kOzoneCrossSection[index],
				kMaxOzoneNumberDensity * kOzoneCrossSection[index+1],s) * 1000;
	}
	
	atmosphereParams.mu_s_min = cos(102.0*D3DX_PI/180.0);
	atmosphereParams.sun_angular_radius = 32.f / 2.f / 60.f * ((2.f * D3DX_PI) / 180); //0.2678 * D3DX_PI / 180.0;
	atmosphereParams.bottom_radius = 6360.f;
	atmosphereParams.top_radius = 6420.f;
	atmosphereParams.mie_g = 0.8f;
	atmosphereParams.ground_albedo = 0.1f;
	atmosphereParams.ozone_width = 25.0f;
	atmosphereParams.nu_power = 1.5;
	atmosphereParams.exposure = exposure;
	atmosphereParams.rayleigh_density = DensityProfileLayer
	{
		1.f, -1.0 / kRayleighScaleHeight * 1000.0, 0.f, 0.f
	};
	atmosphereParams.mie_density = DensityProfileLayer
	{
		1.f, -1.0 / kMieScaleHeight * 1000.0, 0.f, 0.f
	};
	atmosphereParams.ozone_density[0] = DensityProfileLayer
	{
		0.0, 0.0, 1.0 / 15.0, -2.0 / 3.0
	};
	atmosphereParams.ozone_density[1] = DensityProfileLayer
	{
		0.0, 0.0, -1.0 / 15.0, 8.0 / 3.0
	};
	
	scatter_order_num = 4;
	
	m_pCloud = new Cloud();

#if CREATE_TEXTURE_DDS_TEST

	//V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"Texture/Transmittance.dds", nullptr, nullptr, &pTransmittanceSRV.p, nullptr));
	//V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"Texture/", nullptr, nullptr, &pSingleScatterCombinedSRV.p, nullptr));
	//V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"", nullptr, nullptr, &pSingleScatterMieSRV.p, nullptr));
	//V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"", nullptr, nullptr, &pDirectIrradianceSRV.p, nullptr));
	//V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"", nullptr, nullptr, &pMultiScatterCombinedSRV.p, nullptr));

	CreateDirectory(L"Texture", nullptr);
	CreateDirectory(L"Texture/SingleScatter", nullptr);
	CreateDirectory(L"Texture/MultiScatter", nullptr);
	CreateDirectory(L"Texture/Irradiance", nullptr);
	for(UINT order =2;order<= scatter_order_num;order++)
	{
		std::wstringstream ss;
		ss << "Texture/MultiScatter/scatter_order_" << order;
		CreateDirectory(ss.str().c_str(), nullptr);
	}
#endif
	SetView(9000.0, 1.47, -0.1, 1.3, 2.9, 10.0);
}	


void Atmosphere::Release()
{
	pTransmittanceTex2D.Release();
	pTransmittanceSRV.Release();

	pSingleScatterTex3D.Release();
	pSingleScatterSRV.Release();

	pSingleScatterCombinedTex3D.Release();
	pSingleScatterCombinedSRV.Release();

	pSingleScatterMieTex3D.Release();
	pSingleScatterMieSRV.Release();

	pDirectIrradianceTex2D.Release();
	pDirectIrradianceSRV.Release();

	pIndirectIrradianceTex2D.Release();
	pIndirectIrradianceSRV.Release();

	pMultiScatterTex3D.Release();
	pMultiScatterSRV.Release();

	pMultiScatterCombinedTex3D.Release();
	pMultiScatterCombinedSRV.Release();
	if(m_pCloud)
	{
		m_pCloud->Release();
		delete m_pCloud;
		m_pCloud = nullptr;
	}

	pSpaceLinearDepthTex2D.Release();
	pSpaceLinearDepthSRV.Release();

	pSliceEndTex2D.Release();
	pSliceEndSRV.Release();

	pEpipolarSampleTex2D.Release();
	pEpipolarSampleSRV.Release();
	pEpipolarSampleDSV.Release();

	pEpipolarSampleCamDepthTex2D.Release();
	pEpipolarSampleCamDepthSRV.Release();

	pUnshadowedSampleScatterTex2D.Release();
	pUnshadowedSampleScatterSRV.Release();

	pInterpolationSampleTex2D.Release();
	pInterpolationSampleSRV.Release();
	pInterpolationSampleUAV.Release();

	pSliceUVOrigDirTex2D.Release();
	pSliceUVOrigDirSRV.Release();

	for (int i = 0; i < 2; i++)
	{
		pMinMaxMinMapTex2D[i].Release();
		pMinMaxMinMapTexSRV[i].Release();
	}

	pSampleScatterTex2D.Release();
	pSampleScatterSRV.Release();

	pInterpolatedSampleScatterTex2D.Release();
	pInterpolatedSampleScatterSRV.Release();
	
	pApplyScatterDSV.Release();

	GameObject::Release();
}


HRESULT Atmosphere::OnD3D11CreateDevice(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	V_RETURN(GameObject::OnD3D11CreateDevice(pDevice, pContext, L"Atmosphere.fx", TechStr, VarStr, ShaderResourceVarStr));
	V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"Texture/earth.tiff", nullptr, nullptr, &pEarthGroundSRV.p, nullptr));

#if USE_LUT_DDS
	V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"Texture/OpticalLength.dds", nullptr, nullptr, &pOpticalLengthSRV.p, nullptr));
	V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"Texture/DirectIrradiance.dds", nullptr, nullptr, &pDirectIrradianceSRV.p, nullptr));
	V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"Texture/SingleScatterCombined.dds", nullptr, nullptr, &pSingleScatterCombinedSRV.p, nullptr));
	V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"Texture/IndirectIrradiance.dds", nullptr, nullptr, &pIndirectIrradianceSRV.p, nullptr));
	V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"Texture/MultiScatterCombined.dds", nullptr, nullptr, &pMultiScatterCombinedSRV.p, nullptr));

	IsPreComputed = true;
#endif

	SetTextureSize();

	m_pCloud->OnD3D11CreateDevice(pDevice, pContext);

	return hr;
}


void Atmosphere::SetTextureSize()
{
	VarMap["SCREEN_WIDTH"]->SetRawValue(&screen_width, 0, sizeof(UINT));
	VarMap["SCREEN_HEIGHT"]->SetRawValue(&screen_height, 0, sizeof(UINT));

	VarMap["TRANSMITTANCE_TEXTURE_WIDTH"]->SetRawValue(&TRANSMITTANCE_TEXTURE_WIDTH, 0, sizeof(UINT));
	VarMap["TRANSMITTANCE_TEXTURE_HEIGHT"]->SetRawValue(&TRANSMITTANCE_TEXTURE_WIDTH, 0, sizeof(UINT));

	VarMap["SCATTERING_TEXTURE_R_SIZE"]->SetRawValue(&SCATTERING_TEXTURE_R_SIZE, 0, sizeof(UINT));
	VarMap["SCATTERING_TEXTURE_MU_SIZE"]->SetRawValue(&SCATTERING_TEXTURE_MU_SIZE, 0, sizeof(UINT));
	VarMap["SCATTERING_TEXTURE_MU_S_SIZE"]->SetRawValue(&SCATTERING_TEXTURE_MU_S_SIZE, 0, sizeof(UINT));
	VarMap["SCATTERING_TEXTURE_NU_SIZE"]->SetRawValue(&SCATTERING_TEXTURE_NU_SIZE, 0, sizeof(UINT));

	VarMap["SCATTERING_TEXTURE_WIDTH"]->SetRawValue(&SCATTERING_TEXTURE_WIDTH, 0, sizeof(UINT));
	VarMap["SCATTERING_TEXTURE_HEIGHT"]->SetRawValue(&SCATTERING_TEXTURE_HEIGHT, 0, sizeof(UINT));
	VarMap["SCATTERING_TEXTURE_DEPTH"]->SetRawValue(&SCATTERING_TEXTURE_DEPTH, 0, sizeof(UINT));

	VarMap["IRRADIANCE_TEXTURE_WIDTH"]->SetRawValue(&IRRADIANCE_TEXTURE_WIDTH, 0, sizeof(UINT));
	VarMap["IRRADIANCE_TEXTURE_HEIGHT"]->SetRawValue(&IRRADIANCE_TEXTURE_HEIGHT, 0, sizeof(UINT));

	VarMap["EPIPOLAR_SLICE_NUM"]->SetRawValue(&EPIPOLAR_SLICE_NUM, 0, sizeof(UINT));
	VarMap["EPIPOLAR_SAMPLE_NUM"]->SetRawValue(&EPIPOLAR_SAMPLE_NUM, 0, sizeof(UINT));
}


void Atmosphere::SetCameraParams()
{
	//D3DXVECTOR3 EyePos = *m_FirstPersonCamera.GetEyePt();
	//D3DXVECTOR3 LookAt = *m_FirstPersonCamera.GetLookAtPt();

	cameraParams.f3CameraPos = f3CamPos;
	cameraParams.View = camView;
	cameraParams.Proj = camProj;
	cameraParams.ViewProj = cameraParams.View * cameraParams.Proj;
	float det = D3DXMatrixDeterminant(&cameraParams.ViewProj);
	D3DXMatrixInverse(&cameraParams.InvViewProj, &det, &cameraParams.ViewProj);

	D3DXMatrixTranspose(&cameraParams.View, &cameraParams.View);
	D3DXMatrixTranspose(&cameraParams.Proj, &cameraParams.Proj);
	D3DXMatrixTranspose(&cameraParams.ViewProj, &cameraParams.ViewProj);
	D3DXMatrixTranspose(&cameraParams.InvViewProj, &cameraParams.InvViewProj);

	cameraParams.f3CameraDir = f3CamDir;
	D3DXVec3Normalize(&cameraParams.f3CameraDir, &cameraParams.f3CameraDir);

	cameraParams.fNearZ = fCamNear;
	cameraParams.fFarZ = fCamFar * 0.999999f;

	VarMap["camera"]->SetRawValue(&cameraParams, 0, sizeof(CameraParams));
}


void Atmosphere::SetLightParams()
{
	lightParams.f3LightDir = GetSunDir();
	//D3DXMATRIX View = *m_FirstPersonCamera.GetViewMatrix();
	//D3DXMATRIX Proj = *m_FirstPersonCamera.GetProjMatrix();
	D3DXMATRIX camViewProj = camView * camProj;
	D3DXVECTOR4 f4LightPos = D3DXVECTOR4(lightParams.f3LightDir + f3CamPos,1);
	D3DXVec4Transform(&lightParams.f4LightScreenPos, &f4LightPos, &camViewProj);
	if (lightParams.f4LightScreenPos.w < 0)
	{
		fIsLightInSpaceCorrect = 0;
	}
	else
	{
		fIsLightInSpaceCorrect = 1;
	}

	lightParams.f4LightScreenPos /= abs(lightParams.f4LightScreenPos.w);
	float fDistLightToScreen = D3DXVec2Length((D3DXVECTOR2*)(&lightParams.f4LightScreenPos));
	//float fMaxDist = 100;
	//if (fDistLightToScreen > 100)
	//	lightParams.f4LightScreenPos *= fMaxDist / fDistLightToScreen;

	lightParams.View = lightView;
	lightParams.Proj = lightProj;
	lightParams.ViewProj = lightParams.View *lightParams.Proj;
	float det = D3DXMatrixDeterminant(&lightParams.ViewProj);
	D3DXMatrixInverse(&lightParams.InvViewProj, &det, &lightParams.ViewProj);

	D3DXMatrixTranspose(&lightParams.View, &lightParams.View);
	D3DXMatrixTranspose(&lightParams.Proj, &lightParams.Proj);
	D3DXMatrixTranspose(&lightParams.ViewProj, &lightParams.ViewProj);
	D3DXMatrixTranspose(&lightParams.InvViewProj, &lightParams.InvViewProj);

	VarMap["light"]->SetRawValue(&lightParams, 0, sizeof(LightParams));

	fEnableLightShaft = 1;
}


HRESULT Atmosphere::PreCompute(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11RenderTargetView* pRTV)
{
	HRESULT hr = S_OK;

	if(!IsPreComputed)
	{
		V_RETURN(PreComputeOpticalLengthTex2D(pDevice, pContext));
		//V_RETURN(PreComputeTransmittanceTex2D(pDevice, pContext));
		V_RETURN(PreComputeDirectIrradianceTex2D(pDevice, pContext));
		V_RETURN(PreComputeSingleSctrTex3D(pDevice, pContext));

		for (int scatter_order = 2; scatter_order <= scatter_order_num; ++scatter_order)
		{
			V_RETURN(PreComputeInDirectIrradianceTex2D(pDevice, pContext, scatter_order - 1));
			V_RETURN(PreComputeMultiSctrTex3D(pDevice, pContext, scatter_order));
		}
		//pContext->CopyResource(pMultiScatterTex3D, pSingleScatterTex3D);

		//m_pCloud->PreCompute(pDevice, pContext, pRTV);
		{
			V_RETURN(SaveTextureToDDS(pContext, "Texture/OpticalLength.dds", pOpticalLengthTex2D));
			V_RETURN(SaveTextureToDDS(pContext, "Texture/DirectIrradiance.dds", pDirectIrradianceTex2D));
			V_RETURN(SaveTextureToDDS(pContext, "Texture/SingleScatterCombined.dds", pSingleScatterCombinedTex3D));
			V_RETURN(SaveTextureToDDS(pContext, "Texture/IndirectIrradiance.dds", pIndirectIrradianceTex2D));
			V_RETURN(SaveTextureToDDS(pContext, "Texture/MultiScatterCombined.dds", pMultiScatterCombinedTex3D));
		}
		IsPreComputed = true;
	}
	return hr;
}


void Atmosphere::Render(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, 
						ID3D11RenderTargetView* pRTV,
						ID3D11ShaderResourceView *pColorBufferSRV,
						ID3D11ShaderResourceView* pDepthSRV,
						ID3D11ShaderResourceView* pShadowMapSRV,
						UINT shadowMapResolution)
{
	SetCameraParams();
	SetLightParams();
	VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));
	VarMap["SHADOWMAP_TEXTURE_DIM"]->SetRawValue(&shadowMapResolution, 0, sizeof(UINT));

	//RenderBackGround(pDevice, pContext, pRTV, pColorBufferSRV,pDepthSRV);
	ComputeSpaceLinearDepthTex2D(pDevice, pContext, pDepthSRV);
	ComputeSliceEndTex2D(pDevice, pContext);
	ComputeEpipolarCoordTex2D(pDevice, pContext);
	ComputeUnshadowedSampleScatter(pDevice, pContext);
	RefineSampleLocal(pDevice, pContext);
	ComputeSliceUVOrigDirTex2D(pDevice, pContext, pShadowMapSRV);
	Build1DMinMaxMipMap(pDevice, pContext, pShadowMapSRV, shadowMapResolution);
	MarkRayMarchSample(pDevice, pContext);
	if(fEnableLightShaft)
		DoRayMarch(pDevice, pContext, pShadowMapSRV);
	InterpolateScatter(pDevice, pContext);
	ApplyAndFixInterpolateScatter(pDevice, pContext, pRTV, pColorBufferSRV);
}


void Atmosphere::RenderBackGround(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, 
								ID3D11RenderTargetView* pRTV,
								ID3D11ShaderResourceView* pColorBufferSRV,
								ID3D11ShaderResourceView* pDepthSRV)
{
	ID3DX11EffectTechnique* activeTech = TechMap["DrawGroundAndSkyTech"];
	MiscDynamicParams misc;
	misc.scatter_order = 2;
	VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));
	
	ShaderResourceVarMap["g_tex2DTransmittanceLUT"]->SetResource(pTransmittanceSRV);
	ShaderResourceVarMap["g_tex2DOpticalLengthLUT"]->SetResource(pOpticalLengthSRV);
	ShaderResourceVarMap["g_tex2DDirectIrradianceLUT"]->SetResource(pDirectIrradianceSRV);
	ShaderResourceVarMap["g_tex2DIndirectIrradianceLUT"]->SetResource(pIndirectIrradianceSRV);
	ShaderResourceVarMap["g_tex3DSingleScatteringLUT"]->SetResource(pSingleScatterSRV);
	ShaderResourceVarMap["g_tex3DMultiScatteringLUT"]->SetResource(pMultiScatterSRV);
	ShaderResourceVarMap["g_tex3DSingleMieScatteringLUT"]->SetResource(pSingleScatterMieSRV);
	ShaderResourceVarMap["g_tex3DSingleScatteringCombinedLUT"]->SetResource(pSingleScatterCombinedSRV);
	ShaderResourceVarMap["g_tex3DMultiScatteringCombinedLUT"]->SetResource(pMultiScatterCombinedSRV);
	ShaderResourceVarMap["g_tex2DEarthGround"]->SetResource(pEarthGroundSRV);
	ShaderResourceVarMap["g_tex2DSpaceDepth"]->SetResource(pDepthSRV);
	ShaderResourceVarMap["g_tex2DColorBuffer"]->SetResource(pColorBufferSRV);

	pContext->OMSetRenderTargets(1, &pRTV, nullptr);
	RenderQuad(pContext, activeTech, screen_width, screen_height);

	UnbindResources(pContext);
}


HRESULT Atmosphere::PreComputeTransmittanceTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	pTransmittanceTex2D.Release();
	pTransmittanceSRV.Release();
	CComPtr<ID3D11RenderTargetView>	pTransmittanceRTV;
	V_RETURN(CreateTexture2D(pDevice, pContext, TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT, format,
	{ &pTransmittanceTex2D.p }, { &pTransmittanceSRV.p }, { &pTransmittanceRTV.p }));

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeTransmittanceTex2DTech"];

	VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

	ID3D11RenderTargetView* pRTVs[] =
	{
		pTransmittanceRTV.p,
	};
	pContext->OMSetRenderTargets(ARRAYSIZE(pRTVs), pRTVs, nullptr);

	RenderQuad(pContext, activeTech, TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
	//pContext->OMSetRenderTargets(0, nullptr, nullptr);
	UnbindResources(pContext);
#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pTransmittanceTex2D, D3DX11_IFF_DDS, L"Texture/Transmittance.dds"));
#endif
	return hr;
}


HRESULT Atmosphere::PreComputeOpticalLengthTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	pOpticalLengthTex2D.Release();
	pOpticalLengthSRV.Release();
	CComPtr<ID3D11RenderTargetView>	pOpticalLengthRTV;
	V_RETURN(CreateTexture2D(pDevice, pContext, TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT, format,
	{ &pOpticalLengthTex2D.p }, { &pOpticalLengthSRV.p }, { &pOpticalLengthRTV.p }));

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeOpticalLengthTex2DTech"];

	VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

	ID3D11RenderTargetView* pRTVs[] =
	{
		pOpticalLengthRTV.p,
	};
	pContext->OMSetRenderTargets(ARRAYSIZE(pRTVs), pRTVs, nullptr);

	RenderQuad(pContext, activeTech, TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
	//pContext->OMSetRenderTargets(0, nullptr, nullptr);
	UnbindResources(pContext);
#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pTransmittanceTex2D, D3DX11_IFF_DDS, L"Texture/OpticalLength.dds"));
#endif
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
	{ &pDirectIrradianceTex2D.p }, { &pDirectIrradianceSRV.p }, { &pDirectIrradianceRTV.p }));

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeDirectIrradiance2DTech"];

	ID3D11RenderTargetView* pRTVs[] =
	{
		pDirectIrradianceRTV.p,
	};
	pContext->OMSetRenderTargets(ARRAYSIZE(pRTVs), pRTVs, nullptr);

	ShaderResourceVarMap["g_tex2DTransmittanceLUT"]->SetResource(pTransmittanceSRV);

	RenderQuad(pContext, activeTech, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
	UnbindResources(pContext);
#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pDirectIrradianceTex2D, D3DX11_IFF_DDS, L"Texture/DirectIrradiance.dds"));
#endif
	return hr;
}


HRESULT Atmosphere::PreComputeSingleSctrTex3D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	pSingleScatterTex3D.Release();
	pSingleScatterCombinedTex3D.Release();
	pSingleScatterMieTex3D.Release();
	pSingleScatterSRV.Release();
	pSingleScatterCombinedSRV.Release();
	pSingleScatterMieSRV.Release();
	V_RETURN(CreateTexture3D(pDevice, pContext, 
							SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH,format,
							{ &pSingleScatterTex3D.p,&pSingleScatterCombinedTex3D.p,&pSingleScatterMieTex3D.p },
							{ &pSingleScatterSRV.p,&pSingleScatterCombinedSRV.p,&pSingleScatterMieSRV.p }));
	
	ShaderResourceVarMap["g_tex2DTransmittanceLUT"]->SetResource(pTransmittanceSRV);
	ShaderResourceVarMap["g_tex2DOpticalLengthLUT"]->SetResource(pOpticalLengthSRV);

	std::vector<CComPtr<ID3D11RenderTargetView>> pSingleScaterRTVs(SCATTERING_TEXTURE_DEPTH);
	std::vector<CComPtr<ID3D11RenderTargetView>> pSingleScaterCombinedRTVs(SCATTERING_TEXTURE_DEPTH);
	std::vector<CComPtr<ID3D11RenderTargetView>> pSingleScaterMieRTVs(SCATTERING_TEXTURE_DEPTH);	
	MiscDynamicParams misc;
	for(UINT depthSlice = 0;depthSlice<SCATTERING_TEXTURE_DEPTH;++depthSlice)
	{
		D3D11_RENDER_TARGET_VIEW_DESC CurrSliceRTVDesc;
		CurrSliceRTVDesc.Format = format;
		CurrSliceRTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
		CurrSliceRTVDesc.Texture3D.MipSlice = 0;
		CurrSliceRTVDesc.Texture3D.FirstWSlice = depthSlice;
		CurrSliceRTVDesc.Texture3D.WSize = 1;

		V_RETURN(pDevice->CreateRenderTargetView(pSingleScatterTex3D, &CurrSliceRTVDesc, &pSingleScaterRTVs[depthSlice]));
		V_RETURN(pDevice->CreateRenderTargetView(pSingleScatterCombinedTex3D, &CurrSliceRTVDesc, &pSingleScaterCombinedRTVs[depthSlice]));		

		V_RETURN(pDevice->CreateRenderTargetView(pSingleScatterMieTex3D, &CurrSliceRTVDesc, &pSingleScaterMieRTVs[depthSlice]));

		ID3DX11EffectTechnique* activeTech = TechMap["ComputeSingleScatterTex3DTech"];

		VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

	
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
			pSingleScaterCombinedRTVs[depthSlice].p,
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
	UnbindResources(pContext);
#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pSingleScatterTex3D, D3DX11_IFF_DDS, L"Texture/SingleScatter/SingleScatter3D.dds"));
	V_RETURN(D3DX11SaveTextureToFile(pContext, pSingleScatterCombinedTex3D, D3DX11_IFF_DDS, L"Texture/SingleScatter/SingleScatterCombined3D.dds"));
#endif
	return hr;
}


HRESULT Atmosphere::PreComputeInDirectIrradianceTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, UINT scatter_order)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	CComPtr<ID3D11RenderTargetView>	pIndirectIrradianceRTV;
	pIndirectIrradianceTex2D.Release();
	pIndirectIrradianceSRV.Release();
	V_RETURN(CreateTexture2D(pDevice, pContext, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT,
		format, { &pIndirectIrradianceTex2D.p }, { &pIndirectIrradianceSRV.p }, { &pIndirectIrradianceRTV.p }));

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeIndirectIrradiance2DTech"];

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
	UnbindResources(pContext);
#if CREATE_TEXTURE_DDS_TEST
	std::wstringstream ss;
	ss << "Texture/Irradiance/IndirectIrradiance_" << scatter_order << ".dds";
	V_RETURN(D3DX11SaveTextureToFile(pContext, pIndirectIrradianceTex2D, D3DX11_IFF_DDS, ss.str().c_str()));
#endif
	return hr;
}


HRESULT Atmosphere::PreComputeMultiSctrTex3D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, UINT scatter_order)
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
	ShaderResourceVarMap["g_tex2DOpticalLengthLUT"]->SetResource(pOpticalLengthSRV);
	ShaderResourceVarMap["g_tex3DSingleMieScatteringLUT"]->SetResource(pSingleScatterMieSRV);
	ShaderResourceVarMap["g_tex3DSingleScatteringCombinedLUT"]->SetResource(pSingleScatterCombinedSRV);

	CComPtr<ID3D11Texture3D> pTex3D[2];
	CComPtr<ID3D11ShaderResourceView> pSRV[2];
	V_RETURN(CreateTexture3D(pDevice, pContext, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH,
							format,{ &pTex3D[0].p,&pTex3D[1].p },{ &pSRV[0].p,&pSRV[1].p }));

	std::vector<CComPtr<ID3D11RenderTargetView>> pMultiScaterRTVs(SCATTERING_TEXTURE_DEPTH);
	std::vector<CComPtr<ID3D11RenderTargetView>> pMultiScaterCombinedRTVs(SCATTERING_TEXTURE_DEPTH);
	MiscDynamicParams misc;
	for (UINT depthSlice = 0; depthSlice<SCATTERING_TEXTURE_DEPTH; ++depthSlice)
	{
		D3D11_RENDER_TARGET_VIEW_DESC CurrSliceRTVDesc;
		CurrSliceRTVDesc.Format = format;
		CurrSliceRTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
		CurrSliceRTVDesc.Texture3D.MipSlice = 0;
		CurrSliceRTVDesc.Texture3D.FirstWSlice = depthSlice;
		CurrSliceRTVDesc.Texture3D.WSize = 1;

		V_RETURN(pDevice->CreateRenderTargetView(pTex3D[0], &CurrSliceRTVDesc, &pMultiScaterRTVs[depthSlice]));
		V_RETURN(pDevice->CreateRenderTargetView(pTex3D[1], &CurrSliceRTVDesc, &pMultiScaterCombinedRTVs[depthSlice]));
		
		ID3DX11EffectTechnique* activeTech = TechMap["ComputeMultiScatterTex3DTech"];

		VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

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
			pMultiScaterCombinedRTVs[depthSlice].p
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
	UnbindResources(pContext);

#if CREATE_TEXTURE_DDS_TEST
	std::wstringstream ss;
	ss << "Texture/MultiScatter/MultiScatter_" << scatter_order << ".dds";
	V_RETURN(D3DX11SaveTextureToFile(pContext, pTex3D[0], D3DX11_IFF_DDS, ss.str().c_str()));
	ss.str(L"");
	ss.clear();
	ss << "Texture/MultiScatter/MultiScatterCombined_" << scatter_order << ".dds";
	V_RETURN(D3DX11SaveTextureToFile(pContext, pTex3D[1], D3DX11_IFF_DDS, ss.str().c_str()));
#endif
	pMultiScatterTex3D.Release();
	pMultiScatterSRV.Release();

	pMultiScatterCombinedTex3D.Release();
	pMultiScatterCombinedSRV.Release();

	pMultiScatterTex3D = pTex3D[0];
	pMultiScatterSRV = pSRV[0];

	pMultiScatterCombinedTex3D = pTex3D[1];
	pMultiScatterCombinedSRV = pSRV[1];

	return hr;
}


void Atmosphere::RenderSun(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	if (fEnableLightShaft)
	{
		ID3DX11EffectTechnique* activeTech = TechMap["RenderSunTech"];
		RenderQuad(pContext, activeTech, screen_width, screen_height);
	}
}


HRESULT Atmosphere::ComputeSpaceLinearDepthTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11ShaderResourceView* depthSRV)
{
	HRESULT hr = S_OK;
	ID3DX11EffectTechnique* activeTech = TechMap["ComputeSpaceLinearDepthTex2DTech"];

	CComPtr<ID3D11RenderTargetView>	pRTV;

	DXGI_FORMAT format = DXGI_FORMAT_R32_FLOAT;

	pSpaceLinearDepthTex2D.Release();
	pSpaceLinearDepthSRV.Release();
	V_RETURN(CreateTexture2D(pDevice, pContext, screen_width,screen_height,format,
				{ &pSpaceLinearDepthTex2D.p }, { &pSpaceLinearDepthSRV.p }, { &pRTV.p }));

	ShaderResourceVarMap["g_tex2DSpaceDepth"]->SetResource(depthSRV);

	pContext->OMSetRenderTargets(1, &pRTV.p, nullptr);
	RenderQuad(pContext, activeTech, screen_width, screen_height);

	UnbindResources(pContext);
	pContext->OMSetRenderTargets(0, nullptr, nullptr);
#if CREATE_TEXTURE_DDS_TEST
	//V_RETURN(D3DX11SaveTextureToFile(pContext, pSpaceLinearDepthTex2D, D3DX11_IFF_DDS, L"Texture/SpaceLinearDepth.dds"));
#endif

	return hr;
}


HRESULT Atmosphere::ComputeSliceEndTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeSliceEndTex2DTech"];

	CComPtr<ID3D11RenderTargetView>	pRTV;
	DXGI_FORMAT format = DXGI_FORMAT_R32G32B32A32_FLOAT;
	pSliceEndTex2D.Release();
	pSliceEndSRV.Release();
	V_RETURN(CreateTexture2D(pDevice, pContext, EPIPOLAR_SLICE_NUM, 1, format,
		{ &pSliceEndTex2D.p }, { &pSliceEndSRV.p }, { &pRTV.p }));

	pContext->OMSetRenderTargets(1, &pRTV.p, nullptr);
	RenderQuad(pContext, activeTech, EPIPOLAR_SLICE_NUM, 1);

	pContext->OMSetRenderTargets(0, nullptr, nullptr);
#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pSliceEndTex2D, D3DX11_IFF_DDS, L"Texture/SliceEndTex.dds"));
#endif

	return hr;
}


HRESULT Atmosphere::ComputeEpipolarCoordTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeEpipolarCoordTex2DTech"];

	CComPtr<ID3D11RenderTargetView>	pEpipolarSampleRTV,pEpipolarSampleCamDepthRTV;

	pEpipolarSampleTex2D.Release();
	pEpipolarSampleSRV.Release();
	pEpipolarSampleCamDepthTex2D.Release();
	pEpipolarSampleCamDepthSRV.Release();
	V_RETURN(CreateTexture2D(pDevice, pContext, EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM, DXGI_FORMAT_R32G32_FLOAT,
					{ &pEpipolarSampleTex2D.p }, 
					{ &pEpipolarSampleSRV.p }, 
					{ &pEpipolarSampleRTV.p}));
	V_RETURN(CreateTexture2D(pDevice, pContext, EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM, DXGI_FORMAT_R32_FLOAT,
					{ &pEpipolarSampleCamDepthTex2D.p },
					{ &pEpipolarSampleCamDepthSRV.p },
					{ &pEpipolarSampleCamDepthRTV.p }));

	D3D11_TEXTURE2D_DESC TexDes;
	pEpipolarSampleTex2D->GetDesc(&TexDes);
	TexDes.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
	TexDes.BindFlags = D3D11_BIND_DEPTH_STENCIL;
	CComPtr<ID3D11Texture2D> pStencilTex;
	pEpipolarSampleDSV.Release();
	V_RETURN(pDevice->CreateTexture2D(&TexDes, nullptr, &pStencilTex));
	V_RETURN(pDevice->CreateDepthStencilView(pStencilTex, nullptr, &pEpipolarSampleDSV));

	ShaderResourceVarMap["g_tex2DSpaceLinearDepth"]->SetResource(pSpaceLinearDepthSRV);
	ShaderResourceVarMap["g_tex2DSliceEnd"]->SetResource(pSliceEndSRV);
	
	ID3D11RenderTargetView* pRTVs[] =
	{
		pEpipolarSampleRTV.p,
		pEpipolarSampleCamDepthRTV.p
	};
	int num = ARRAYSIZE(pRTVs);
	static const float fInvalidCoordinate = -1e+30f; // Both coord texture and epipolar CamSpaceZ are 32-bit float
	float InvalidCoords[] = { fInvalidCoordinate, fInvalidCoordinate, fInvalidCoordinate, fInvalidCoordinate };
	for(int i = 0; i<num;++i)
	{
		pContext->ClearRenderTargetView(pRTVs[i], InvalidCoords);
	}
	pContext->OMSetRenderTargets(num, pRTVs, pEpipolarSampleDSV);
	pContext->ClearDepthStencilView(pEpipolarSampleDSV, D3D11_CLEAR_STENCIL, 1.0f, 0);
	RenderQuad(pContext, activeTech,EPIPOLAR_SAMPLE_NUM , EPIPOLAR_SLICE_NUM);
	UnbindResources(pContext);
	pContext->OMSetRenderTargets(0, nullptr, nullptr);
#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pEpipolarSampleTex2D, D3DX11_IFF_DDS, L"Texture/EpipolarSampleTex.dds"));
	V_RETURN(D3DX11SaveTextureToFile(pContext, pEpipolarSampleCamDepthTex2D, D3DX11_IFF_DDS, L"Texture/EpipolarSampleDepth.dds"));
#endif
	return hr;
}


HRESULT Atmosphere::ComputeUnshadowedSampleScatter(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;
	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;
	pUnshadowedSampleScatterTex2D.Release();
	pUnshadowedSampleScatterSRV.Release();
	CComPtr<ID3D11RenderTargetView> pUnshadowedSampleScatterRTV;
	V_RETURN(CreateTexture2D(pDevice, pContext, EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM, format,
			{ &pUnshadowedSampleScatterTex2D.p },
			{ &pUnshadowedSampleScatterSRV.p },
			{ &pUnshadowedSampleScatterRTV.p }));
	misc.scatter_order = 2;
	VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));
	ShaderResourceVarMap["g_tex2DEpipolarSample"]->SetResource(pEpipolarSampleSRV);
	ShaderResourceVarMap["g_tex2DEpipolarSampleCamDepth"]->SetResource(pEpipolarSampleCamDepthSRV);
	ShaderResourceVarMap["g_tex2DOpticalLengthLUT"]->SetResource(pOpticalLengthSRV);
	ShaderResourceVarMap["g_tex3DMultiScatteringCombinedLUT"]->SetResource(pMultiScatterCombinedSRV);

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeUnshadowedSampleScatterTech"];
	pContext->OMSetRenderTargets(1, &pUnshadowedSampleScatterRTV.p, nullptr);
	RenderQuad(pContext, activeTech, EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM);
	UnbindResources(pContext);
	return hr;
}


HRESULT Atmosphere::RefineSampleLocal(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	pInterpolationSampleTex2D.Release();
	pInterpolationSampleUAV.Release();
	pInterpolationSampleSRV.Release();

	D3D11_TEXTURE2D_DESC PreCompute2DTexDesc;
	ZeroMemory(&PreCompute2DTexDesc, sizeof(PreCompute2DTexDesc));
	PreCompute2DTexDesc.Width = EPIPOLAR_SAMPLE_NUM;
	PreCompute2DTexDesc.Height = EPIPOLAR_SLICE_NUM;
	PreCompute2DTexDesc.MipLevels = 1;
	PreCompute2DTexDesc.ArraySize = 1;
	PreCompute2DTexDesc.Format = DXGI_FORMAT_R16G16_UINT;
	PreCompute2DTexDesc.SampleDesc.Count = 1;
	PreCompute2DTexDesc.SampleDesc.Quality = 0;
	PreCompute2DTexDesc.Usage = D3D11_USAGE_DEFAULT;
	PreCompute2DTexDesc.BindFlags = D3D11_BIND_UNORDERED_ACCESS | D3D11_BIND_SHADER_RESOURCE;
	PreCompute2DTexDesc.CPUAccessFlags = 0;
	PreCompute2DTexDesc.MiscFlags = 0;

	V_RETURN(pDevice->CreateTexture2D(&PreCompute2DTexDesc, nullptr, &pInterpolationSampleTex2D));
	V_RETURN(pDevice->CreateShaderResourceView(pInterpolationSampleTex2D, nullptr, &pInterpolationSampleSRV));
	V_RETURN(pDevice->CreateUnorderedAccessView(pInterpolationSampleTex2D, nullptr, &pInterpolationSampleUAV));
	
	RefineSampleCSThreadGroupSize = max(RefineSampleCSThreadGroupSize, InterpolationSampleStep);
	RefineSampleCSThreadGroupSize = min(RefineSampleCSThreadGroupSize, EPIPOLAR_SAMPLE_NUM);

	ID3DX11EffectTechnique* activeTech = TechMap["RefineSampleTech"];

	//CComPtr<ID3DBlob> pShaderByteCode;
	//CComPtr<ID3D11ComputeShader> pComputeShader;
	//std::string sThreadGroupSize = std::to_string(RefineSampleCSThreadGroupSize);
	//std::string sSampleStep = std::to_string(InterpolationSampleStep);
	//std::vector<D3D_SHADER_MACRO> Macros
	//{
	//	{ "THREAD_GROUP_SIZE", sThreadGroupSize.c_str() },
	//	{ "SAMPLE_STEP", sSampleStep.c_str() },
	//	{ nullptr,nullptr }
	//};
	//std::vector<D3D_SHADER_MACRO> Macros;
	//AddMarco(Macros, "THREAD_GROUP_SIZE", RefineSampleCSThreadGroupSize);
	//AddMarco(Macros, "SAMPLE_STEP", InterpolationSampleStep);
	//FinishMarco(Macros);

	//V_RETURN(CompileShaderFromFile(L"Atmosphere.fx", "RefineSampleCS", &Macros[0], "cs_5_0", &pShaderByteCode));
	//V_RETURN(pDevice->CreateComputeShader(pShaderByteCode->GetBufferPointer(), pShaderByteCode->GetBufferSize(), nullptr, &pComputeShader));

	ShaderResourceVarMap["g_tex2DEpipolarSample"]->SetResource(pEpipolarSampleSRV);
	ShaderResourceVarMap["g_tex2DEpipolarSampleCamDepth"]->SetResource(pEpipolarSampleCamDepthSRV);
	ShaderResourceVarMap["g_tex2DUnshadowedSampleScatter"]->SetResource(pUnshadowedSampleScatterSRV);
	//ShaderResourceVarMap["g_rwtex2DInterpolationSource"]->SetResource(pInterpolationSampleUAV);
	CComPtr<ID3DX11EffectUnorderedAccessViewVariable> pUAV = pEffect->GetVariableByName("g_rwtex2DInterpolationSource")->AsUnorderedAccessView();
	pUAV->SetUnorderedAccessView(pInterpolationSampleUAV);

	D3DX11_TECHNIQUE_DESC techDesc;
	activeTech->GetDesc(&techDesc);
	ID3D11RenderTargetView *pDummyRTV = nullptr;
	pContext->OMSetRenderTargets(1, &pDummyRTV, nullptr);
	for (UINT p = 0; p<techDesc.Passes; ++p)
	{
		CComPtr<ID3DX11EffectPass> pass = activeTech->GetPassByIndex(p);
		//pContext->CSSetShader(pComputeShader, nullptr, 0);
		pass->Apply(0, pContext);
		//pContext->CSSetUnorderedAccessViews(0, 1, &pInterpolationSampleUAV.p, nullptr);
		//pContext->CSSetShaderResources(0, 1, &pEpipolarSampleSRV.p);
		//pContext->CSSetShaderResources(1, 1, &pEpipolarSampleCamDepthSRV.p);
		pContext->Dispatch(EPIPOLAR_SAMPLE_NUM / RefineSampleCSThreadGroupSize,
			EPIPOLAR_SLICE_NUM, 1);
	}

	
#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pInterpolationSampleTex2D, D3DX11_IFF_DDS, L"Texture/InterpolationSampleTex.dds"));
#endif
	UnbindResources(pContext);
	return hr;
}


HRESULT Atmosphere::ComputeSliceUVOrigDirTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11ShaderResourceView* pShadowMapSRV)
{
	HRESULT hr = S_OK;
	DXGI_FORMAT format = DXGI_FORMAT_R32G32B32A32_FLOAT;
	pSliceUVOrigDirTex2D.Release();
	pSliceUVOrigDirSRV.Release();
	CComPtr<ID3D11RenderTargetView> pSliceUVOrigDirRTV;
	V_RETURN(CreateTexture2D(pDevice, pContext, EPIPOLAR_SLICE_NUM, 1, format,
							{ &pSliceUVOrigDirTex2D.p }, { &pSliceUVOrigDirSRV.p }, 
							{ &pSliceUVOrigDirRTV.p }));

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeSliceUVOrigDirTex2DTech"];
	ShaderResourceVarMap["g_tex2DSliceEnd"]->SetResource(pSliceEndSRV);

	pContext->OMSetRenderTargets(1, &pSliceUVOrigDirRTV.p, nullptr);
	RenderQuad(pContext, activeTech, EPIPOLAR_SLICE_NUM, 1);
	
	UnbindResources(pContext);
	return hr;
}


HRESULT Atmosphere::Build1DMinMaxMipMap(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11ShaderResourceView* pShadowMapSRV,UINT shadowMapResolution)
{
	HRESULT hr = S_OK;
	DXGI_FORMAT format = DXGI_FORMAT_R32G32_FLOAT;
	UINT MIN_MAX_TEXTURE_DIM = shadowMapResolution;
	for (int i = 0; i < 2; i++)
	{
		pMinMaxMinMapTex2D[i].Release();
		pMinMaxMinMapTexSRV[i].Release();
	}
	CComPtr<ID3D11RenderTargetView> pMinMaxMinMapTexRTV[2];
	VarMap["MIN_MAX_TEXTURE_DIM"]->SetRawValue(&MIN_MAX_TEXTURE_DIM, 0, sizeof(UINT));
	V_RETURN(CreateTexture2D(pDevice, pContext, MIN_MAX_TEXTURE_DIM, EPIPOLAR_SLICE_NUM, format,
							{ &pMinMaxMinMapTex2D[0].p,&pMinMaxMinMapTex2D[1].p }, 
							{ &pMinMaxMinMapTexSRV[0].p,&pMinMaxMinMapTexSRV[1].p },
							{ &pMinMaxMinMapTexRTV[0].p,&pMinMaxMinMapTexRTV[1].p }));


	ID3DX11EffectTechnique* activeTech = TechMap["Initial1DMinMaxMipMapTech"];
	ShaderResourceVarMap["g_tex2DShadowMap"]->SetResource(pShadowMapSRV);
	ShaderResourceVarMap["g_tex2DSliceUVOrigDir"]->SetResource(pSliceUVOrigDirSRV);
	pContext->OMSetRenderTargets(1, &pMinMaxMinMapTexRTV[0].p, nullptr);
	float clearColor[] = { 1.0,0.0,0.0,1.0 };
	pContext->ClearRenderTargetView(pMinMaxMinMapTexRTV[0], clearColor);
	RenderQuad(pContext, activeTech, shadowMapResolution / 2, EPIPOLAR_SLICE_NUM);
	UnbindResources(pContext);
	CComPtr<ID3D11RenderTargetView> pDummyRTV = nullptr;
	pContext->OMSetRenderTargets(1, &pDummyRTV, nullptr);

	activeTech = TechMap["Compute1DMinMaxMipMapLevelTech"];
	CComPtr<ID3D11Resource> pDstResource, pSrcResource;
	pMinMaxMinMapTexRTV[0]->GetResource(&pDstResource.p);
	pMinMaxMinMapTexRTV[1]->GetResource(&pSrcResource.p);

	UINT offsetX = shadowMapResolution / 2;
	UINT preOffsetX = 0;
	for (UINT level = 2, step = 4; level <= MaxMinMaxMapLevel; level++,step <<= 1)
	{
		ShaderResourceVarMap["g_tex2DMinMaxMipMap"]->SetResource(pMinMaxMinMapTexSRV[level % 2]);
		misc.uiSrcMinMaxOffsetX = preOffsetX;
		misc.uiSrcMinMaxOffsetY = 0;
		misc.uiDstMinMaxOffsetX = offsetX;
		misc.uiDstMinMaxOffsetY = 0;
		VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));
		pContext->OMSetRenderTargets(1, &pMinMaxMinMapTexRTV[(level + 1) % 2].p, nullptr);
		RenderQuad(pContext, activeTech, shadowMapResolution / step, EPIPOLAR_SLICE_NUM, offsetX, 0);

		preOffsetX = offsetX;
		offsetX += shadowMapResolution / step;
		if((level + 1) % 2)
		{
			D3D11_BOX SrcBox;
			SrcBox.left = preOffsetX;
			SrcBox.right = offsetX;
			SrcBox.top = 0;
			SrcBox.bottom = EPIPOLAR_SLICE_NUM;
			SrcBox.front = 0;
			SrcBox.back = 1;
			pContext->CopySubresourceRegion(pDstResource, 0, preOffsetX, 0, 0, pSrcResource, 0, &SrcBox);
		}
		pContext->OMSetRenderTargets(1, &pDummyRTV, nullptr);
		UnbindResources(pContext);
	}

	UnbindResources(pContext);
	return hr;
}


HRESULT Atmosphere::MarkRayMarchSample(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	ID3DX11EffectTechnique* activeTech = TechMap["MarkRayMarchSampleTech"];
	ShaderResourceVarMap["g_tex2DInterpolationSample"]->SetResource(pInterpolationSampleSRV);

	ID3D11RenderTargetView *pDummyRTV = nullptr;
	pContext->OMSetRenderTargets(1, &pDummyRTV, pEpipolarSampleDSV);
	RenderQuad(pContext, activeTech, EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM);
	pContext->OMSetRenderTargets(1, &pDummyRTV, nullptr);
	UnbindResources(pContext);
	return hr;
}


HRESULT Atmosphere::DoRayMarch(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11ShaderResourceView* pShadowMapSRV)
{
	HRESULT hr = S_OK;
	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;
	pSampleScatterTex2D.Release();
	pSampleScatterSRV.Release();
	CComPtr<ID3D11RenderTargetView> pSampleScatterRTV;
	V_RETURN(CreateTexture2D(pDevice, pContext, EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM, format,
			{ &pSampleScatterTex2D.p }, { &pSampleScatterSRV.p }, { &pSampleScatterRTV.p }));

	ID3DX11EffectTechnique* activeTech = TechMap["DoRayMarchTech"];
	ShaderResourceVarMap["g_tex2DEpipolarSample"]->SetResource(pEpipolarSampleSRV);
	ShaderResourceVarMap["g_tex2DEpipolarSampleCamDepth"]->SetResource(pEpipolarSampleCamDepthSRV);
	ShaderResourceVarMap["g_tex2DOpticalLengthLUT"]->SetResource(pOpticalLengthSRV);
	ShaderResourceVarMap["g_tex3DMultiScatteringCombinedLUT"]->SetResource(pMultiScatterCombinedSRV);

	misc.fEnableLightShaft = fEnableLightShaft;
	misc.fIsLightInSpaceCorrect = fIsLightInSpaceCorrect;
	misc.uiMinMaxLevelMax = MaxMinMaxMapLevel;
	misc.scatter_order = 2;
	VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));

	pContext->OMSetRenderTargets(1, &pSampleScatterRTV.p, pEpipolarSampleDSV);
	RenderQuad(pContext, activeTech, EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM);
	ID3D11RenderTargetView *pDummyRTV = nullptr;
	pContext->OMSetRenderTargets(1, &pDummyRTV, nullptr);
	UnbindResources(pContext);

	return hr;
}


HRESULT Atmosphere::InterpolateScatter(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;
	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;
	pInterpolatedSampleScatterTex2D.Release();
	pInterpolatedSampleScatterSRV.Release();
	CComPtr<ID3D11RenderTargetView> pInterpolateSampleScatterRTV;
	V_RETURN(CreateTexture2D(pDevice, pContext, EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM, format,
							{ &pInterpolatedSampleScatterTex2D.p }, { &pInterpolatedSampleScatterSRV.p }, 
							{ &pInterpolateSampleScatterRTV.p }));

	ID3DX11EffectTechnique* activeTech = TechMap["InterpolateScatterTech"];
	ShaderResourceVarMap["g_tex2DSampleScatter"]->SetResource(pSampleScatterSRV);
	ShaderResourceVarMap["g_tex2DInterpolationSample"]->SetResource(pInterpolationSampleSRV);

	misc.fEnableLightShaft = fEnableLightShaft;
	misc.fIsLightInSpaceCorrect = fIsLightInSpaceCorrect;
	VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));

	pContext->OMSetRenderTargets(1, &pInterpolateSampleScatterRTV.p, nullptr);
	RenderQuad(pContext, activeTech, EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM);
	ID3D11RenderTargetView *pDummyRTV = nullptr;
	pContext->OMSetRenderTargets(1, &pDummyRTV, nullptr);
	UnbindResources(pContext);
	return hr;
}


HRESULT Atmosphere::ApplyAndFixInterpolateScatter(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, 
												ID3D11RenderTargetView* pRTV, ID3D11ShaderResourceView* pColorBufferSRV)
{
	HRESULT hr = S_OK;

	D3D11_TEXTURE2D_DESC TexDes;
	ZeroMemory(&TexDes, sizeof(TexDes));
	TexDes.Width = screen_width;
	TexDes.Height = screen_height;
	TexDes.MipLevels = 1;
	TexDes.ArraySize = 1;
	TexDes.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
	TexDes.SampleDesc.Count = 1;
	TexDes.SampleDesc.Quality = 0;
	TexDes.BindFlags = D3D11_BIND_DEPTH_STENCIL;
	TexDes.Usage = D3D11_USAGE_DEFAULT;
	TexDes.CPUAccessFlags = 0;
	TexDes.MiscFlags = 0;
	CComPtr<ID3D11Texture2D> pStencilTex;
	pApplyScatterDSV.Release();
	V_RETURN(pDevice->CreateTexture2D(&TexDes, nullptr, &pStencilTex));
	V_RETURN(pDevice->CreateDepthStencilView(pStencilTex, nullptr, &pApplyScatterDSV));

	ID3DX11EffectTechnique* activeTech = TechMap["ApplyInterpolateScatterTech"];
	ShaderResourceVarMap["g_tex2DColorBuffer"]->SetResource(pColorBufferSRV);
	ShaderResourceVarMap["g_tex2DSliceEnd"]->SetResource(pSliceEndSRV);
	ShaderResourceVarMap["g_tex2DSpaceLinearDepth"]->SetResource(pSpaceLinearDepthSRV);
	ShaderResourceVarMap["g_tex2DEpipolarSampleCamDepth"]->SetResource(pEpipolarSampleCamDepthSRV);
	ShaderResourceVarMap["g_tex2DInterpolatedScatter"]->SetResource(pInterpolatedSampleScatterSRV);
	ShaderResourceVarMap["g_tex2DEarthGround"]->SetResource(pEarthGroundSRV);
	ShaderResourceVarMap["g_tex2DDirectIrradianceLUT"]->SetResource(pDirectIrradianceSRV);
	ShaderResourceVarMap["g_tex2DIndirectIrradianceLUT"]->SetResource(pIndirectIrradianceSRV);

	MiscDynamicParams misc;
	misc.fEnableLightShaft = fEnableLightShaft;
	misc.fIsLightInSpaceCorrect = fIsLightInSpaceCorrect;
	VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));

	pContext->ClearDepthStencilView(pApplyScatterDSV, D3D11_CLEAR_STENCIL, 1.0f, 0);
	pContext->OMSetRenderTargets(1, &pRTV, pApplyScatterDSV);
	RenderQuad(pContext, activeTech, screen_width, screen_height);

	activeTech = TechMap["FixInterpolateScatterTech"];
	RenderQuad(pContext, activeTech, screen_width, screen_height);

	UnbindResources(pContext);
	return hr;
}


void Atmosphere::SetView(float view_distance_meters, float view_zenith_angle_radians, float view_azimuth_angle_radians,
							float sun_zenith_angle_radians, float sun_azimuth_angle_radians,float exposure)
{
	this->view_distance_meters = view_distance_meters;
	this->view_zenith_angle_radians = view_zenith_angle_radians;
	this->view_azimuth_angle_radians = view_azimuth_angle_radians;
	this->sun_zenith_angle_radians = sun_zenith_angle_radians;
	this->sun_azimuth_angle_radians = sun_azimuth_angle_radians;
	this->exposure = exposure;

	float cos_z = cos(view_zenith_angle_radians);
	float sin_z = sin(view_zenith_angle_radians);
	float cos_a = cos(view_azimuth_angle_radians);
	float sin_a = sin(view_azimuth_angle_radians);

	D3DXVECTOR3 EyePos = D3DXVECTOR3(sin_z*cos_a, cos_z, sin_z*sin_a) * view_distance_meters / 1000;
	D3DXVECTOR3 LookAt = D3DXVECTOR3(0.0f, 0.0f, 0.0f);
	//m_FirstPersonCamera.SetViewParams(&EyePos, &LookAt);

}


void Atmosphere::MsgProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	//m_FirstPersonCamera.HandleMessages(hWnd, uMsg, wParam, lParam);

	switch (uMsg)
	{
		case WM_KEYDOWN:
		{
			switch(wParam)
			{
				case 0x31:
				{
					SetView(9000.0, 1.47, -0.1, 1.3, 2.9, 10.0);
				}break;
				case 0x32:
				{
					SetView(9000.0, 1.47, 0.0, 1.564, -3.0, 10.0);
				}break;
				case 0x33:
				{
					SetView(7000.0, 1.57, 0.0, 1.54, -2.96, 10.0);
				}break;
				case 0x34:
				{
					SetView(7000.0, 1.57, 0.0, 1.328, -3.044, 10.0);
				}break;
				case 0x35:
				{
					SetView(9000.0, 1.39, 0.0, 1.2, 0.7, 10.0);
				}break;
				case 0x36:
				{
					SetView(9000.0, 1.5, 0.0, 1.628, 1.05, 200.0);
				}break;
				case 0x37:
				{
					SetView(7000.0, 1.43, 0.0, 1.57, 1.34, 40.0);
				}break;
				case 0x38:
				{
					SetView(2.7e6, 0.81, 0.0, 1.57, 2.0, 10.0);
				}break;
				case 0x39:
				{
					SetView(1.2e7, 0.0, 0.0, 0.93, -2.0, 10.0);
				}break;
			}
		}
	}
}


void Atmosphere::OnFrameMove(double fTime, float fElapsedTime, float fScale)
{
	//m_FirstPersonCamera.FrameMove(fElapsedTime);
	if(GetAsyncKeyState('I') & 0x8000)
	{
		sun_zenith_angle_radians -= fElapsedTime / fScale;
	}
	if (GetAsyncKeyState('K') & 0x8000)
	{
		sun_zenith_angle_radians += fElapsedTime / fScale;
	}
	if (GetAsyncKeyState('J') & 0x8000)
	{
		sun_azimuth_angle_radians += fElapsedTime / fScale;
	}
	if (GetAsyncKeyState('L') & 0x8000)
	{
		sun_azimuth_angle_radians -= fElapsedTime / fScale;
	}
}


void Atmosphere::Resize(int screen_width, int screen_height, float fFOV, float fAspect, float fNear, float fFar)
{
	this->screen_width = screen_width;
	this->screen_height = screen_height;

	VarMap["SCREEN_WIDTH"]->SetRawValue(&screen_width, 0, sizeof(int));
	VarMap["SCREEN_HEIGHT"]->SetRawValue(&screen_height, 0, sizeof(int));

	//m_FirstPersonCamera.SetProjParams(fFOV, fAspect, fNear, fFar);
}


D3DXVECTOR3 Atmosphere::GetSunDir()
{
	float cos_sun_z = cos(sun_zenith_angle_radians);
	float sin_sun_z = sin(sun_zenith_angle_radians);
	float cos_sun_a = cos(sun_azimuth_angle_radians);
	float sin_sun_a = sin(sun_azimuth_angle_radians);
	return D3DXVECTOR3(sin_sun_z*cos_sun_a, cos_sun_z, sin_sun_z*sin_sun_a);
}

void Atmosphere::SetLightParam(const D3DXMATRIX& lightView, 
								const D3DXMATRIX& lightProj)
{
	this->lightView = lightView;
	this->lightProj = lightProj;
	this->lightProj.m[0][0] *= 1000;
	this->lightProj.m[1][1] *= 1000;
	this->lightProj.m[2][2] *= 1000;
}

void Atmosphere::SetCamParam(const D3DXVECTOR3& f3CamPos, const D3DXVECTOR3& f3CamDir, 
							const D3DXMATRIX& camView, const D3DXMATRIX& camProj,
							float fCamNear,float fCamFar)
{
	this->f3CamPos = f3CamPos / 1000;
	this->f3CamDir = f3CamDir;
	this->camView = camView;
	D3DXVECTOR3* pos = (D3DXVECTOR3*)(&this->camView.m[3][0]);
	*pos = *pos / 1000;

	this->camProj = camProj;
	this->camProj.m[3][2] = this->camProj.m[3][2] / 1000;

	this->fCamNear = fCamNear / 1000;
	this->fCamFar = fCamFar / 1000;

	//this->f3CamPos = f3CamPos;
	//this->f3CamDir = f3CamDir;
	//this->camView = camView;

	//this->camProj = camProj;
	//this->camProj.m[3][2] = this->camProj.m[3][2];

	//this->fCamNear = fCamNear;
	//this->fCamFar = fCamFar;
}