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
	atmosphereParams.sun_angular_radius = 0.2678 * D3DX_PI / 180.0;
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
	GameObject::Release();
}


HRESULT Atmosphere::OnD3D11CreateDevice(ID3D11Device* pDevice, ID3D11DeviceContext* pContext)
{
	HRESULT hr = S_OK;

	V_RETURN(GameObject::OnD3D11CreateDevice(pDevice, pContext, L"Atmosphere.fx", TechStr, VarStr, ShaderResourceVarStr));
	V_RETURN(D3DX11CreateShaderResourceViewFromFile(pDevice, L"Texture/earth.tiff", nullptr, nullptr, &pEarthGroundSRV.p, nullptr));

	SetTextureSize();

	m_pCloud->OnD3D11CreateDevice(pDevice, pContext);

	return hr;
}


void Atmosphere::SetTextureSize()
{
	VarMap["SCREEN_WIDTH"]->SetRawValue(&screen_width, 0, sizeof(int));
	VarMap["SCREEN_HEIGHT"]->SetRawValue(&screen_height, 0, sizeof(int));

	VarMap["TRANSMITTANCE_TEXTURE_WIDTH"]->SetRawValue(&TRANSMITTANCE_TEXTURE_WIDTH, 0, sizeof(int));
	VarMap["TRANSMITTANCE_TEXTURE_HEIGHT"]->SetRawValue(&TRANSMITTANCE_TEXTURE_WIDTH, 0, sizeof(int));

	VarMap["SCATTERING_TEXTURE_R_SIZE"]->SetRawValue(&SCATTERING_TEXTURE_R_SIZE, 0, sizeof(int));
	VarMap["SCATTERING_TEXTURE_MU_SIZE"]->SetRawValue(&SCATTERING_TEXTURE_MU_SIZE, 0, sizeof(int));
	VarMap["SCATTERING_TEXTURE_MU_S_SIZE"]->SetRawValue(&SCATTERING_TEXTURE_MU_S_SIZE, 0, sizeof(int));
	VarMap["SCATTERING_TEXTURE_NU_SIZE"]->SetRawValue(&SCATTERING_TEXTURE_NU_SIZE, 0, sizeof(int));

	VarMap["SCATTERING_TEXTURE_WIDTH"]->SetRawValue(&SCATTERING_TEXTURE_WIDTH, 0, sizeof(int));
	VarMap["SCATTERING_TEXTURE_HEIGHT"]->SetRawValue(&SCATTERING_TEXTURE_HEIGHT, 0, sizeof(int));
	VarMap["SCATTERING_TEXTURE_DEPTH"]->SetRawValue(&SCATTERING_TEXTURE_DEPTH, 0, sizeof(int));

	VarMap["IRRADIANCE_TEXTURE_WIDTH"]->SetRawValue(&IRRADIANCE_TEXTURE_WIDTH, 0, sizeof(int));
	VarMap["IRRADIANCE_TEXTURE_HEIGHT"]->SetRawValue(&IRRADIANCE_TEXTURE_HEIGHT, 0, sizeof(int));
}


HRESULT Atmosphere::PreCompute(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11RenderTargetView* pRTV)
{
	HRESULT hr = S_OK;

	if(!IsPreComputed)
	{
		V_RETURN(PreComputeTransmittanceTex2D(pDevice, pContext));
		V_RETURN(PreComputeDirectIrradianceTex2D(pDevice, pContext));
		V_RETURN(PreComputeSingleSctrTex3D(pDevice, pContext));

		for (int scatter_order = 2; scatter_order <= scatter_order_num; ++scatter_order)
		{
			V_RETURN(PreComputeInDirectIrradianceTex2D(pDevice, pContext, scatter_order - 1));
			V_RETURN(PreComputeMultiSctrTex3D(pDevice, pContext, scatter_order));
		}
		pContext->CopyResource(pMultiScatterTex3D, pSingleScatterTex3D);

		m_pCloud->PreCompute(pDevice, pContext, pRTV);
		
		IsPreComputed = true;
	}
	return hr;
}


void Atmosphere::Render(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11RenderTargetView* pRTV,ID3D11ShaderResourceView* depthSRV)
{
	//SetView(2.7e6, 0.81, 0.0, 1.57, 2.0, 10.0);
	ID3DX11EffectTechnique* activeTech = TechMap["DrawGroundAndSkyTech"];
	MiscDynamicParams misc;

	D3DXVECTOR3 EyePos = *m_FirstPersonCamera.GetEyePt();
	D3DXVECTOR3 LookAt = *m_FirstPersonCamera.GetLookAtPt();

	cameraParams.f3CameraPos = *m_FirstPersonCamera.GetEyePt();
	D3DXMATRIX View = *m_FirstPersonCamera.GetViewMatrix();
	float det = D3DXMatrixDeterminant(&View);
	D3DXMATRIX InvView;
	D3DXMatrixInverse(&InvView, &det,&View);

	D3DXMATRIX InvViewProj;
	InvViewProj = InvProj * InvView;
	D3DXMatrixTranspose(&InvViewProj, &InvViewProj);
	cameraParams.InvViewProj = InvViewProj;
	
	float cos_sun_z = cos(sun_zenith_angle_radians);
	float sin_sun_z = sin(sun_zenith_angle_radians);
	float cos_sun_a = cos(sun_azimuth_angle_radians);
	float sin_sun_a = sin(sun_azimuth_angle_radians);
	lightParams.f3LightDir = D3DXVECTOR3(sin_sun_z*cos_sun_a, cos_sun_z, sin_sun_z*sin_sun_a);
	
	cameraParams.f3CameraDir = LookAt - EyePos;
	D3DXVec3Normalize(&cameraParams.f3CameraDir, &cameraParams.f3CameraDir);

	misc.scatter_order = 1;

	VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));
	VarMap["misc"]->SetRawValue(&misc, 0, sizeof(MiscDynamicParams));
	VarMap["camera"]->SetRawValue(&cameraParams, 0, sizeof(CameraParams));
	VarMap["light"]->SetRawValue(&lightParams, 0, sizeof(LightParams));
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

	pContext->OMSetRenderTargets(1, &pRTV, nullptr);
	RenderQuad(pContext, activeTech, screen_width, screen_height);
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

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeTransmittanceTex2DTech"];

	VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

	ID3D11RenderTargetView* pRTVs[] =
	{
		pTransmittanceRTV.p,
	};
	pContext->OMSetRenderTargets(ARRAYSIZE(pRTVs), pRTVs, nullptr);

	RenderQuad(pContext, activeTech, TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
	//pContext->OMSetRenderTargets(0, nullptr, nullptr);
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
		&pOpticalLengthTex2D.p, &pOpticalLengthSRV.p, &pOpticalLengthRTV.p));

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeOpticalLengthTex2DTech"];

	VarMap["atmosphere"]->SetRawValue(&atmosphereParams, 0, sizeof(AtmosphereParameters));

	ID3D11RenderTargetView* pRTVs[] =
	{
		pOpticalLengthRTV.p,
	};
	pContext->OMSetRenderTargets(ARRAYSIZE(pRTVs), pRTVs, nullptr);

	RenderQuad(pContext, activeTech, TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
	//pContext->OMSetRenderTargets(0, nullptr, nullptr);
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
							&pDirectIrradianceTex2D.p, &pDirectIrradianceSRV.p, &pDirectIrradianceRTV.p));

	ID3DX11EffectTechnique* activeTech = TechMap["ComputeDirectIrradiance2DTech"];

	ID3D11RenderTargetView* pRTVs[] =
	{
		pDirectIrradianceRTV.p,
	};
	pContext->OMSetRenderTargets(ARRAYSIZE(pRTVs), pRTVs, nullptr);

	ShaderResourceVarMap["g_tex2DTransmittanceLUT"]->SetResource(pTransmittanceSRV);

	RenderQuad(pContext, activeTech, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
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
#if CREATE_TEXTURE_DDS_TEST
	V_RETURN(D3DX11SaveTextureToFile(pContext, pSingleScatterTex3D, D3DX11_IFF_DDS, L"Texture/SingleScatter/SingleScatter3D.dds"));
	V_RETURN(D3DX11SaveTextureToFile(pContext, pSingleScatterCombinedTex3D, D3DX11_IFF_DDS, L"Texture/SingleScatter/SingleScatterCombined3D.dds"));
#endif
	return hr;
}


HRESULT Atmosphere::PreComputeInDirectIrradianceTex2D(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, int scatter_order)
{
	HRESULT hr = S_OK;

	DXGI_FORMAT format = DXGI_FORMAT_R16G16B16A16_FLOAT;

	CComPtr<ID3D11RenderTargetView>	pIndirectIrradianceRTV;
	pIndirectIrradianceTex2D.Release();
	pIndirectIrradianceSRV.Release();
	V_RETURN(CreateTexture2D(pDevice, pContext, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT,
		format,&pIndirectIrradianceTex2D.p, &pIndirectIrradianceSRV.p, &pIndirectIrradianceRTV.p));

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
#if CREATE_TEXTURE_DDS_TEST
	std::wstringstream ss;
	ss << "Texture/Irradiance/IndirectIrradiance_" << scatter_order << ".dds";
	V_RETURN(D3DX11SaveTextureToFile(pContext, pIndirectIrradianceTex2D, D3DX11_IFF_DDS, ss.str().c_str()));
#endif
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
	m_FirstPersonCamera.SetViewParams(&EyePos, &LookAt);

}


void Atmosphere::MsgProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	m_FirstPersonCamera.HandleMessages(hWnd, uMsg, wParam, lParam);

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


void Atmosphere::OnFrameMove(double fTime, float fElapsedTime)
{
	m_FirstPersonCamera.FrameMove(fElapsedTime);
	if(GetAsyncKeyState('I') & 0x8000)
	{
		sun_zenith_angle_radians -= fElapsedTime / 10;
	}
	if (GetAsyncKeyState('K') & 0x8000)
	{
		sun_zenith_angle_radians += fElapsedTime / 10;
	}
	if (GetAsyncKeyState('J') & 0x8000)
	{
		sun_azimuth_angle_radians += fElapsedTime / 10;
	}
	if (GetAsyncKeyState('L') & 0x8000)
	{
		sun_azimuth_angle_radians -= fElapsedTime / 10;
	}
}


void Atmosphere::Resize(int screen_width, int screen_height, float fFOV, float fAspect, float fNear, float fFar)
{
	this->screen_width = screen_width;
	this->screen_height = screen_height;

	VarMap["SCREEN_WIDTH"]->SetRawValue(&screen_width, 0, sizeof(int));
	VarMap["SCREEN_HEIGHT"]->SetRawValue(&screen_height, 0, sizeof(int));

	m_FirstPersonCamera.SetProjParams(fFOV, fAspect, fNear, fFar);
	D3DXMATRIX Proj = *m_FirstPersonCamera.GetProjMatrix();
	float det = D3DXMatrixDeterminant(&Proj);
	D3DXMatrixInverse(&InvProj, &det, &Proj);
	//InvProj = D3DXMATRIX(fFOV * fAspect, 0.0, 0.0, 0.0,
	//	0.0, fFOV, 0.0, 0.0,
	//	0.0, 0.0, 0.0, -1,
	//	0.0, 0.0, 1.0, 1);
}
