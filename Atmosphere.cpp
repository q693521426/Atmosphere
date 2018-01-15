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

	std::vector<double> solar_irradiance;
	std::vector<double> rayleigh_scatter;
	std::vector<double> mie_scatter;
	std::vector<double> mie_extinction;
	std::vector<double> absorption_extinction;	//ozone

	auto lerp = [](double x,double y,double s)
	{
		return x*(1 - s) + y*s;
	};

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
					kMieAngstromBeta / kMieScaleHeight * pow(lambda_y_um, -kMieAngstromAlpha), s);

		solar_irradiance.push_back(
			lerp(kSolarIrradiance[index], kSolarIrradiance[index + 1], s));
		rayleigh_scatter.push_back(
			lerp(kRayleigh * pow(lambda_x_um, -4), kRayleigh *pow(lambda_y_um, -4), s));
		mie_scatter.push_back(mie);
		mie_extinction.push_back(mie*kMieSingleScatteringAlbedo);
		absorption_extinction.push_back(
			lerp(kMaxOzoneNumberDensity * kOzoneCrossSection[index],
				kMaxOzoneNumberDensity * kOzoneCrossSection[index+1],s));

	}
	
	AtmosphereParams.bottom_radius = 6360.f;
	AtmosphereParams.top_radius = 6420.f;
	std::move(rayleigh_scatter.begin(), rayleigh_scatter.end(), 
				AtmosphereParams.rayleigh_scattering);
	AtmosphereParams.rayleigh_density = DensityProfileLayer
	{
		0.f, 1.f, -1.0 / kMieScaleHeight * 1000.0, 0.f, 0.f
	};
	std::move(mie_scatter.begin(), mie_scatter.end(), 
				AtmosphereParams.mie_scattering);
	AtmosphereParams.mie_density = DensityProfileLayer
	{
		0.f, 1.f, -1.0 / kRayleighScaleHeight * 1000.0, 0.f, 0.f
	};
	AtmosphereParams.mie_g = 0.8;
	std::move(absorption_extinction.begin(), absorption_extinction.end(),
				AtmosphereParams.absorption_extinction);
	AtmosphereParams.ozone_density[0] = DensityProfileLayer
	{
		25.0, 0.0, 0.0, 1.0 / 15.0, -2.0 / 3.0
	};
	AtmosphereParams.ozone_density[1] = DensityProfileLayer
	{
		0.0, 0.0, 0.0, -1.0 / 15.0, 8.0 / 3.0
	};
	AtmosphereParams.ground_albedo = 0.1f;
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


void Atmosphere::Render(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, ID3D11RenderTargetView* pRTV)
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