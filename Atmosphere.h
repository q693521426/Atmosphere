#pragma once

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include "DXUT.h"
#include <atlcomcli.h>
#include <unordered_map>
#include <string>
#include <vector>

struct DensityProfileLayer
{
	float width;
	float exp_term;
	float exp_scale;
	float linear_term;
	float const_term;
};

struct DensityProfile
{
	DensityProfileLayer layer[2];
};

struct AtmosphereParameters
{
	float bottom_radius;
	float top_radius;

	DensityProfile rayleigh_density;
	float rayleigh_scattering[3];

	DensityProfile mie_density;
	float mie_scattering[3];
};

class Atmosphere
{
public:
	Atmosphere();
	~Atmosphere();
	
	void Initialize();
	void Release();

	HRESULT OnD3D11CreateDevice(ID3D11Device*, ID3D11DeviceContext*);

	void Render(ID3D11Device*, ID3D11DeviceContext*, ID3D11RenderTargetView*);
private:
	bool IsPreComputed = false;
	AtmosphereParameters AtmosphereParams;
	
	HRESULT PreComputeTransmittanceTex2D(ID3D11Device*, ID3D11DeviceContext*, ID3D11RenderTargetView*);
	HRESULT PreComputeSingleSctrTex3D(ID3D11Device*, ID3D11DeviceContext*);

	const int TRANSMITTANCE_TEXTURE_WIDTH = 256;    //mu
	const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;    //r

	const int SCATTERING_TEXTURE_R_SIZE = 32;
	const int SCATTERING_TEXTURE_MU_SIZE = 128;
	const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
	const int SCATTERING_TEXTURE_NU_SIZE = 8;

	const int SCATTERING_TEXTURE_WIDTH =
		SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE;
	const int SCATTERING_TEXTURE_HEIGHT = SCATTERING_TEXTURE_MU_SIZE;
	const int SCATTERING_TEXTURE_DEPTH = SCATTERING_TEXTURE_R_SIZE;

	const int IRRADIANCE_TEXTURE_WIDTH = 64;
	const int IRRADIANCE_TEXTURE_HEIGHT = 16;

	// The conversion factor between watts and lumens.
	const double MAX_LUMINOUS_EFFICACY = 683.0;

	CComPtr<ID3DX11Effect>		pAtmosphereEffect;

	std::unordered_map<std::string, CComPtr<ID3DX11EffectTechnique>>				AtmosphereTechMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectMatrixVariable>>			MatrixVarMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectVectorVariable>>			VectorVarMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectScalarVariable>>			ScalarVarMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectVariable>>					VarMap;
	std::unordered_map<std::string, CComPtr<ID3DX11EffectShaderResourceVariable>>	ShaderResourceVarMap;

	CComPtr<ID3D11Texture2D>							pTransmittanceTex2D;
	CComPtr<ID3D11ShaderResourceView>					pTransmittanceSRV;
	CComPtr<ID3D11RenderTargetView>						pTransmittanceRTV;
	
	std::vector<std::string> TechStr
	{
		"ComputeTransmittanceTex2DTech"
	};

	std::vector<std::string> MatrixVarStr
	{
		"World",
		"WorldViewProj",
		"WorldInvTranspose"
	};

	std::vector<std::string> VectorVarStr
	{
		"v3CameraPos",
		"v3LightPos"
	};

	std::vector<std::string> ScalarVarStr
	{
	};

	std::vector<std::string> VarStr
	{
		"AtmosphereParams"
	};

	std::vector<std::string> ShaderResourceVarStr
	{
		"g_tex2DTransmittanceLUT",
		"g_tex3DSingleScatteringLUT"
	};
};

#endif