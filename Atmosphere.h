#pragma once

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include "DXUT.h"
#include <atlcomcli.h>
#include <unordered_map>
#include <string>
#include <vector>

#define CREATE_TEXTURE_DDS_TEST 1

struct DensityProfileLayer
{
	float exp_term;
	float exp_scale;
	float linear_term;
	float const_term;
};

struct AtmosphereParameters
{
	D3DXVECTOR3 solar_irradiance;
	float bottom_radius;

	D3DXVECTOR3 rayleigh_scattering;
	float top_radius;

	D3DXVECTOR3 mie_scattering;
	float mie_g;

	D3DXVECTOR3 mie_extinction;
	float ground_albedo;

	D3DXVECTOR3 absorption_extinction;
	float ozone_width;

	float sun_angular_radius;
	float mu_s_min;
	float padding[2];

	DensityProfileLayer rayleigh_density;
	DensityProfileLayer mie_density;
	DensityProfileLayer ozone_density[2];
};

struct MiscDynamicParams
{
	D3DXVECTOR2 f2WQ;
	int scatter_order;
	float exposure;

	D3DXVECTOR3 f3CameraPos;
	D3DXVECTOR3 f3EarthCenter;
	D3DXVECTOR3 f3SunDir;
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
	int scatter_order_num;

	bool IsPreComputed = false;

	float view_distance_meters = 9000.f;
	float view_zenith_angle_radians = 1.47f;
	float view_azimuth_angle_radians = -0.1f;
	float sun_zenith_angle_radians = 1.3f;
	float sun_azimuth_angle_radians = 2.9f;
	float exposure = 10.f;


	AtmosphereParameters atmosphereParams;

	HRESULT PreComputeTransmittanceTex2D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT PreComputeDirectIrradianceTex2D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT PreComputeSingleSctrTex3D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT PreComputeInDirectIrradianceTex2D(ID3D11Device*, ID3D11DeviceContext*,int);
	HRESULT PreComputeMultiSctrTex3D(ID3D11Device*, ID3D11DeviceContext*,int);

	HRESULT PreComputeSingleSctrTex3D_Test(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT PreComputeMultiSctrTex3D_Test(ID3D11Device*, ID3D11DeviceContext*, int);

	const int TRANSMITTANCE_TEXTURE_WIDTH = 256;    //mu
	const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;    //r

	const int SCATTERING_TEXTURE_R_SIZE = 32;
	const int SCATTERING_TEXTURE_MU_SIZE = 128;
	const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
	const int SCATTERING_TEXTURE_NU_SIZE = 8;

	const int SCATTERING_TEXTURE_WIDTH = SCATTERING_TEXTURE_R_SIZE;
	const int SCATTERING_TEXTURE_HEIGHT = SCATTERING_TEXTURE_MU_SIZE;
	const int SCATTERING_TEXTURE_DEPTH = SCATTERING_TEXTURE_MU_S_SIZE*SCATTERING_TEXTURE_NU_SIZE;

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
	
	CComPtr<ID3D11Texture3D>							pSingleScatterTex3D;
	CComPtr<ID3D11ShaderResourceView>					pSingleScatterSRV;

	CComPtr<ID3D11Texture3D>							pSingleScatterCombinedTex3D;
	CComPtr<ID3D11ShaderResourceView>					pSingleScatterCombinedSRV;
	
	CComPtr<ID3D11Texture3D>							pSingleScatterMieTex3D;
	CComPtr<ID3D11ShaderResourceView>					pSingleScatterMieSRV;

	CComPtr<ID3D11Texture2D>							pDirectIrradianceTex2D;
	CComPtr<ID3D11ShaderResourceView>					pDirectIrradianceSRV;

	CComPtr<ID3D11Texture2D>							pIndirectIrradianceTex2D;
	CComPtr<ID3D11ShaderResourceView>					pIndirectIrradianceSRV;

	CComPtr<ID3D11Texture3D>							pMultiScatterTex3D;
	CComPtr<ID3D11ShaderResourceView>					pMultiScatterSRV;

	CComPtr<ID3D11Texture3D>							pMultiScatterCombinedTex3D;
	CComPtr<ID3D11ShaderResourceView>					pMultiScatterCombinedSRV;

	std::vector<std::string> TechStr
	{
		"ComputeTransmittanceTex2DTech",
		"ComputeDirectIrradiance2DTech",
		"ComputeSingleScatterTex3DTech",
		"ComputeIndirectIrradiance2DTech",
		"ComputeMultiScatterTex3DTech"
	};

	std::vector<std::string> MatrixVarStr
	{
		"World",
		"WorldViewProj",
		"WorldInvTranspose"
	};

	std::vector<std::string> VectorVarStr
	{
		"solar_irradiance",
		"rayleigh_scattering",
		"mie_scattering",
		"mie_extinction",
		"absorption_extinction",
		"v3CameraPos",
		"v3LightPos"
	};

	std::vector<std::string> ScalarVarStr
	{
		"bottom_radius",
		"top_radius",
		"mie_g",
		"ground_albedo",
		"ozone_width"
	};

	std::vector<std::string> VarStr
	{
		"atmosphere",
		"rayleigh_density",
		"mie_density",
		"ozone_density",
		"misc"
	};

	std::vector<std::string> ShaderResourceVarStr
	{
		"g_tex2DTransmittanceLUT",
		"g_tex2DDirectIrradianceLUT",
		"g_tex3DSingleScatteringLUT",
		"g_tex2DDirectIrradianceLUT",
		"g_tex3DMultiScatteringLUT",

		"g_tex3DSingleMieScatteringLUT",
		"g_tex3DSingleScatteringCombinedLUT",
		"g_tex3DMultiScatteringCombinedLUT",

		"g_tex2DEarthGround"
	};
};

#endif