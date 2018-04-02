#pragma once

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include "DXUT.h"
#include "DXUTcamera.h"
#include "Common.h"
#include "Cloud.h"

#define CREATE_TEXTURE_DDS_TEST 0

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
	float nu_power;
	float exposure;

	DensityProfileLayer rayleigh_density;
	DensityProfileLayer mie_density;
	DensityProfileLayer ozone_density[2];
};

class Atmosphere : public GameObject
{
public:
	Atmosphere();
	~Atmosphere();

	void Initialize();
	void Release();

	void SetView(float, float, float, float, float, float);

	void Resize(int, int, float, float, float, float);

	HRESULT OnD3D11CreateDevice(ID3D11Device*, ID3D11DeviceContext*);
	void SetTextureSize();
	void SetCameraParams();
	void SetLightParams();

	HRESULT PreCompute(ID3D11Device*, ID3D11DeviceContext*, ID3D11RenderTargetView*);
	void Render(ID3D11Device*, ID3D11DeviceContext*, ID3D11RenderTargetView*, ID3D11ShaderResourceView* depthSRV=nullptr);

	void MsgProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

	void OnFrameMove(double fTime, float fElapsedTime);

	D3DXVECTOR3 GetSunDir();

	float GetCameraHeight();
private:
	CFirstPersonCamera					m_FirstPersonCamera;
	Cloud*								m_pCloud;

	int scatter_order_num;

	bool IsPreComputed = false;

	float view_distance_meters = 9000.f;
	float view_zenith_angle_radians = 1.47f;
	float view_azimuth_angle_radians = -0.1f;
	float sun_zenith_angle_radians = 1.3f;
	float sun_azimuth_angle_radians = 2.9f;
	float exposure = 10.f;

	D3DXMATRIX InvView, InvProj;
	int screen_width, screen_height;

	AtmosphereParameters atmosphereParams;
	CameraParams cameraParams;
	LightParams lightParams;

	HRESULT PreComputeTransmittanceTex2D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT PreComputeOpticalLengthTex2D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT PreComputeDirectIrradianceTex2D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT PreComputeSingleSctrTex3D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT PreComputeInDirectIrradianceTex2D(ID3D11Device*, ID3D11DeviceContext*,int);
	HRESULT PreComputeMultiSctrTex3D(ID3D11Device*, ID3D11DeviceContext*,int);

	HRESULT ComputeSpaceLinearDepthTex2D(ID3D11Device*, ID3D11DeviceContext*, ID3D11ShaderResourceView*);
	HRESULT ComputeSliceEndTex2D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT ComputeEpipolarCoordTex2D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT RefineSampleLocal(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT ComputeSliceUVOrigDirTex2D(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT Build1DMinMaxMipMap(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT MarkRayMarchSample(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT DoRayMarch(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT InterpolateScatter(ID3D11Device*, ID3D11DeviceContext*);
	HRESULT ApplyInterpolateScatter(ID3D11Device*, ID3D11DeviceContext*);

	int TRANSMITTANCE_TEXTURE_WIDTH = 256;    //mu
	int TRANSMITTANCE_TEXTURE_HEIGHT = 64;    //r

	int SCATTERING_TEXTURE_R_SIZE = 32;
	int SCATTERING_TEXTURE_MU_SIZE = 128;
	int SCATTERING_TEXTURE_MU_S_SIZE = 32;
	int SCATTERING_TEXTURE_NU_SIZE = 8;

	int SCATTERING_TEXTURE_WIDTH = SCATTERING_TEXTURE_R_SIZE;
	int SCATTERING_TEXTURE_HEIGHT = SCATTERING_TEXTURE_MU_SIZE;
	int SCATTERING_TEXTURE_DEPTH = SCATTERING_TEXTURE_MU_S_SIZE*SCATTERING_TEXTURE_NU_SIZE;

	int IRRADIANCE_TEXTURE_WIDTH = 64;
	int IRRADIANCE_TEXTURE_HEIGHT = 16;

	// The conversion factor between watts and lumens.
	double MAX_LUMINOUS_EFFICACY = 683.0;

	int EPIPOLAR_SLICE_NUM = 512;
	int EPIPOLAR_SAMPLE_NUM = 256;

	int RefineSampleCSThreadGroupSize = 128;
	int InterpolationSampleStep = 16;

	CComPtr<ID3D11Texture2D>							pTransmittanceTex2D;
	CComPtr<ID3D11ShaderResourceView>					pTransmittanceSRV;

	CComPtr<ID3D11Texture2D>							pOpticalLengthTex2D;
	CComPtr<ID3D11ShaderResourceView>					pOpticalLengthSRV;
	
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

	CComPtr<ID3D11ShaderResourceView>					pEarthGroundSRV;

	CComPtr<ID3D11Texture2D>							pSpaceLinearDepthTex2D;
	CComPtr<ID3D11ShaderResourceView>					pSpaceLinearDepthSRV;
	
	CComPtr<ID3D11Texture2D>							pSliceEndTex2D;
	CComPtr<ID3D11ShaderResourceView>					pSliceEndSRV;
	
	CComPtr<ID3D11Texture2D>							pEpipolarSampleTex2D;
	CComPtr<ID3D11ShaderResourceView>					pEpipolarSampleSRV;
	CComPtr<ID3D11DepthStencilView>						pEpipolarSampleDSV;

	CComPtr<ID3D11Texture2D>							pEpipolarSampleCamDepthTex2D;
	CComPtr<ID3D11ShaderResourceView>					pEpipolarSampleCamDepthSRV;

	CComPtr<ID3D11Texture2D>							pInterpolationSampleTex2D;
	CComPtr<ID3D11ShaderResourceView>					pInterpolationSampleSRV;
	CComPtr<ID3D11UnorderedAccessView>					pInterpolationSampleUAV;

	CComPtr<ID3D11Texture2D>							pSliceUVOrigDirTex2D;
	CComPtr<ID3D11ShaderResourceView>					pSliceUVOrigDirSRV;

	CComPtr<ID3D11Texture2D>							pScatterTex2D;
	CComPtr<ID3D11ShaderResourceView>					pScatterSRV;

	CComPtr<ID3D11Texture2D>							pInterpolatedScatterTex2D;
	CComPtr<ID3D11ShaderResourceView>					pInterpolatedScatterSRV;

	std::vector<std::string> TechStr
	{
		"ComputeTransmittanceTex2DTech",
		"ComputeOpticalLengthTex2DTech",
		"ComputeDirectIrradiance2DTech",
		"ComputeSingleScatterTex3DTech",
		"ComputeIndirectIrradiance2DTech",
		"ComputeMultiScatterTex3DTech",
		"DrawGroundAndSkyTech",

		"ComputeSpaceLinearDepthTex2DTech",
		"ComputeSliceEndTex2DTech",
		"ComputeEpipolarCoordTex2DTech",
		"RefineSampleTech",
		"ComputeSliceUVOrigDirTex2DTech",
		"Build1DMinMaxMipMapTech",
		"MarkRayMarchSampleTech",
		"DoRayMarchTech",
		"InterpolateScatterTech",
		"ApplyInterpolateScatterTech"
	};

	std::vector<std::string> VarStr
	{
		"atmosphere",
		"misc",
		"camera",
		"light",

		"SCREEN_WIDTH",
		"SCREEN_HEIGHT",

		"TRANSMITTANCE_TEXTURE_WIDTH",
		"TRANSMITTANCE_TEXTURE_HEIGHT",

		"SCATTERING_TEXTURE_R_SIZE",
		"SCATTERING_TEXTURE_MU_SIZE",
		"SCATTERING_TEXTURE_MU_S_SIZE",
		"SCATTERING_TEXTURE_NU_SIZE",

		"SCATTERING_TEXTURE_WIDTH",
		"SCATTERING_TEXTURE_HEIGHT",
		"SCATTERING_TEXTURE_DEPTH",

		"IRRADIANCE_TEXTURE_WIDTH",
		"IRRADIANCE_TEXTURE_HEIGHT",

		"EPIPOLAR_SLICE_NUM",
		"EPIPOLAR_SAMPLE_NUM"
	};

	std::vector<std::string> ShaderResourceVarStr
	{
		"g_tex2DTransmittanceLUT",
		"g_tex2DOpticalLengthLUT",
		"g_tex2DDirectIrradianceLUT",
		"g_tex2DIndirectIrradianceLUT",
		"g_tex3DSingleScatteringLUT",
		"g_tex3DMultiScatteringLUT",

		"g_tex3DSingleMieScatteringLUT",
		"g_tex3DSingleScatteringCombinedLUT",
		"g_tex3DMultiScatteringCombinedLUT",

		"g_tex2DEarthGround",

		"g_tex2DSpaceDepth",
		"g_tex2DSpaceLinearDepth",
		"g_tex2DSliceEnd",
		"g_tex2DEpipolarSample",
		"g_tex2DEpipolarSampleCamDepth",
		"g_tex2DInterpolationSample",
		"g_tex2DSliceUVOrigDir",
		"g_tex2DScatter",
		"g_tex2DInterpolatedScatter"
	};
};

#endif