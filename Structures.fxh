
#ifndef STRUCTURES_FXH_
#define STRUCTURES_FXH_

#ifdef __cplusplus

#	define float2 D3DXVECTOR2
#	define float3 D3DXVECTOR3
#	define float4 D3DXVECTOR4
#	define matrix D3DXMATRIX

#	define uint UINT
#	define uint4 UINT4

#else

#endif

#ifdef __cplusplus
#	define CHECK_STRUCT_ALIGNMENT(s) static_assert(sizeof(s)%16 == 0,"sizeof("#s") is not multiple of 16")
#else
#	define CHECK_STRUCT_ALIGNMENT(s)
#endif

struct DensityProfileLayer
{
	float exp_term;
	float exp_scale;
	float linear_term;
	float const_term;
};
CHECK_STRUCT_ALIGNMENT(DensityProfileLayer);

struct AtmosphereParams
{
	float3 solar_irradiance;
	float bottom_radius;

	float3 rayleigh_scattering;
	float top_radius;

	float3 mie_scattering;
	float mie_g;

	float3 mie_extinction;
	float ground_albedo;

	float3 absorption_extinction;
	float ozone_width;

	float sun_angular_radius;
	float mu_s_min;
	float nu_power;
	float radius_scale;

	DensityProfileLayer rayleigh_density;
	DensityProfileLayer mie_density;
	DensityProfileLayer ozone_density[2];
};
CHECK_STRUCT_ALIGNMENT(AtmosphereParams);

struct MiscDynamicParams
{
	float2 f2WQ;
	float scatter_order;
	uint uiMinMaxLevelMax;
#ifdef __cplusplus
	UINT uiSrcMinMaxOffsetX;
	UINT uiSrcMinMaxOffsetY;
	UINT uiDstMinMaxOffsetX;
	UINT uiDstMinMaxOffsetY;
#else
	uint4 ui4SrcDstMinMaxOffset;
#endif
	float fEnableLightShaft;
	float fIsLightInSpaceCorrect;
	float fTime;
	float padding;
};
CHECK_STRUCT_ALIGNMENT(MiscDynamicParams);

struct CameraParams
{
	float3 f3CameraPos;
	float fNearZ;

	float3 f3CameraDir;
	float fFarZ;

	matrix View;
	matrix Proj;
	matrix ViewProj;
	matrix InvViewProj;
};
CHECK_STRUCT_ALIGNMENT(CameraParams);

struct LightParams
{
	float3 f3LightDir;
	float padding;

	float4 f4LightScreenPos;

	matrix View;
	matrix Proj;
	matrix ViewProj;
	matrix InvViewProj;
};
CHECK_STRUCT_ALIGNMENT(LightParams);

#endif