#pragma once

#include "DXUT.h"
#include <string>
#include "d3dx11effect.h"
#include "GeometryGenerator.h"
#include <unordered_map>
#include <vector>
#include <atlcomcli.h>

class Atmosphere
{
public:
	Atmosphere(void);
	~Atmosphere(void);

	void Initialize(ID3D11Device*,ID3D11DeviceContext*);
	void Release();

	HRESULT OnD3D11CreateDevice();
	void ResetInputView(GeometryGenerator::MeshData&);

	void PreCompute();

	void Render(D3DXMATRIX&,D3DXVECTOR3&,GeometryGenerator::MeshData&);
private:
	void PreComputeTransmittance();

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

	float Kr;
	float Km;
	float fKr4PI;
	float fKm4PI;
	float ESun;
	float fKrESun;
	float fKmESun;
	float fMieG;
	float fMieG2;
	float fInnerRadius;
	float fInnerRadius2;
	float fOutRadius;
	float fOutRadius2;
	float fScale;
	float fScaleOverScaleDepth;
	float fScaleDepth;
	float fInvScaleDepth;
	float fCameraHeight;
	float fCameraHeight2;
	D3DXVECTOR3 fWavelength;
	D3DXVECTOR3 fInvWavelength4;
	D3DXVECTOR3 SunPos;
	
	ID3DX11Effect*									pAtmosphereEffect;

	std::unordered_map<std::string,ID3DX11EffectTechnique*>					AtmosphereTechMap;	
	std::unordered_map<std::string,ID3DX11EffectMatrixVariable*>			MatrixVarMap;
	std::unordered_map<std::string,ID3DX11EffectVectorVariable*>			VectorVarMap;
	std::unordered_map<std::string,ID3DX11EffectScalarVariable*>			ScalarVarMap;
	std::unordered_map<std::string,ID3DX11EffectShaderResourceVariable*>	ShaderResourceVarMap;

	ID3D11InputLayout*								pInputLayout;
	ID3D11Buffer*									pGroundVB;
	ID3D11Buffer*									pGroundIB;
	ID3D11Buffer*									pAtmosphereVB;
	ID3D11Buffer*									pAtmosphereIB;
	ID3D11ShaderResourceView*						pGroundSRV;

	ID3D11Texture2D*								pTransmittanceTex2D;
	ID3D11ShaderResourceView*						pTransmittanceSRV;
	ID3D11RenderTargetView*							pTransmittanceRTV;

	ID3D11Device* pd3dDevice;
	ID3D11DeviceContext* pd3dImmediateContext;
	UINT groundIndexNum;

	template<typename T1, typename T2>
	void MapRelease(std::unordered_map<T1,T2>& m);

	void ResetShaderParam();

	std::vector<std::string> TechStr
	{
		"GroundFromAtmosphere",
		"GroundFromSpace",
		"SkyFromAtmosphere",
		"SkyFromSpace",
		"SpaceFromAtmosphere",
		"SpaceFromSpace",
		"Test",
		"mSkyfromAtmosphere",
		"PreComputeTransmittanceTexture"
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
		"v3LightPos",
		"v3InvWavelength"
	};

	std::vector<std::string> ScalarVarStr
	{
		"fCameraHeight",
		"fCameraHeight2",
		"fOuterRadius",
		"fOuterRadius2",
		"fInnerRadius",
		"fInnerRadius2",
		"fKrESun",
		"fKmESun",
		"fKr4PI",
		"fKm4PI",
		"fScale",
		"fScaleOverScaleDepth",
		"fScaleDepth",
		"fInvScaleDepth",
		"fMieG",
		"fMieG2"
	};

	std::vector<std::string> ShaderResourceVarStr
	{
		"GroundMap"
	};
};

