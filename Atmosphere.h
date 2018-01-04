#pragma once

#include "DXUT.h"
#include <string>
#include "d3dx11effect.h"
#include "GeometryGenerator.h"
#include <unordered_map>
#include <vector>

class Atmosphere
{
public:
	Atmosphere(void);
	~Atmosphere(void);

	void Initialize(ID3D11Device*,ID3D11DeviceContext*);
	void Release();

	HRESULT OnD3D11CreateDevice();
	void ResetInputView(GeometryGenerator::MeshData&);
	
	void Render(D3DXMATRIX&,D3DXVECTOR3&,GeometryGenerator::MeshData&);
private:
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

	ID3D11Device* pd3dDevice;
	ID3D11DeviceContext* pd3dImmediateContext;
	UINT groundIndexNum;

	template<typename T1, typename T2>
	void MapRelease(std::unordered_map<T1,T2>& m);

	void ResetShaderParam();

};

