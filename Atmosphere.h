#pragma once

#include "DXUT.h"
#include "Common.h"
#include <string>
#include <map>
#include "d3dx11effect.h"

class Atmosphere
{
public:
	Atmosphere(void);
	~Atmosphere(void);

	void Initialize();
	void Release();

	HRESULT OnD3D11CreateDevice(ID3D11Device* pd3dDevice);
	
	void Render();
private:
	float Kr;
	float Km;
	float fKr4PI;
	float fKm4PI;
	float ESun;
	float fKrESun;
	float fKmESun;
	float g;
	float g2;
	float fInnerRadius;
	float fOutRadius;
	float fScale;
	float fCameraHeight;
	float fCameraHeight2;
	D3DXVECTOR3 fWavelength;
	D3DXVECTOR3 fInvWavelength4;
	D3DXVECTOR3 SunPos;
	
	D3DXMATRIX mView;
	D3DXMATRIX mProj;

	ID3DX11Effect*									pAtmosphereEffect;

	std::map<std::string,ID3DX11EffectTechnique*>	AtmosphereTechMap;	
	std::map<std::string,ID3DX11EffectMatrixVariable*>	MatrixVarMap;
	std::map<std::string,ID3DX11EffectVectorVariable*>	VectorVarMap;
	std::map<std::string,ID3DX11EffectScalarVariable*>	ScalarVarMap;
	std::map<std::string,ID3DX11EffectShaderResourceVariable*>	ShaderResourceVarMap;

};

