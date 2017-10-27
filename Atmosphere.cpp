#include "DXUT.h"
#include "Atmosphere.h"

char* TechStr[]= 
{
	"GroundFromAtmosphere",
	"GroundFromSpace",
	"SkyFromAtmosphere",
	"SkyFromSpace",
	"SpaceFromAtmosphere",
	"SpaceFromSpace"
};

char* MatrixVarStr[]=
{
	"World",
	"WorldViewProj",
	"WorldInvTranspose"
};

char* VectorVarStr[]=
{
	"v3CameraPos",
	"v3LightPos",
	"v3InvWavelength"
};

char* ScalarVarStr[]=
{
	"fCameraHeight",
	"fCameraHeight2",
	"fOuterRadius",
	"fOuterRadius2",
	"fInnerRadius",
	"fInnerRadius2",
	"fKrESun",
	"fKmESun",
	"fScale",			
	"fScaleOverScaleDepth",
	"fMieG",
	"fMieG2"
};

char* ShaderResourceVarStr[]=
{
	"GroundMap"
};

Atmosphere::Atmosphere(void):pAtmosphereEffect(nullptr)
{
}


Atmosphere::~Atmosphere(void)
{
}

void Atmosphere::Initialize()
{
	Kr = 0.0025f;
	Km = 0.0010f;
	fKr4PI = Kr*4.0f*D3DX_PI;
	fKm4PI = Km*4.0f*D3DX_PI;
	ESun = 20.0f;
	fKrESun = ESun*Kr;
	fKmESun = ESun*Km;
	g = -0.990;
	g2 = g*g;
	fInnerRadius;
	fOutRadius;
	fScale = 1/(fOutRadius - fInnerRadius);
	fCameraHeight;
	fCameraHeight2;
	fWavelength = D3DXVECTOR3(0.650f,0.570f,0.475f);
	fInvWavelength4 = D3DXVECTOR3(1/pow(fWavelength[0],4),1/pow(fWavelength[1],4),1/pow(fWavelength[2],4));
}

void Atmosphere::Release()
{
	SAFE_RELEASE(pAtmosphereEffect);
}

HRESULT Atmosphere::OnD3D11CreateDevice(ID3D11Device* pd3dDevice)
{
	HRESULT hr = S_OK;
	V_RETURN(CompileEffectFromFile(pd3dDevice,pAtmosphereEffect,L"Atmosphere.fx"));
	
	for(int i=0;i<ARRAYSIZE(TechStr);i++)
	{
		auto pTech = pAtmosphereEffect->GetTechniqueByName(TechStr[i]);
		AtmosphereTechMap.emplace(TechStr[i],pTech);
	}
	for(int i=0;i<ARRAYSIZE(MatrixVarStr);i++)
	{
		auto pMatrixVar = pAtmosphereEffect->GetVariableByName(MatrixVarStr[i])->AsMatrix();
		MatrixVarMap.emplace(MatrixVarStr[i],pMatrixVar);
	}
	for(int i=0;i<ARRAYSIZE(VectorVarStr);i++)
	{
		auto pVectorVar = pAtmosphereEffect->GetVariableByName(MatrixVarStr[i])->AsVector();
		VectorVarMap.emplace(VectorVarStr[i],pVectorVar);
	}
	for(int i=0;i<ARRAYSIZE(ScalarVarStr);i++)
	{
		auto pScalarVar = pAtmosphereEffect->GetVariableByName(ScalarVarStr[i])->AsScalar();
		ScalarVarMap.emplace(ScalarVarStr[i],pScalarVar);
	}
	for(int i=0;i<ARRAYSIZE(ShaderResourceVarStr);i++)
	{
		auto pShaderResourceVar = pAtmosphereEffect->GetVariableByName(ShaderResourceVarStr[i])->AsShaderResource();
		ShaderResourceVarMap.emplace(ShaderResourceVarStr[i],pShaderResourceVar);
	}
	return hr;
}

void Atmosphere::Render()
{
		
}