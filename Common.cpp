#include "DXUT.h"
#include "Common.h"

HRESULT GameObject::OnD3D11CreateDevice(ID3D11Device* pDevice, ID3D11DeviceContext* pContext, 
				std::wstring EffectStr, const std::vector<std::string>& TechStr, 
				const std::vector<std::string>& VarStr, const std::vector<std::string>& ShaderResourceVarStr)
{
	HRESULT hr = S_OK;

	V_RETURN(CompileEffectFromFile(pDevice, &pEffect, const_cast<wchar_t*>(EffectStr.c_str())));

	for (int i = 0; i < TechStr.size(); i++)
	{
		CComPtr<ID3DX11EffectTechnique> pTech = pEffect->GetTechniqueByName(TechStr[i].c_str());
		TechMap.emplace(TechStr[i], pTech);
	}

	//for (int i = 0; i < MatrixVarStr.size(); i++)
	//{
	//	CComPtr<ID3DX11EffectMatrixVariable> pMatrixVar = pEffect->GetVariableByName(MatrixVarStr[i].c_str())->AsMatrix();
	//	MatrixVarMap.emplace(MatrixVarStr[i], pMatrixVar);
	//}

	//for (int i = 0; i < VectorVarStr.size(); i++)
	//{
	//	CComPtr<ID3DX11EffectVectorVariable> pVectorVar = pAtmosphereEffect->GetVariableByName(VectorVarStr[i].c_str())->AsVector();
	//	VectorVarMap.emplace(VectorVarStr[i], pVectorVar);
	//}

	//for (int i = 0; i < ScalarVarStr.size(); i++)
	//{
	//	CComPtr<ID3DX11EffectScalarVariable> pScalarVar = pAtmosphereEffect->GetVariableByName(ScalarVarStr[i].c_str())->AsScalar();
	//	ScalarVarMap.emplace(ScalarVarStr[i], pScalarVar);
	//}

	for (int i = 0; i < VarStr.size(); i++)
	{
		CComPtr<ID3DX11EffectVariable> pVar = pEffect->GetVariableByName(VarStr[i].c_str());
		VarMap.emplace(VarStr[i], pVar);
	}

	for (int i = 0; i < ShaderResourceVarStr.size(); i++)
	{
		CComPtr<ID3DX11EffectShaderResourceVariable> pShaderResourceVar = pEffect->GetVariableByName(ShaderResourceVarStr[i].c_str())->AsShaderResource();
		ShaderResourceVarMap.emplace(ShaderResourceVarStr[i], pShaderResourceVar);
	}
}

void GameObject::Release()
{
	MapRelease(TechMap);
	MapRelease(MatrixVarMap);
	MapRelease(VectorVarMap);
	MapRelease(ScalarVarMap);
	MapRelease(VarMap);
	MapRelease(ShaderResourceVarMap);
	pEffect.Release();
}