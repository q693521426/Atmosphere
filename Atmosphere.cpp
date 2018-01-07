#include "DXUT.h"
#include "Atmosphere.h"
#include "D3dBuffurDesc.h"
#include "RenderStates.h"
#include <algorithm>
#include <iostream>

Atmosphere::Atmosphere(void) :pAtmosphereEffect(nullptr), pInputLayout(nullptr),
pGroundVB(nullptr), pGroundIB(nullptr), pAtmosphereVB(nullptr),
pAtmosphereIB(nullptr), pGroundSRV(nullptr),
pd3dDevice(nullptr), pd3dImmediateContext(nullptr)
{
}


Atmosphere::~Atmosphere(void)
{
}

void Atmosphere::Initialize(ID3D11Device* d3dDevice, ID3D11DeviceContext* d3dImmediateContext)
{
	Kr = 0.0025f;
	Km = 0.0010f;
	fKr4PI = Kr*4.0f*D3DX_PI;
	fKm4PI = Km*4.0f*D3DX_PI;
	ESun = 20.0f;
	fKrESun = ESun*Kr;
	fKmESun = ESun*Km;
	fMieG = -0.990;
	fMieG2 = fMieG*fMieG;
	fInnerRadius = 6357.0f;
	fInnerRadius2 = fInnerRadius*fInnerRadius;
	fOutRadius = 6357.0f + 100.0f;
	fOutRadius2 = fOutRadius*fOutRadius;
	fScale = 1 / (fOutRadius - fInnerRadius);
	fScaleDepth = 0.25f;
	fInvScaleDepth = 1 / fScaleDepth;
	fScaleOverScaleDepth = fScale / fScaleDepth;
	fCameraHeight = 6357.0f + 25.0f;
	fCameraHeight2 = fCameraHeight * fCameraHeight;
	fWavelength = D3DXVECTOR3(0.650f, 0.570f, 0.475f);
	fInvWavelength4 = D3DXVECTOR3(1 / pow(fWavelength[0], 4), 1 / pow(fWavelength[1], 4), 1 / pow(fWavelength[2], 4));
	SunPos = D3DXVECTOR3(0, 1, 0);

	pd3dDevice = d3dDevice;
	pd3dImmediateContext = d3dImmediateContext;
}
 
template<typename T1, typename T2>
void Atmosphere::MapRelease(std::unordered_map<T1, T2>& m)
{
	std::for_each(m.begin(), m.end(), [](std::pair<const T1, T2>& v) { SAFE_RELEASE(v.second)});
	m.clear();
}

void Atmosphere::Release()
{
	pd3dDevice = nullptr;
	pd3dImmediateContext = nullptr;
	MapRelease(AtmosphereTechMap);
	MapRelease(MatrixVarMap);
	MapRelease(VectorVarMap);
	MapRelease(ScalarVarMap);
	MapRelease(ShaderResourceVarMap);
	SAFE_RELEASE(pAtmosphereEffect);
	SAFE_RELEASE(pInputLayout);
	SAFE_RELEASE(pGroundVB);
	SAFE_RELEASE(pGroundIB);
	SAFE_RELEASE(pAtmosphereVB);
	SAFE_RELEASE(pAtmosphereIB);
	SAFE_RELEASE(pGroundSRV);
}

HRESULT Atmosphere::OnD3D11CreateDevice()
{
	HRESULT hr = S_OK;
	V_RETURN(CompileEffectFromFile(pd3dDevice, &pAtmosphereEffect, L"Atmosphere.fx"));

	for (int i = 0; i < TechStr.size(); i++)
	{
		auto pTech = pAtmosphereEffect->GetTechniqueByName(TechStr[i].c_str());
		AtmosphereTechMap.emplace(TechStr[i], pTech);
	}

	for (int i = 0; i < MatrixVarStr.size(); i++)
	{
		auto pMatrixVar = pAtmosphereEffect->GetVariableByName(MatrixVarStr[i].c_str())->AsMatrix();
		MatrixVarMap.emplace(MatrixVarStr[i], pMatrixVar);
	}

	for (int i = 0; i < VectorVarStr.size(); i++)
	{
		auto pVectorVar = pAtmosphereEffect->GetVariableByName(VectorVarStr[i].c_str())->AsVector();
		VectorVarMap.emplace(VectorVarStr[i], pVectorVar);
	}

	for (int i = 0; i < ScalarVarStr.size(); i++)
	{
		auto pScalarVar = pAtmosphereEffect->GetVariableByName(ScalarVarStr[i].c_str())->AsScalar();
		ScalarVarMap.emplace(ScalarVarStr[i], pScalarVar);
	}

	for (int i = 0; i < ShaderResourceVarStr.size(); i++)
	{
		auto pShaderResourceVar = pAtmosphereEffect->GetVariableByName(ShaderResourceVarStr[i].c_str())->AsShaderResource();
		ShaderResourceVarMap.emplace(ShaderResourceVarStr[i], pShaderResourceVar);
	}

	D3D11_INPUT_ELEMENT_DESC LayoutDesc[] =
	{
		{"POSITION", 0,	DXGI_FORMAT_R32G32B32_FLOAT,0,0,D3D11_INPUT_PER_VERTEX_DATA,0},
		{"NORMAL",	 0,	DXGI_FORMAT_R32G32B32_FLOAT,0,D3D11_APPEND_ALIGNED_ELEMENT,D3D11_INPUT_PER_VERTEX_DATA,0},
		{"TANGENT",  0, DXGI_FORMAT_R32G32B32_FLOAT,0,D3D11_APPEND_ALIGNED_ELEMENT,D3D11_INPUT_PER_VERTEX_DATA,0},
		{"TEXCOORD", 0,	DXGI_FORMAT_R32G32B32_FLOAT,0,D3D11_APPEND_ALIGNED_ELEMENT,D3D11_INPUT_PER_VERTEX_DATA,0}
	};
	UINT numElements = ARRAYSIZE(LayoutDesc);
	D3DX11_PASS_DESC passDesc;
	AtmosphereTechMap[TechStr[0]]->GetPassByIndex(0)->GetDesc(&passDesc);
	V_RETURN(pd3dDevice->CreateInputLayout(LayoutDesc, numElements, passDesc.pIAInputSignature,
		passDesc.IAInputSignatureSize, &pInputLayout));

	GeometryGenerator geometryGenerator;
	GeometryGenerator::MeshData groundMesh = geometryGenerator.CreateSphere(fInnerRadius, 100, 100);
	groundIndexNum = groundMesh.Indices32.size();
	V_RETURN(D3dBufferDesc::BufferDescVertex(pd3dDevice, &pGroundVB, groundMesh.Vertices.data(), groundMesh.Vertices.size() * sizeof(Vertex)));
	V_RETURN(D3dBufferDesc::BufferDescIndex(pd3dDevice, &pGroundIB, groundMesh.Indices32.data(), groundMesh.Indices32.size() * sizeof(GeometryGenerator::uint32)));

	V_RETURN(D3DX11CreateShaderResourceViewFromFile(pd3dDevice, L"Resource/earthmap1k.dds",
		nullptr, nullptr, &pGroundSRV, nullptr));

	D3D11_TEXTURE2D_DESC PreCompute2DTexDesc;
	ZeroMemory(&PreCompute2DTexDesc, sizeof(PreCompute2DTexDesc));
	PreCompute2DTexDesc.Width = TRANSMITTANCE_TEXTURE_WIDTH;
	PreCompute2DTexDesc.Height = TRANSMITTANCE_TEXTURE_HEIGHT;
	PreCompute2DTexDesc.MipLevels = 1;
	PreCompute2DTexDesc.ArraySize = 1;
	PreCompute2DTexDesc.Format = DXGI_FORMAT_R32G32B32_FLOAT;
	PreCompute2DTexDesc.SampleDesc.Count = 1;
	PreCompute2DTexDesc.SampleDesc.Quality = 0;
	PreCompute2DTexDesc.Usage = D3D11_USAGE_DEFAULT;
	PreCompute2DTexDesc.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;;
	PreCompute2DTexDesc.CPUAccessFlags = 0;
	PreCompute2DTexDesc.MiscFlags = 0;
	V_RETURN(pd3dDevice->CreateTexture2D(&PreCompute2DTexDesc, nullptr, &pTransmittanceTex2D));
	
	D3D11_RENDER_TARGET_VIEW_DESC PreComputeRTVDesc;
	PreComputeRTVDesc.Format = PreCompute2DTexDesc.Format;
	PreComputeRTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
	PreComputeRTVDesc.Texture2D.MipSlice = 0;
	V_RETURN(pd3dDevice->CreateRenderTargetView(pTransmittanceTex2D, &PreComputeRTVDesc, &pTransmittanceRTV));


	return hr;
}

void Atmosphere::ResetInputView(GeometryGenerator::MeshData& viewQuad)
{
	SAFE_RELEASE(pAtmosphereVB);
	SAFE_RELEASE(pAtmosphereIB);
	D3dBufferDesc::BufferDescVertex(pd3dDevice, &pAtmosphereVB, viewQuad.Vertices.data(), viewQuad.Vertices.size() * sizeof(Vertex));
	D3dBufferDesc::BufferDescIndex(pd3dDevice, &pAtmosphereIB, viewQuad.Indices32.data(), viewQuad.Indices32.size() * sizeof(GeometryGenerator::uint32));
}

void Atmosphere::ResetShaderParam()
{
	
}

void Atmosphere::PreCompute()
{

}

void Atmosphere::PreComputeTransmittance()
{

}

void Atmosphere::Render(D3DXMATRIX& viewProj, D3DXVECTOR3& v3CameraPos, GeometryGenerator::MeshData& viewQuad)
{
	ResetInputView(viewQuad);
	fCameraHeight = D3DXVec3Length(&v3CameraPos);
	fCameraHeight2 = fCameraHeight*fCameraHeight;
	pd3dImmediateContext->IASetInputLayout(pInputLayout);
	pd3dImmediateContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	UINT stride = sizeof(Vertex);
	UINT offset = 0;
	D3DXMATRIX world, worldViewProj, worldInvTranspose;
	D3DXMatrixIdentity(&world);
	worldInvTranspose = InverseTranspose(&world);
	D3DXMatrixMultiply(&worldViewProj, &world, &viewProj);

	D3DX11_TECHNIQUE_DESC techDesc;
	ID3DX11EffectTechnique* activeTech;

	activeTech = AtmosphereTechMap["mSkyfromAtmosphere"];
	activeTech->GetDesc(&techDesc);
	for (UINT p = 0; p < techDesc.Passes; ++p)
	{
		ID3DX11EffectPass* pass = activeTech->GetPassByIndex(p);

		pd3dImmediateContext->IASetVertexBuffers(0, 1, &pAtmosphereVB, &stride, &offset);
		pd3dImmediateContext->IASetIndexBuffer(pAtmosphereIB, DXGI_FORMAT_R32_UINT, 0);

		MatrixVarMap["World"]->SetMatrix(world);
		MatrixVarMap["WorldViewProj"]->SetMatrix(worldViewProj);
		MatrixVarMap["WorldInvTranspose"]->SetMatrix(worldInvTranspose);

		VectorVarMap["v3CameraPos"]->SetRawValue(v3CameraPos, 0, sizeof(D3DXVECTOR3));
		VectorVarMap["v3LightPos"]->SetRawValue(SunPos, 0, sizeof(D3DXVECTOR3));
		VectorVarMap["v3InvWavelength"]->SetRawValue(fInvWavelength4, 0, sizeof(D3DXVECTOR3));

		ScalarVarMap["fCameraHeight"]->SetFloat(fCameraHeight);
		ScalarVarMap["fCameraHeight2"]->SetFloat(fCameraHeight2);
		ScalarVarMap["fOuterRadius"]->SetFloat(fOutRadius);
		ScalarVarMap["fOuterRadius2"]->SetFloat(fOutRadius2);
		ScalarVarMap["fInnerRadius"]->SetFloat(fInnerRadius);
		ScalarVarMap["fInnerRadius2"]->SetFloat(fInnerRadius2);
		ScalarVarMap["fKrESun"]->SetFloat(fKrESun);
		ScalarVarMap["fKmESun"]->SetFloat(fKmESun);
		ScalarVarMap["fKr4PI"]->SetFloat(fKr4PI);
		ScalarVarMap["fKm4PI"]->SetFloat(fKm4PI);
		ScalarVarMap["fScale"]->SetFloat(fScale);
		ScalarVarMap["fScaleOverScaleDepth"]->SetFloat(fScaleOverScaleDepth);
		ScalarVarMap["fScaleDepth"]->SetFloat(fScaleDepth);
		ScalarVarMap["fInvScaleDepth"]->SetFloat(fInvScaleDepth);
		ScalarVarMap["fMieG"]->SetFloat(fMieG);
		ScalarVarMap["fMieG2"]->SetFloat(fMieG2);

		pd3dImmediateContext->RSSetState(RenderStates::NoCullRS);

		pass->Apply(0, pd3dImmediateContext);
		pd3dImmediateContext->DrawIndexed(6, 0, 0);
	}

	//activeTech = AtmosphereTechMap["GroundFromAtmosphere"];
	//activeTech->GetDesc(&techDesc);

	//for (UINT p = 0; p < techDesc.Passes; ++p)
	//{
	//	ID3DX11EffectPass* pass = activeTech->GetPassByIndex(p);

	//	pd3dImmediateContext->IASetVertexBuffers(0, 1, &pGroundVB, &stride, &offset);
	//	pd3dImmediateContext->IASetIndexBuffer(pGroundIB, DXGI_FORMAT_R32_UINT, 0);

	//	MatrixVarMap["World"]->SetMatrix(world);
	//	MatrixVarMap["WorldViewProj"]->SetMatrix(worldViewProj);
	//	MatrixVarMap["WorldInvTranspose"]->SetMatrix(worldInvTranspose);

	//	VectorVarMap["v3CameraPos"]->SetRawValue(v3CameraPos, 0, sizeof(D3DXVECTOR3));
	//	VectorVarMap["v3LightPos"]->SetRawValue(SunPos, 0, sizeof(D3DXVECTOR3));
	//	VectorVarMap["v3InvWavelength"]->SetRawValue(fInvWavelength4, 0, sizeof(D3DXVECTOR3));

	//	ScalarVarMap["fCameraHeight"]->SetFloat(fCameraHeight);
	//	ScalarVarMap["fCameraHeight2"]->SetFloat(fCameraHeight2);
	//	ScalarVarMap["fOuterRadius"]->SetFloat(fOutRadius);
	//	ScalarVarMap["fOuterRadius2"]->SetFloat(fOutRadius2);
	//	ScalarVarMap["fInnerRadius"]->SetFloat(fInnerRadius);
	//	ScalarVarMap["fInnerRadius2"]->SetFloat(fInnerRadius2);
	//	ScalarVarMap["fKrESun"]->SetFloat(fKrESun);
	//	ScalarVarMap["fKmESun"]->SetFloat(fKmESun);
	//	ScalarVarMap["fKr4PI"]->SetFloat(fKr4PI);
	//	ScalarVarMap["fKm4PI"]->SetFloat(fKm4PI);
	//	ScalarVarMap["fScale"]->SetFloat(fScale);
	//	ScalarVarMap["fScaleOverScaleDepth"]->SetFloat(fScaleOverScaleDepth);
	//	ScalarVarMap["fScaleDepth"]->SetFloat(fScaleDepth);
	//	ScalarVarMap["fInvScaleDepth"]->SetFloat(fInvScaleDepth);
	//	ScalarVarMap["fMieG"]->SetFloat(fMieG);
	//	ScalarVarMap["fMieG2"]->SetFloat(fMieG2);

	//	ShaderResourceVarMap["GroundMap"]->SetResource(pGroundSRV);

	//	pd3dImmediateContext->RSSetState(RenderStates::CullCounterClockWiseRS);

	//	pass->Apply(0, pd3dImmediateContext);
	//	pd3dImmediateContext->DrawIndexed(groundIndexNum, 0, 0);
	//}
}
