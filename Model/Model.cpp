#include "DXUT.h"
#include "Model.h"

Model::Model():
	m_pVertexShader(nullptr),
	m_pPixelShader(nullptr),
	m_pVertexBuffer(nullptr),
	m_pIndexBuffer(nullptr),
	m_pCBChangesEveryFrame(nullptr),
	m_pLightBuffer (nullptr),
	m_pFrustumBuffer (nullptr),
	m_pTextureRV(nullptr),
	m_pSamplerLinear(nullptr),
	ModelHeight(1800.f),
	ModelScaling(100.f)
{
}


Model::~Model()
{
}

bool Model::Initialize()
{
#if LOAD_MODEL
	m_files[0] = new char[256];
	strcpy_s(m_files[0],256,"Resource/sulan2.fbx");
#endif

#if LIGHT_SPHERE
	m_LightSphere = new Sphere();
	m_LightSphere->Initialize(D3DXVECTOR3(m_Light.LightPos.x, m_Light.LightPos.y, m_Light.LightPos.z), 0.2);
#endif
	return true;
}

void Model::Release()
{
#if LIGHT_SPHERE
	SAFE_RELEASE(m_LightSphere);
#endif

#if LOAD_MODEL
	SAFE_RELEASE(m_pVertexShader);
	SAFE_RELEASE(m_pPixelShader);
	SAFE_RELEASE(m_pVertexBuffer);
	SAFE_RELEASE(m_pIndexBuffer);
	SAFE_RELEASE(m_pTextureRV);
	SAFE_RELEASE(m_pSamplerLinear);
	SAFE_RELEASE(m_pCBChangesEveryFrame);
	SAFE_RELEASE(m_pLightBuffer);
	SAFE_RELEASE(m_pFrustumBuffer);
	for (int i = 0; i < NUMBER_OF_MODELS; i++)
	{
		if(m_pFbxDX11[i])
		{
			m_pFbxDX11[i]->Release();
			delete m_pFbxDX11[i];
			m_pFbxDX11[i] = nullptr;
		}

	}
#endif
}

HRESULT Model::OnD3D11CreateDevice(ID3D11Device* pd3dDevice, ID3D11DeviceContext* pd3dImmediateContext)
{
#if LIGHT_SPHERE
	m_LightSphere->OnD3D11CreateDevice(pd3dDevice, pd3dImmediateContext);
#endif

#if LOAD_MODEL
	HRESULT hr = S_OK;

	for (DWORD i = 0; i < NUMBER_OF_MODELS; i++)
	{
		m_pFbxDX11[i] = new FBX_LOADER::CFBXRenderDX11;
		hr = m_pFbxDX11[i]->LoadFBX(m_files[i], pd3dDevice, pd3dImmediateContext);
	}
	if (FAILED(hr))
	{
		MessageBox(NULL,
			L"FBX Error", L"Error", MB_OK);
		return hr;
	}

	for (int i = 0; i < NUMBER_OF_MODELS; i++)
	{
		if (m_files[i])
			delete m_files[i];
		m_files[i] = nullptr;
	}
	ID3DBlob* pVSBlob = nullptr;
	V_RETURN(CompileShader(L"Model/LoadModel.hlsl", "VS", "vs_5_0", &pVSBlob));

	hr = pd3dDevice->CreateVertexShader(pVSBlob->GetBufferPointer(), pVSBlob->GetBufferSize(), nullptr, &m_pVertexShader);
	if (FAILED(hr))
	{
		SAFE_RELEASE(pVSBlob);
		return hr;
	}

	D3D11_INPUT_ELEMENT_DESC layout[] =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, D3D11_APPEND_ALIGNED_ELEMENT, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, D3D11_APPEND_ALIGNED_ELEMENT, D3D11_INPUT_PER_VERTEX_DATA, 0 },
	};

	UINT numElements = ARRAYSIZE(layout);
	for (DWORD i = 0; i < NUMBER_OF_MODELS; ++i)
	{
		hr = m_pFbxDX11[i]->CreateInputLayout(pd3dDevice, pVSBlob->GetBufferPointer(), pVSBlob->GetBufferSize(), layout, numElements);
	}
	SAFE_RELEASE(pVSBlob);
	if (FAILED(hr))
		return hr;

	ID3DBlob* pPSBlob = nullptr;
	V_RETURN(CompileShader(L"Model/LoadModel.hlsl","PS", "ps_5_0",&pPSBlob));

	// Create the pixel shader
	hr = pd3dDevice->CreatePixelShader(pPSBlob->GetBufferPointer(), pPSBlob->GetBufferSize(), nullptr, &m_pPixelShader);
	SAFE_RELEASE(pPSBlob);
	if (FAILED(hr))
		return hr;

	// Create the constant buffers
	D3D11_BUFFER_DESC bd;
	ZeroMemory(&bd, sizeof(bd));
	bd.Usage = D3D11_USAGE_DYNAMIC;
	bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	bd.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
	bd.ByteWidth = sizeof(CBChangesEveryFrame);
	V_RETURN(pd3dDevice->CreateBuffer(&bd, nullptr, &m_pCBChangesEveryFrame));

	bd.ByteWidth = sizeof(DirectionalLight);
	V_RETURN(pd3dDevice->CreateBuffer(&bd, nullptr, &m_pLightBuffer));

	bd.ByteWidth = sizeof(FrustumBuffer);
	V_RETURN(pd3dDevice->CreateBuffer(&bd, nullptr, &m_pFrustumBuffer));

	// Create the sample state
	D3D11_SAMPLER_DESC sampDesc;
	ZeroMemory(&sampDesc, sizeof(sampDesc));
	sampDesc.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
	sampDesc.AddressU = D3D11_TEXTURE_ADDRESS_WRAP;
	sampDesc.AddressV = D3D11_TEXTURE_ADDRESS_WRAP;
	sampDesc.AddressW = D3D11_TEXTURE_ADDRESS_WRAP;
	sampDesc.ComparisonFunc = D3D11_COMPARISON_NEVER;
	sampDesc.MinLOD = 0;
	sampDesc.MaxLOD = D3D11_FLOAT32_MAX;
	V_RETURN(pd3dDevice->CreateSamplerState(&sampDesc, &m_pSamplerLinear));
#endif
	return hr;
}

void Model::RenderModel(ID3D11Device* pd3dDevice, ID3D11DeviceContext* pd3dImmediateContext,ID3D11ShaderResourceView* pShadowMap,const D3DXMATRIX& viewProj,bool IsOnlyZPass)
{
#if LOAD_MODEL
	for (DWORD i = 0; i < NUMBER_OF_MODELS; i++)
	{
		size_t nodeCount = m_pFbxDX11[i]->GetNodeCount();

		if(IsOnlyZPass)
		{
			pd3dImmediateContext->VSSetShader(m_pVertexShader, nullptr, 0);
			pd3dImmediateContext->PSSetShader(m_pPixelShader, nullptr, 0);
		}
		else
		{
			pd3dImmediateContext->VSSetShader(m_pVertexShader, nullptr, 0);
			pd3dImmediateContext->PSSetShader(m_pPixelShader, nullptr, 0);
		}

		HRESULT hr;
		D3D11_MAPPED_SUBRESOURCE MappedResource;

		V(pd3dImmediateContext->Map(m_pLightBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &MappedResource));
		auto  LB = reinterpret_cast<DirectionalLight*>(MappedResource.pData);
		*LB = m_DirectionalLight;
		D3DXMatrixTranspose(&LB->mLightViewProj, &m_LightViewProj);
		pd3dImmediateContext->Unmap(m_pLightBuffer, 0);

		V(pd3dImmediateContext->Map(m_pFrustumBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &MappedResource));
		auto  FB = reinterpret_cast<FrustumBuffer*>(MappedResource.pData);
		FB->ViewPos = D3DXVECTOR4(m_ViewPos, 1.0f);
		pd3dImmediateContext->Unmap(m_pFrustumBuffer, 0);

		pd3dImmediateContext->PSSetConstantBuffers(2, 1, &m_pLightBuffer);
		pd3dImmediateContext->PSSetConstantBuffers(3, 1, &m_pFrustumBuffer);

		pd3dImmediateContext->PSSetShaderResources(0, 1, &pShadowMap);

		for (size_t j = 0; j < nodeCount; j++)
		{
			D3DXMATRIX m_Local;
			m_pFbxDX11[i]->GetNodeMatrix(j, m_Local);

			V(pd3dImmediateContext->Map(m_pCBChangesEveryFrame, 0, D3D11_MAP_WRITE_DISCARD, 0, &MappedResource));
			auto pCB = reinterpret_cast<CBChangesEveryFrame*>(MappedResource.pData);

			D3DXMATRIX change(1.f, 0.f, 0.f, 0.f,
				0.f, 0.f, 1.f, 0.f,
				0.f, 1.f, 0.f, 0.f,
				0.f, 0.f, 0.f, 1.f);
			D3DXMATRIX translation;
			D3DXMatrixTranslation(&translation, 40 * ModelScaling, ModelHeight - 1 * ModelScaling, 565 * ModelScaling);
			D3DXMATRIX scale;
			D3DXMatrixScaling(&scale, ModelScaling, ModelScaling, ModelScaling);

			m_Local = m_Local * scale;
			m_Local = m_Local * change;
			m_Local = m_Local * translation;

			auto WVP = m_Local * viewProj;
			D3DXMatrixTranspose(&pCB->mWorld, &m_Local);
			D3DXMatrixTranspose(&pCB->mWVP, &WVP);
			pd3dImmediateContext->Unmap(m_pCBChangesEveryFrame, 0);

			pd3dImmediateContext->VSSetConstantBuffers(0, 1, &m_pCBChangesEveryFrame);

			//SetModelInstancingMatrix(pd3dImmediateContext);

			FBX_LOADER::MATERIAL_DATA material = m_pFbxDX11[i]->GetNodeMaterial(j);

			if (material.pMaterialCb)
				pd3dImmediateContext->UpdateSubresource(material.pMaterialCb, 0, nullptr, &material.materialConstantData, 0, 0);

			pd3dImmediateContext->PSSetShaderResources(1, 1, &material.pSRV);
			pd3dImmediateContext->PSSetConstantBuffers(1, 1, &material.pMaterialCb);
			pd3dImmediateContext->PSSetSamplers(0, 1, &material.pSampler);

			m_pFbxDX11[i]->RenderNode(pd3dImmediateContext, j);
		}
	}
#endif
}

void Model::Render(ID3D11Device* pd3dDevice, ID3D11DeviceContext* pd3dImmediateContext, ID3D11ShaderResourceView* pShadowMap)
{
	RenderModel(pd3dDevice, pd3dImmediateContext, pShadowMap,m_ViewProj, false);
#if LIGHT_SPHERE
	m_LightSphere->SetWVP(m_WVP);
	m_LightSphere->Render(pd3dDevice, pd3dImmediateContext);
#endif
}

void Model::RenderShadowMap(ID3D11Device* pd3dDevice, ID3D11DeviceContext* pd3dImmediateContext, const D3DXVECTOR3& SunDir,UINT shadowMapDim)
{
	D3DXVECTOR3 lightDir = -SunDir;
	D3DXVECTOR3 lightSpaceZ = lightDir;
	D3DXVECTOR3 lightSpaceX = D3DXVECTOR3(1, 0, 0);
	D3DXVECTOR3 lightSpaceY;
	D3DXVec3Cross(&lightSpaceY, &lightSpaceX, &lightSpaceZ);
	D3DXVec3Cross(&lightSpaceX, &lightSpaceZ, &lightSpaceY);

	D3DXVec3Normalize(&lightSpaceX, &lightSpaceX);
	D3DXVec3Normalize(&lightSpaceY, &lightSpaceY);
	D3DXVec3Normalize(&lightSpaceZ, &lightSpaceZ);

	D3DXMatrixIdentity(&m_LightView);
	m_LightView._11 = lightSpaceX.x;
	m_LightView._21 = lightSpaceX.y;
	m_LightView._31 = lightSpaceX.z;

	m_LightView._12 = lightSpaceY.x;
	m_LightView._22 = lightSpaceY.y;
	m_LightView._32 = lightSpaceY.z;

	m_LightView._13 = lightSpaceZ.x;
	m_LightView._23 = lightSpaceZ.y;
	m_LightView._33 = lightSpaceZ.z;

	D3DXMATRIX camInvViewProj;
	float det = D3DXMatrixDeterminant(&m_ViewProj);
	D3DXMatrixInverse(&camInvViewProj, &det, &m_ViewProj);

	D3DXVECTOR3 f3ViewPosInLightSpace;
	D3DXVec3TransformCoord(&f3ViewPosInLightSpace, &m_ViewPos, &m_LightView);

	D3DXMATRIX camProjSpaceToLightSpace = camInvViewProj * m_LightView;

	//D3DXVECTOR3 f3MinXYZ(FLT_MAX, FLT_MAX, FLT_MAX);
	//D3DXVECTOR3 f3MaxXYZ(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	D3DXVECTOR3 f3MinXYZ(f3ViewPosInLightSpace);
	D3DXVECTOR3 f3MaxXYZ(f3ViewPosInLightSpace);

	for (int iClipPlaneCorner = 0; iClipPlaneCorner < 8; ++iClipPlaneCorner)
	{
		D3DXVECTOR3 f3PlaneCornerProjSpace((iClipPlaneCorner & 0x01) ? +1.f : -1.f,
										(iClipPlaneCorner & 0x02) ? +1.f : -1.f,
										(iClipPlaneCorner & 0x04) ? 1.f : 0.f);
		D3DXVECTOR3 f3PlaneCornerLightSpace;
		D3DXVec3TransformCoord(&f3PlaneCornerLightSpace, &f3PlaneCornerProjSpace, &camProjSpaceToLightSpace);
		D3DXVec3Minimize(&f3MinXYZ, &f3MinXYZ, &f3PlaneCornerLightSpace);
		D3DXVec3Maximize(&f3MaxXYZ, &f3MaxXYZ, &f3PlaneCornerLightSpace);
	}
	float fEarthRadius = 6360.f;
	f3MinXYZ.z -= fEarthRadius * sqrt(2.f);

	//float fShadowMapDim = (float)shadowMapDim;
	//float fXExt = (f3MaxXYZ.x - f3MinXYZ.x) * (1 + 1.f / fShadowMapDim);
	//float fYExt = (f3MaxXYZ.y - f3MinXYZ.y) * (1 + 1.f / fShadowMapDim);
	//const float fExtStep = 2.f;
	//fXExt = pow(fExtStep, ceil(log(fXExt) / log(fExtStep)));
	//fYExt = pow(fExtStep, ceil(log(fYExt) / log(fExtStep)));
	//// Align cascade center with the shadow map texels to alleviate temporal aliasing
	//float fXCenter = (f3MaxXYZ.x + f3MinXYZ.x) / 2.f;
	//float fYCenter = (f3MaxXYZ.y + f3MinXYZ.y) / 2.f;
	//float fTexelXSize = fXExt / fShadowMapDim;
	//float fTexelYSize = fXExt / fShadowMapDim;
	//fXCenter = floor(fXCenter / fTexelXSize) * fTexelXSize;
	//fYCenter = floor(fYCenter / fTexelYSize) * fTexelYSize;
	//// Compute new cascade min/max xy coords
	//f3MaxXYZ.x = fXCenter + fXExt / 2.f;
	//f3MinXYZ.x = fXCenter - fXExt / 2.f;
	//f3MaxXYZ.y = fYCenter + fYExt / 2.f;
	//f3MinXYZ.y = fYCenter - fYExt / 2.f;

	D3DXVECTOR3 f3LightSpaceScale, f3LightSpaceScaledBias;

	f3LightSpaceScale.x = 2.f / (f3MaxXYZ.x - f3MinXYZ.x);
	f3LightSpaceScale.y = 2.f / (f3MaxXYZ.y - f3MinXYZ.y);
	f3LightSpaceScale.z = 1.f / (f3MaxXYZ.z - f3MinXYZ.z);
	// Apply bias to shift the extent to [-1,1]x[-1,1]x[0,1]
	f3LightSpaceScaledBias.x = -f3MinXYZ.x * f3LightSpaceScale.x - 1.f;
	f3LightSpaceScaledBias.y = -f3MinXYZ.y * f3LightSpaceScale.y - 1.f;
	f3LightSpaceScaledBias.z = -f3MinXYZ.z * f3LightSpaceScale.z + 0.f;
	D3DXMATRIX ScaleMatrix;
	D3DXMatrixScaling(&ScaleMatrix, f3LightSpaceScale.x, f3LightSpaceScale.y, f3LightSpaceScale.z);
	D3DXMATRIX ScaledBiasMatrix;
	D3DXMatrixTranslation(&ScaledBiasMatrix, f3LightSpaceScaledBias.x, f3LightSpaceScaledBias.y, f3LightSpaceScaledBias.z);

	// Note: bias is applied after scaling!
	m_LightProj = ScaleMatrix * ScaledBiasMatrix;
	//D3DXMatrixOrthoOffCenterLH(&lightProj, f3MinXYZ.x, f3MaxXYZ.x, f3MinXYZ.y, f3MaxXYZ.y, f3MinXYZ.z, f3MaxXYZ.z);
	m_LightViewProj = m_LightView * m_LightProj;
	RenderModel(pd3dDevice, pd3dImmediateContext, nullptr,m_LightViewProj, true);
}

void Model::Resize(const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc)
{
#if LIGHT_SPHERE
	m_LightSphere->SetWVP(m_WVP);
#endif
}

void Model::SetModelHeight(float h)
{
	ModelHeight = h;
}

void Model::SetModelScaling(float s)
{
	ModelScaling = s;
}

void Model::SetViewProj(const D3DXMATRIX& viewProj)
{
	this->m_ViewProj = viewProj;
}


void Model::SetViewPos(const D3DXVECTOR3& viewPos)
{
	m_ViewPos = viewPos;
}

void Model::SetLight(const DirectionalLight* l)
{
	m_DirectionalLight = *l;
}

D3DXMATRIX Model::GetLightView()const
{
	return m_LightView;
}

D3DXMATRIX Model::GetLightProj()const
{
	return m_LightProj;
}

D3DXMATRIX Model::GetLightViewProj()const
{
	return m_LightViewProj;
}