#pragma once

#define LOAD_MODEL		1
#define LIGHT_SPHERE	0

#include "FBXLoader/CFBXRendererDX11.h"
#include "ModelConstants.h"

#if LIGHT_SPHERE
#include "Sphere.h"
#endif

using namespace DirectX;

class Model
{

public:
	Model();
	~Model();
	bool Initialize();
	HRESULT OnD3D11CreateDevice(ID3D11Device*,ID3D11DeviceContext*);
	void Release();

	void Render(ID3D11Device*, ID3D11DeviceContext*, ID3D11ShaderResourceView*);
	void RenderShadowMap(ID3D11Device*, ID3D11DeviceContext*, const D3DXVECTOR3&, UINT);
	void Model::RenderModel(ID3D11Device*, ID3D11DeviceContext*,ID3D11ShaderResourceView*, const D3DXMATRIX&, bool);

	void Resize(const DXGI_SURFACE_DESC*);

	void SetViewProj(const D3DXMATRIX&);
	void SetViewPos(const D3DXVECTOR3&);
	void SetLight(const DirectionalLight*);

	D3DXMATRIX GetLightView()const;
	D3DXMATRIX GetLightProj()const;
	D3DXMATRIX GetLightViewProj()const;

	void SetModelHeight(float h);
	void SetModelScaling(float s);
private:
	float							ModelHeight;
	float							ModelScaling;
	D3DXMATRIX						m_World;
	D3DXMATRIX						m_ViewProj;
	D3DXMATRIX						m_LightView;
	D3DXMATRIX						m_LightProj;
	D3DXMATRIX						m_LightViewProj;
	D3DXVECTOR3						m_ViewPos;
	const static DWORD				NUMBER_OF_MODELS = 1;
	ID3D11VertexShader*				m_pVertexShader;
	ID3D11PixelShader*				m_pPixelShader;
	ID3D11Buffer*					m_pVertexBuffer;
	ID3D11Buffer*					m_pIndexBuffer;
	ID3D11Buffer*					m_pCBChangesEveryFrame;
	ID3D11Buffer*					m_pLightBuffer;
	ID3D11Buffer*					m_pFrustumBuffer;
	ID3D11ShaderResourceView*		m_pTextureRV;
	ID3D11SamplerState*				m_pSamplerLinear;

	D3DXVECTOR3						m_f3MinXYZ;
	D3DXVECTOR3						m_f3MaxXYZ;
#if LIGHT_SPHERE
	Sphere*							m_LightSphere;
#endif
	FBX_LOADER::CFBXRenderDX11*		m_pFbxDX11[NUMBER_OF_MODELS];
	char*							m_files[NUMBER_OF_MODELS];
	DirectionalLight				m_DirectionalLight;

};

