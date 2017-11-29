#include "d3dx11effect.h"
#include "SDKmisc.h"

#define D3D_COMPILE_STANDARD_FILE_INCLUDE ((ID3DInclude*)(UINT_PTR)1)

namespace Vertex
{
	struct Basic32
	{
		Basic32() : Pos(0.0f, 0.0f, 0.0f), Normal(0.0f, 0.0f, 0.0f), Tex(0.0f, 0.0f) {}
		Basic32(const D3DXVECTOR3& p, const D3DXVECTOR3& n, const D3DXVECTOR2& uv)
			: Pos(p), Normal(n), Tex(uv) {}
		Basic32(float px, float py, float pz, float nx, float ny, float nz, float u, float v)
			: Pos(px, py, pz), Normal(nx, ny, nz), Tex(u,v) {}
		D3DXVECTOR3 Pos;
		D3DXVECTOR3 Normal;
		D3DXVECTOR2 Tex;
	};	
}

static HRESULT CompileEffectFromFile(ID3D11Device* pd3dDevice,ID3DX11Effect** pEffect,wchar_t* FileName)
{
	HRESULT hr = S_OK;
	DWORD dwShaderFlags = D3DCOMPILE_ENABLE_STRICTNESS;
#ifdef _DEBUG
    dwShaderFlags |= D3DCOMPILE_DEBUG;
    dwShaderFlags |= D3DCOMPILE_SKIP_OPTIMIZATION;
#endif

	ID3DBlob* pEffectBuffer = nullptr;

    WCHAR szShaderPath[MAX_PATH];
    V_RETURN( DXUTFindDXSDKMediaFileCch( szShaderPath, MAX_PATH, FileName ) );

	ID3DBlob* pErrorBlob = nullptr;
	V_RETURN(D3DX11CompileFromFile(szShaderPath, nullptr, nullptr, nullptr, "fx_5_0", dwShaderFlags, 0, nullptr, &pEffectBuffer, &pErrorBlob, nullptr));
	if (pErrorBlob)
	{
		OutputDebugStringA(reinterpret_cast<const char*>(pErrorBlob->GetBufferPointer()));
		pErrorBlob->Release();
		pErrorBlob = nullptr;
	}
	hr = D3DX11CreateEffectFromMemory(pEffectBuffer->GetBufferPointer(), pEffectBuffer->GetBufferSize(), 0, pd3dDevice, pEffect);
	//hr = D3DX11CompileEffectFromFile( szShaderPath, nullptr, D3D_COMPILE_STANDARD_FILE_INCLUDE, dwShaderFlags, 0, pd3dDevice, pEffect, &pErrorBlob );

    if ( pErrorBlob )
    {
        OutputDebugStringA( reinterpret_cast<const char*>( pErrorBlob->GetBufferPointer() ) );
        pErrorBlob->Release();
    }

    if ( FAILED(hr) )
        return hr;

    SAFE_RELEASE( pEffectBuffer );

	return hr;
}

static D3DXMATRIX InverseTranspose(D3DXMATRIX* m)
{
	D3DXMATRIX out;
	float det = D3DXMatrixDeterminant(m);
	D3DXMatrixInverse(&out,&det,m);
	D3DXMatrixTranspose(&out,&out);
	return out;
}