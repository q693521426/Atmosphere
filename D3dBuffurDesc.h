#include <DXUT.h>

class D3dBufferDesc
{
public:
	static HRESULT BufferDesc(ID3D11Device* pd3dDevice,ID3D11Buffer** ppbuf,
								D3D11_USAGE Usage,UINT BindFlags,
								UINT CPUAccessFlags,UINT MiscFlags,
								UINT ByteWidth,const void* pSysMem = nullptr)
	{
		HRESULT hr = S_OK;
		D3D11_BUFFER_DESC bd;

		ZeroMemory(&bd, sizeof(bd));
		bd.Usage = Usage;
		bd.BindFlags = BindFlags;
		bd.ByteWidth = ByteWidth;
		bd.CPUAccessFlags = CPUAccessFlags;
		bd.MiscFlags = MiscFlags;
		if(pSysMem)
		{
			D3D11_SUBRESOURCE_DATA subsource;
			subsource.pSysMem = pSysMem;
			V_RETURN(pd3dDevice->CreateBuffer(&bd, &subsource, ppbuf));
		}
		else
			V_RETURN(pd3dDevice->CreateBuffer(&bd, nullptr, ppbuf));
		return hr;
	}

	static HRESULT BufferDescDynamicConst(ID3D11Device* pd3dDevice,ID3D11Buffer** ppbuf,UINT ByteWidth)
	{
		return BufferDesc(pd3dDevice,ppbuf,
			D3D11_USAGE_DYNAMIC,
			D3D11_BIND_CONSTANT_BUFFER,
			D3D11_CPU_ACCESS_WRITE,
			0,ByteWidth);
	}

	static HRESULT BufferDescVertex(ID3D11Device* pd3dDevice,ID3D11Buffer** ppbuf,const void* pSysMem,UINT ByteWidth)
	{
		return BufferDesc(pd3dDevice,ppbuf,
			D3D11_USAGE_IMMUTABLE,
			D3D11_BIND_VERTEX_BUFFER,
			0,0,ByteWidth,pSysMem);
	}

	static HRESULT BufferDescIndex(ID3D11Device* pd3dDevice,ID3D11Buffer** ppbuf,const void* pSysMem,UINT ByteWidth)
	{
		return BufferDesc(pd3dDevice,ppbuf,
			D3D11_USAGE_IMMUTABLE,
			D3D11_BIND_INDEX_BUFFER,
			0,0,ByteWidth,pSysMem);
	}
};