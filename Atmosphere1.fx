#include "Common.fx"

struct VertexIn
{
	float3 PosL		: POSITION;
	float3 Normal	: NORMAL;
	float3 Tangent	: TANGENT;
	float2 Tex		: TEXCOORD;
};

struct VertexOut
{
	float4 PosH		: SV_POSITION;
	float4 C0		: COLOR0; // The Rayleigh color
	float4 C1		: COLOR1; // The Mie color
	float2 Tex		: TEXCOORD0;
	float3 DirToCam	: POSITION0;
	float3 PosW		: POSITION1;
};

cbuffer cbPerFrame
{
	float3 v3CameraPos;
	float3 v3LightPos;
	float3 v3InvWavelength;
	float fCameraHeight;
	float fCameraHeight2;
	float fOuterRadius;
	float fOuterRadius2;
	float fInnerRadius;
	float fInnerRadius2;
	float fKrESun;
	float fKmESun;
	float fKr4PI;
	float fKm4PI;
	float fScale;					// 1 / (fOuterRadius - fInnerRadius)
	float fScaleOverScaleDepth;		// fScale / fScaleDepth
	float fScaleDepth;
	float fInvScaleDepth;
	float fMieG;
	float fMieG2;
};

cbuffer cbPerObject
{
	float4x4 World;
	float4x4 WorldInvTranspose;
	float4x4 WorldViewProj;
};

Texture2D GroundMap;

SamplerState samAnisotropic
{
	Filter = ANISOTROPIC;
	MaxAnisotropy = 4;

	AddressU = WRAP;
	AddressV = WRAP;
};


VertexOut GroundFromAtmosphereVS(VertexIn vin)
{
	float3 v3Pos = vin.PosL;
	float3 v3Ray = v3Pos - v3CameraPos;
	v3Pos = normalize(v3Pos);
	float fFar = length(v3Ray);
	v3Ray /= fFar;

	float3 v3Start = v3CameraPos;
	//float fDepth = exp((fInnerRadius-fCameraHeight)*fInvScaleDepth); 
	float fDepth = exp((fInnerRadius-fCameraHeight)*fScaleOverScaleDepth); 
	float fCameraAngle = dot(-v3Ray,v3Pos);
	float fLightAngle = dot(v3Pos,v3LightPos);
	float fCameraScale = scale(fCameraAngle, fScaleDepth);
	float fLightScale = scale(fLightAngle, fScaleDepth);
	float fCameraOffset = fDepth * fCameraScale;
	float fTemp = (fLightScale + fCameraScale);

	float fSampleLength = fFar / fSamples;
	float fScaledLength = fSampleLength * fScale;
	float3 v3SampleRay = v3Ray * fSampleLength;
	float3 v3SamplePoint = v3Start + v3SampleRay * 0.5;

	float3 v3FrontColor = float3(0.f,0.f,0.f);
	float3 v3Attenuate;
	for(int i=0;i<nSamples;i++)
	{
		float fHeight = length(v3SamplePoint);
		float fDepth = exp(fScaleOverScaleDepth * (fInnerRadius - fHeight));
		float fScatter = fDepth * fTemp - fCameraOffset;
		v3Attenuate = exp(-fScatter * (v3InvWavelength * fKr4PI + fKm4PI));
		v3FrontColor += v3Attenuate * (fDepth * fScaledLength); 
		v3SamplePoint += v3SampleRay;
	}

	VertexOut vout;
	vout.PosH = mul(float4(vin.PosL,1.0f),WorldViewProj);
	vout.C0.xyz = v3FrontColor * (v3InvWavelength * fKrESun + fKmESun);
	vout.C1.xyz = v3Attenuate;
	vout.C0.w = 1;
	vout.C1.w = 1;
	vout.Tex = vin.Tex;
	vout.DirToCam = v3CameraPos - v3Pos;

	return vout;
}

float4 GroundFromAtmospherePS(VertexOut pin) :SV_Target
{
	float4 texColor = GroundMap.Sample(samAnisotropic,pin.Tex);
	return pin.C0 + pin.C1*texColor;
}

VertexOut SkyFromAtmosphereVS(VertexIn vin)
{
	float3 v3Pos = vin.PosL;
	float3 v3Ray = v3Pos - v3CameraPos;
	float fFar = length(v3Ray);
	v3Ray /= fFar;

	float3 v3Start = v3CameraPos;
	float fHeight = length(v3Start);
	float fDepth = exp((fInnerRadius-fCameraHeight)*fScaleOverScaleDepth);
	float fStartAngle = dot(v3Ray,v3Start) / fHeight;
	float fStartOffset = fDepth*scale(fStartAngle, fScaleDepth);

	float fSampleLength = fFar / fSamples;
	float fScaledLength = fSampleLength * fScale;
	float3 v3SampleRay = v3Ray * fSampleLength;
	float3 v3SamplePoint = v3Start + v3SampleRay * 0.5f;

	float3 v3FrontColor = float3(0.f,0.f,0.f);
	for(int i=0;i<nSamples;i++)
	{
		float fHeight = length(v3SamplePoint);
		float fDepth = exp((fInnerRadius-fCameraHeight)*fScaleOverScaleDepth);
		float fLightAngle = dot(v3LightPos,v3SamplePoint)/fHeight;
		float fCameraAngle = dot(v3Ray,v3SamplePoint)/fHeight;
		float fScatter = (fStartOffset + fDepth*(scale(fLightAngle, fScaleDepth)-scale(fCameraAngle, fScaleDepth)));
		float3 v3Attenuate = exp(-fScatter * (v3InvWavelength * fKr4PI + fKm4PI));
		v3FrontColor += v3Attenuate * (fDepth * fScaledLength);
		v3SamplePoint += v3SampleRay;
	}

	VertexOut vout;
	vout.PosH = mul(float4(vin.PosL,1.0f),WorldViewProj);
	vout.C0.xyz = v3FrontColor * (v3InvWavelength * fKrESun);
	vout.C1.xyz = v3FrontColor * fKmESun;
	vout.C0.w = 1;
	vout.C1.w = 1;
//  vout.Tex = v3CameraPos - v3Pos;	
	vout.Tex = vin.Tex;
	vout.DirToCam = v3CameraPos - v3Pos;
	return vout;
}

float4 SkyFromAtmospherePS(VertexOut pin) :SV_Target
{
	float3 v3Direction = pin.DirToCam;
	float fCos = dot(v3LightPos,v3Direction)/length(v3Direction);
	float fCos2 = fCos*fCos;
	float4 color = getRayleighPhase(fCos2) * pin.C0 + 
		getMiePhase(fCos,fCos2,fMieG,fMieG2) * pin.C1;
	color.w = color.z;
	return color;
}

VertexOut TestVS(VertexIn vin)
{
	VertexOut vout;
	vout.PosH = mul(float4(vin.PosL,1.0f),WorldViewProj);
	vout.Tex = vin.Tex;
	vout.C0 = float4(0.f, 0.f, 0.f, 0.f);
	vout.C1 = float4(0.f, 0.f, 0.f, 0.f);
	vout.DirToCam = v3CameraPos - vin.PosL;
	return vout;
}

float4 TestPS(VertexOut pin) :SV_Target
{
	//return float4(0,0,0,1);
	float4 texColor = GroundMap.Sample(samAnisotropic,pin.Tex);
	return texColor;
}

VertexOut mSkyfromAtmosphereVS(VertexIn vin)
{
	VertexOut vout;
	vout.PosH = mul(float4(vin.PosL, 1.0f), WorldViewProj);
	vout.Tex = vin.Tex;
	vout.C0 = float4(0.f, 0.f, 0.f, 0.f);
	vout.C1 = float4(0.f, 0.f, 0.f, 0.f);
	vout.DirToCam = v3CameraPos - vin.PosL;
	vout.PosW = vin.PosL;
	return vout;
}

float4 mSkyfromAtmospherePS(VertexOut pin) :SV_Target
{
	float3 v3Pos = pin.PosW;
	float3 v3Ray = v3Pos - v3CameraPos;
	float fFar = length(v3Ray);
	v3Ray /= fFar;

	float3 v3Start = v3CameraPos;
	float fHeight = length(v3Start);
	float fDepth = exp((fInnerRadius - fCameraHeight)*fScaleOverScaleDepth);
	float fStartAngle = dot(v3Ray,v3Start) / fHeight;
	float fStartOffset = fDepth*scale(fStartAngle, fScaleDepth);

	float fSampleLength = fFar / fSamples;
	float fScaledLength = fSampleLength * fScale;
	float3 v3SampleRay = v3Ray * fSampleLength;
	float3 v3SamplePoint = v3Start + v3SampleRay * 0.5f;

	float3 v3FrontColor = float3(0.f,0.f,0.f);
	for (int i = 0; i<nSamples; i++)
	{
		float fHeight = length(v3SamplePoint);
		float fDepth = exp((fInnerRadius - fCameraHeight)*fScaleOverScaleDepth);
		float fLightAngle = dot(v3LightPos,v3SamplePoint) / fHeight;
		float fCameraAngle = dot(v3Ray,v3SamplePoint) / fHeight;
		float fScatter = (fStartOffset + fDepth*(scale(fLightAngle, fScaleDepth) - scale(fCameraAngle, fScaleDepth)));
		float3 v3Attenuate = exp(-fScatter * (v3InvWavelength * fKr4PI + fKm4PI));
		v3FrontColor += v3Attenuate * (fDepth * fScaledLength);
		v3SamplePoint += v3SampleRay;
	}

	pin.C0.xyz = v3FrontColor * (v3InvWavelength * fKrESun);
	pin.C1.xyz = v3FrontColor * fKmESun;
	pin.C0.w = 1;
	pin.C1.w = 1;

	float3 v3Direction = pin.DirToCam;
	float fCos = dot(v3LightPos, v3Direction) / length(v3Direction);
	float fCos2 = fCos*fCos;
	float4 color = getRayleighPhase(fCos2) * pin.C0 +
		getMiePhase(fCos, fCos2, fMieG, fMieG2) * pin.C1;
	color.w = color.z;
	return color;
}

technique11 GroundFromAtmosphere
{
	pass P0
    {
        SetVertexShader( CompileShader( vs_5_0, GroundFromAtmosphereVS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, GroundFromAtmospherePS() ) );
    }
}

technique11 GroundFromSpace
{
	pass P0
    {
        SetVertexShader( CompileShader( vs_5_0, GroundFromAtmosphereVS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, GroundFromAtmospherePS() ) );
    }
}

technique11 SkyFromAtmosphere
{
	pass P0
    {
        SetVertexShader( CompileShader( vs_5_0, SkyFromAtmosphereVS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, SkyFromAtmospherePS() ) );
    }
}

technique11 SkyFromSpace
{
	pass P0
    {
        SetVertexShader( CompileShader( vs_5_0, GroundFromAtmosphereVS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, GroundFromAtmospherePS() ) );
    }
}

technique11 SpaceFromAtmosphere
{
	pass P0
    {
        SetVertexShader( CompileShader( vs_5_0, GroundFromAtmosphereVS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, GroundFromAtmospherePS() ) );
    }
}

technique11 SpaceFromSpace
{
	pass P0
    {
        SetVertexShader( CompileShader( vs_5_0, GroundFromAtmosphereVS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, GroundFromAtmospherePS() ) );
    }
}

technique11 Test
{
	pass P0
    {
        SetVertexShader( CompileShader( vs_5_0, TestVS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, TestPS() ) );
    }
}

technique11 mSkyfromAtmosphere
{
	pass P0
	{
		SetVertexShader(CompileShader(vs_5_0, mSkyfromAtmosphereVS()));
		SetGeometryShader(NULL);
		SetPixelShader(CompileShader(ps_5_0, mSkyfromAtmospherePS()));
	}
};
