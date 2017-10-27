struct VertexIn
{
	float3 PosL		: POITION;
	float2 Tex		: TEXCOORD;
};

struct VertexOut
{
	float4 PosH		: SV_POSITION;
	float4 C0		: COLOR0; // The Rayleigh color
	float4 C1		: COLOR1; // The Mie color
	float2 Tex		: TEXCOORD0;
	float3 DirToCam	: POSITION0;
};

static const int nSamples = 2;
static const float fSamples = (float)nSamples;

static const float fScaleDepth = 0.25;
static const float fInvScaleDepth = 1.0 / fScaleDepth;

float scale(float fCos)
{
	float x = 1.0 - fCos;
	return fScaleDepth * exp(-0.00287 + x*(0.459 + x*(3.83 + x*(-6.80 + x*5.25))));
}

float getMiePhase(float fCos, float fCos2, float g, float g2)
{
	return 1.5 * ((1.0 - g2) / (2.0 + g2)) * (1.0 + fCos2) / pow(1.0 + g2 - 2.0*g*fCos, 1.5);
}

float getRayleighPhase(float fCos2)
{
	return 0.75 + 0.75*fCos2;
}

float2 getIntersection(float3 v3Pos,float3 v3Ray,float fDistance2, float fRadius2)
{
	float b = dot(v3Pos,v3Ray);
	float c = fDistance2 - fRadius2;
	float fDet = sqrt(max(0.0,b*b - c));
	return float2(-b-fDet,-b+fDet);
}