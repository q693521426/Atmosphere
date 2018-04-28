//#include "LightHelper.hlsl"
struct DirectionalLight
{
	float4 Ambient;
	float4 Diffuse;
	float4 Specular;
	float3 Direction;
	float pad;
    matrix LightViewProj;
};

struct PointLight
{
	float4 Ambient;
	float4 Diffuse;
	float4 Specular;

	float3 Position;
	float Range;

	float3 Att;
	float pad;
};

struct SpotLight
{
	float4 Ambient;
	float4 Diffuse;
	float4 Specular;

	float3 Position;
	float Range;

	float3 Direction;
	float Spot;

	float3 Att;
	float pad;
};

struct Material
{
	float4 Ambient;
	float4 Diffuse;
	float4 Specular; // w = SpecPower;
	float4 Reflect;
};

void ComputeDirectionalLight(Material mat, DirectionalLight L,
							 float3 normal, float3 toEye,
							 out float4 ambient,
							 out float4 diffuse,
							 out float4 spec)
{
	ambient = float4(0.0f, 0.0f, 0.0f, 0.0f);
	diffuse = float4(0.0f, 0.0f, 0.0f, 0.0f);
	spec	= float4(0.0f, 0.0f, 0.0f, 0.0f);

	float3 lightVec = -L.Direction;

	ambient = mat.Ambient * L.Ambient;

	float diffuseFactor = dot(lightVec, normal);

	// Flatten to avoid dynamic branching
	[flatten]
	if(diffuseFactor > 0.0f)
	{
		float3 v		 = normalize(toEye + lightVec);
		float specFactor = pow(max(dot(v,normal),0.0f),mat.Specular.w);

		diffuse = diffuseFactor * mat.Diffuse * L.Diffuse;
		spec	= specFactor * mat.Specular * L.Specular;
	}
}

void ComputePointLight(Material mat, PointLight L,
					   float3 pos,float3 normal,float3 toEye,
					   out float4 ambient,
					   out float4 diffuse,
					   out float4 spec)
{
	ambient = float4(0.0f, 0.0f, 0.0f, 0.0f);
	diffuse = float4(0.0f, 0.0f, 0.0f, 0.0f);
	spec	= float4(0.0f, 0.0f, 0.0f, 0.0f);

	float3 lightVec = L.Position - pos;
	float d = length(lightVec);

	if(d>L.Range)
		return;

	lightVec /=d;

	ambient = mat.Ambient * L.Ambient;

	float diffuseFactor = dot(lightVec,normal);

	[flatten]
	if(diffuseFactor > 0.0f)
	{
		float3 v = normalize(toEye + lightVec);
		float specFactor = pow(max(dot(v,normal),0.0f),mat.Specular.w);

		diffuse = diffuseFactor * mat.Diffuse * L.Diffuse;
		spec	= specFactor * mat.Specular * L.Specular;
	}

	float att = 1.0f / dot(L.Att,float3(1.0f,d,d*d));

	diffuse *= att;
	spec *= att;
}

void ComputeSpotLight(Material mat, SpotLight L,
					   float3 pos,float3 normal,float3 toEye,
					   out float4 ambient,
					   out float4 diffuse,
					   out float4 spec)
{
	ambient = float4(0.0f,0.0f,0.0f,0.0f);
	diffuse = float4(0.0f,0.0f,0.0f,0.0f);
	spec = float4(0.0f,0.0f,0.0f,0.0f);

	float3 lightVec = L.Position - pos;

	float d = length(lightVec);

	if(d>L.Range)
		return;

	lightVec /= d;

	ambient = mat.Ambient * L.Ambient;
	
	float diffuseFactor = dot(lightVec,normal);

	[flatten]
	if(diffuseFactor>0.f)
	{
		float3 v = normalize(lightVec+toEye);
		float specFactor = pow(max(dot(v,normal),0.0f),mat.Specular.w);

		diffuse = diffuseFactor * mat.Diffuse * L.Diffuse;
		spec = specFactor * mat.Specular * L.Specular;
	}

	float spot = pow(max(dot(L.Direction,-lightVec),0.0f),L.Spot);

	float att = spot/dot(L.Att,float3(1.0f,d,d*d));

	ambient *= spot;
	diffuse *= att;
	spec	*= att;
}

Texture2D<float> g_tex2DShadowMap : register(t0);
Texture2D txDiffuse : register(t1);
SamplerState samLinear : register(s0);

cbuffer cbChangesEveryFrame : register(b0)
{
    matrix World;
    matrix WVP;
};

cbuffer cbMaterial : register(b1)
{
    Material mat;
};

cbuffer LightBuffer : register(b2)
{
    DirectionalLight light;
};

cbuffer FrustumBuffer : register(b3)
{
    float4 ViewPos;
};

//--------------------------------------------------------------------------------------
struct VS_INPUT
{
    float3 Pos : POSITION;
    float3 Nor : NORMAL;
    float2 Tex : TEXCOORD;
};

struct PS_INPUT
{
    float4 Pos : SV_POSITION;
    float3 Nor : NORMAL;
    float2 Tex : TEXCOORD0;
    float4 WorldPos : POSITION;
};

//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------
PS_INPUT VS(VS_INPUT input)
{
    PS_INPUT output = (PS_INPUT) 0;
    float4 pos = float4(input.Pos, 1.f);
    output.WorldPos = mul(pos, World);
    output.Pos = mul(pos, WVP);
    output.Tex = input.Tex;
    output.Nor = normalize(mul(input.Nor, (float3x3) World));
    return output;
}

//float4 CalculateLight(LightPos light, float4 worldPos, float3 normal);

float4 PS(PS_INPUT input) : SV_Target
{
    float3 toEye = ViewPos.xyz - input.WorldPos.xyz;

    float distToEye = length(toEye);

    toEye /= distToEye;

    float4 ambient = float4(0.0f, 0.0f, 0.0f, 0.0f);
    float4 diffuse = float4(0.0f, 0.0f, 0.0f, 0.0f);
    float4 spec = float4(0.0f, 0.0f, 0.0f, 0.0f);

    float4 A, D, S;
    ComputeDirectionalLight(mat, light, input.Nor, toEye,
			A, D, S);

    ambient += A;
    diffuse += D;
    spec += S;
    
    float4 f4PosInLightSpace = mul(input.WorldPos, light.LightViewProj);
    f4PosInLightSpace /= f4PosInLightSpace.w;
    f4PosInLightSpace.xy = f4PosInLightSpace.xy * float2(0.5, -0.5) + float2(0.5, 0.5);
    float shadow = 0;
    float bias = max(0.05, 0.5 * (1 - dot(light.Direction, input.Nor)));
    float2 f2TexelSize = 1.f / 1024.f;
    for (int x = -1; x <= 1;++x)
    {
        for (int y = -1; y <= 1;++y)
        {
            float pcfDepth = g_tex2DShadowMap.SampleLevel(samLinear, f4PosInLightSpace.xy + float2(x, y) * f2TexelSize, 0);
            shadow += f4PosInLightSpace.z - bias > pcfDepth ? 1 : 0;
        }
    }
    shadow /= 9;
    if (f4PosInLightSpace.z>1)
        shadow = 0;

    float4 texColor = float4(1.f, 1.f, 1.f, 1.f);
    float4 color = (ambient + (1.0 - shadow) * (diffuse + spec)) * texColor;
    return saturate(color);
}