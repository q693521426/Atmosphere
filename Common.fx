static const float PI = 3.141592654f;

#define USE_LUT_PARAMETERIZATION 1
#define USE_TRANSMITTANCE_ANALYTIC 0
#define USE_OZONE_DENSITY 1
#define USE_SCATTER_COMBINED    1
#define USE_OPTICAL_LUT    0
#define USE_DEPTH_WEIGHT 0
#define USE_SHADOW_OBJECT_TO_EARTH 0
#define FLT_MAX 3.402823466e+38f


SamplerState samLinearClamp
{
    Filter = MIN_MAG_MIP_LINEAR;
    AddressU = Clamp;
    AddressV = Clamp;
};

SamplerState samPointClamp
{
    Filter = MIN_MAG_MIP_POINT;
    AddressU = Clamp;
    AddressV = Clamp;
};


SamplerState samLinearBorder0
{
    Filter = MIN_MAG_MIP_LINEAR;
    AddressU = Border;
    AddressV = Border;
    BorderColor = float4(0, 0, 0, 0);
};

SamplerState samLinearBorder1
{
    Filter = MIN_MAG_MIP_LINEAR;
    AddressU = Border;
    AddressV = Border;
    BorderColor = float4(1, 1, 1, 1);
};

SamplerState samComparison
{
    Filter = COMPARISON_MIN_MAG_LINEAR_MIP_POINT;
    AddressU = Border;
    AddressV = Border;
    ComparisonFunc = GREATER;
};

DepthStencilState DSS_EnableDepthEqTest
{
    DepthEnable = true;
    DepthWriteMask = ZERO;
    DepthFunc = EQUAL;
};

// Depth stencil state disabling depth test
DepthStencilState DSS_NoDepthTest
{
    DepthEnable = false;
    DepthWriteMask = ZERO;
};

DepthStencilState DSS_NoDepthTest_IncrStencil
{
    DepthEnable = false;
    DepthWriteMask = ZERO;
    StencilEnable = true;
    FrontFaceStencilFunc = ALWAYS;
    BackFaceStencilFunc = ALWAYS;
    FrontFaceStencilPass = INCR;
    BackFaceStencilPass = INCR;
};

DepthStencilState DSS_NoDepthTest_StEqual_IncrStencil
{
    DepthEnable = false;
    DepthWriteMask = ZERO;
    StencilEnable = true;
    FrontFaceStencilFunc = EQUAL;
    BackFaceStencilFunc = EQUAL;
    FrontFaceStencilPass = INCR;
    BackFaceStencilPass = INCR;
    FrontFaceStencilFail = KEEP;
    BackFaceStencilFail = KEEP;
};

DepthStencilState DSS_NoDepthTest_StEqual_KeepStencil
{
    DepthEnable = false;
    DepthWriteMask = ZERO;
    StencilEnable = true;
    FrontFaceStencilFunc = EQUAL;
    BackFaceStencilFunc = EQUAL;
    FrontFaceStencilPass = KEEP;
    BackFaceStencilPass = KEEP;
    FrontFaceStencilFail = KEEP;
    BackFaceStencilFail = KEEP;
};

RasterizerState RS_SolidFill_NoCull
{
    FILLMODE = Solid;
    CullMode = NONE;
};

BlendState NoBlending
{
    BlendEnable[0] = FALSE;
    BlendEnable[1] = FALSE;
    BlendEnable[2] = FALSE;
};

struct DensityProfileLayer
{
    float exp_term;
    float exp_scale;
    float linear_term;
    float const_term;
};

struct AtmosphereParams
{
    float3 solar_irradiance;
    float bottom_radius;
    
    float3 rayleigh_scattering;
    float top_radius;
    
    float3 mie_scattering;
    float mie_g;
    
    float3 mie_extinction;
    float ground_albedo;
    
    float3 absorption_extinction;
    float ozone_width;

    float sun_angular_radius;
    float mu_s_min;
    float nu_power;
    float exposure;

    DensityProfileLayer rayleigh_density;
    DensityProfileLayer mie_density;
    DensityProfileLayer ozone_density[2];
};

cbuffer cbAtmosphereParams
{
    AtmosphereParams atmosphere;
};

struct MiscDynamicParams
{
    float2 f2WQ;
    float scatter_order;
    uint uiMinMaxLevelMax;

    uint4 ui4SrcDstMinMaxOffset;

    float fEnableLightShaft;
    float fIsLightInSpaceCorrect;
    float2 padding;
};

cbuffer cbMiscDynamicParams
{
    MiscDynamicParams misc;
};

struct CameraParams
{
    float3 f3CameraPos;
    float fNearZ;

    float3 f3CameraDir;
    float fFarZ;

    float4x4 View;
    float4x4 Proj;
    float4x4 ViewProj;
    float4x4 InvViewProj;
};

cbuffer cbCameraParams
{
    CameraParams camera;
};

struct LightParams
{
    float3 f3LightDir;
    float padding;

    float4 f4LightScreenPos;

    matrix View;
    matrix Proj;
    matrix ViewProj;
    matrix InvViewProj;
};

cbuffer cbLightParams
{
    LightParams light;
};

struct VertexIn
{
    float3 PosL : POSITION;
    float3 Normal : NORMAL;
    float3 Tangent : TANGENT;
    float2 Tex : TEXCOORD;
};

struct VertexOut
{
    float4 PosH : SV_POSITION;
    float2 Tex : TEXCOORD0;
    float3 PosW : POSITION1;
};

struct QuadVertexOut
{
    float4 m_f4Pos : SV_Position;
    float2 m_f2PosPS : PosPS; // Position in projection space [-1,1]x[-1,1]
    float m_fInstID : InstanceID;
};

cbuffer cbTextureDim
{
    uint SCREEN_WIDTH;
    uint SCREEN_HEIGHT;

    uint TRANSMITTANCE_TEXTURE_WIDTH; //mu
    uint TRANSMITTANCE_TEXTURE_HEIGHT; //r

    uint SCATTERING_TEXTURE_R_SIZE;
    uint SCATTERING_TEXTURE_MU_SIZE;
    uint SCATTERING_TEXTURE_MU_S_SIZE;
    uint SCATTERING_TEXTURE_NU_SIZE;

    uint SCATTERING_TEXTURE_WIDTH;
    uint SCATTERING_TEXTURE_HEIGHT;
    uint SCATTERING_TEXTURE_DEPTH;

    uint IRRADIANCE_TEXTURE_WIDTH;
    uint IRRADIANCE_TEXTURE_HEIGHT;

    uint EPIPOLAR_SLICE_NUM;
    uint EPIPOLAR_SAMPLE_NUM;

    uint SHADOWMAP_TEXTURE_DIM;
    uint MIN_MAX_TEXTURE_DIM;
};

Texture2D<float3> g_tex2DTransmittanceLUT;
Texture2D<float3> g_tex2DOpticalLengthLUT;
Texture2D<float3> g_tex2DDirectIrradianceLUT;
Texture2D<float3> g_tex2DIndirectIrradianceLUT;

Texture3D<float3> g_tex3DSingleScatteringLUT;
Texture3D<float3> g_tex3DMultiScatteringLUT;

Texture3D<float3> g_tex3DSingleMieScatteringLUT;
Texture3D<float4> g_tex3DSingleScatteringCombinedLUT;
Texture3D<float4> g_tex3DMultiScatteringCombinedLUT;

Texture2D<float4> g_tex2DEarthGround;

Texture2D<float>  g_tex2DSpaceDepth;
Texture2D<float>  g_tex2DShadowMap;
Texture2D<float>  g_tex2DSpaceLinearDepth;
Texture2D<float> g_tex2DEpipolarSampleCamDepth;

Texture2D<float2>  g_tex2DMinMaxMipMap;
Texture2D<float2> g_tex2DEpipolarSample;

Texture2D<uint2> g_tex2DInterpolationSample;

Texture2D<float4> g_tex2DUnshadowedSampleScatter;

Texture2D<float4> g_tex2DSliceEnd;
Texture2D<float4> g_tex2DSliceUVOrigDir;
Texture2D<float4> g_tex2DSampleScatter;
Texture2D<float4> g_tex2DInterpolatedScatter;

Texture2D<float3> g_tex2DColorBuffer;

Texture3D<float4> g_tex3DPerlinWorleyNoise;
Texture3D<float3> g_tex3DWorleyNoise;

QuadVertexOut GenerateScreenSizeQuadVS(in uint VertexId : SV_VertexID,
                                                 in uint InstID : SV_InstanceID)
{
    float4 MinMaxUV = float4(-1, -1, 1, 1);
    
    QuadVertexOut Verts[4] =
    {
        { float4(MinMaxUV.xy, 1.0, 1.0), MinMaxUV.xy, InstID },
        { float4(MinMaxUV.xw, 1.0, 1.0), MinMaxUV.xw, InstID },
        { float4(MinMaxUV.zy, 1.0, 1.0), MinMaxUV.zy, InstID },
        { float4(MinMaxUV.zw, 1.0, 1.0), MinMaxUV.zw, InstID }
    };

    return Verts[VertexId];
}

float SafeSqrt(float x)
{
    return sqrt(max(x, 1e-20));
}

float CosClamp(float cos)
{
    return clamp(cos, 0.0, 1.0);
}

float2 ProjToUV(in float2 f2ProjSpaceXY)
{
    return float2(0.5, 0.5) + float2(0.5, -0.5) * f2ProjSpaceXY;
}

float2 UVToProj(in float2 f2UV)
{
    return float2(-1.0, 1.0) + float2(2.0, -2.0) * f2UV;
}

float GetTextureCoordFromUnitRange(float x, int texture_size)
{
    return 0.5 / float(texture_size) + x * (1.0 - 1.0 / float(texture_size));
}

float GetUnitRangeFromTextureCoord(float u, int texture_size)
{
    return (u - 0.5 / float(texture_size)) / (1.0 - 1.0 / float(texture_size));
}

float3 Uncharted2Tonemap(float3 x)
{
    // http://www.gdcvault.com/play/1012459/Uncharted_2__HDR_Lighting
    // http://filmicgames.com/archives/75 - the coefficients are from here
    float A = 0.15; // Shoulder Strength
    float B = 0.50; // Linear Strength
    float C = 0.10; // Linear Angle
    float D = 0.20; // Toe Strength
    float E = 0.02; // Toe Numerator
    float F = 0.30; // Toe Denominator
    return ((x * (A * x + C * B) + D * E) / (x * (A * x + B) + D * F)) - E / F; // E/F = Toe Angle
}

#define RGB_TO_LUMINANCE float3(0.212671, 0.715160, 0.072169)

float3 ToneMap(in float3 f3Color)
{
    //float fAveLogLum = GetAverageSceneLuminance();
    float fAveLogLum = 0.1;
    
    const float middleGray = 1.03 - 2 / (2 + log10(fAveLogLum+1));
    //const float middleGray = g_PPAttribs.m_fMiddleGray;

    float fLumScale = middleGray / fAveLogLum;

    f3Color = max(f3Color, 0);
    float fInitialPixelLum = max(dot(RGB_TO_LUMINANCE, f3Color), 1e-10);
    float fScaledPixelLum = fInitialPixelLum * fLumScale;
    float3 f3ScaledColor = f3Color * fLumScale;

    float whitePoint = 3.0;
    float m_fLuminanceSaturation = 1.0;

    //float fToneMappedLum = 1.0 - exp(-fScaledPixelLum);
    //return fToneMappedLum * pow(f3Color / fInitialPixelLum, m_fLuminanceSaturation);

    float ExposureBias = 10.0f;
    float3 curr = Uncharted2Tonemap(ExposureBias * f3ScaledColor);
    float3 whiteScale = 1.0f / Uncharted2Tonemap(whitePoint);
    return curr * whiteScale;
}

bool IsValidScreenLocation(in float2 f2XY)
{
    const float SAFETY_EPSILON = 0.2f;
    return all(abs(f2XY) <= 1.f - (1.f - SAFETY_EPSILON) / float2(SCREEN_WIDTH, SCREEN_HEIGHT));
}

float ReMap(float fVal,float fMin_old,float fMax_old,float fMin_new,float fMax_new)
{
    return (fVal - fMin_old) / (fMax_old - fMin_old) * (fMax_new - fMin_new) + fMin_new;
}

float3 GetRMuMuS(float3 f3Pos, float3 f3ViewRay)
{
    float fHeight = length(f3Pos);
    float fCosZenithAngle = dot(f3Pos, f3ViewRay) / fHeight;
    float fCosSunZenithAngle = dot(f3Pos, light.f3LightDir) / fHeight;

    return float3(fHeight, fCosZenithAngle, fCosSunZenithAngle);
}

float GetRayMarchLen(float4 f4RMuMuSNu,
                        float fRayEndCamDepth,
                        float fViewRayLenInWorldSpace,
                        out bool bIsNoScatter, 
                        out bool bIsMarchToAtmosphere,
                        out bool bIsMarchToEarth,
                        out bool bIsIntersectEarth,
                        out float fDistToAtmosphereNear,
                        out float fDistToAtmosphereFar,
                        out float fDistToEarthNear)
{
    float fRayMarchLen = fViewRayLenInWorldSpace;
    float fRMu = f4RMuMuSNu.x * f4RMuMuSNu.y;
    float fRSqr = f4RMuMuSNu.x * f4RMuMuSNu.x;
    float fRMuSqr = fRMu * fRMu;

    float fIntersectAtmosphereDiscriminant = fRMuSqr - fRSqr + atmosphere.top_radius * atmosphere.top_radius;
    bIsNoScatter = f4RMuMuSNu.x > atmosphere.top_radius &&
                      (fIntersectAtmosphereDiscriminant < 0 || (fIntersectAtmosphereDiscriminant >= 0 && f4RMuMuSNu.y > 0));
    //if (bIsNoScatter)
    //{
    //    return fRayMarchLen;
    //}

    float fIntersectAtmosphereSqrt = SafeSqrt(fIntersectAtmosphereDiscriminant);
    fDistToAtmosphereNear = -fRMu - fIntersectAtmosphereSqrt;
    fDistToAtmosphereFar = -fRMu + fIntersectAtmosphereSqrt;

    float fIntersectEarthDiscriminant = fRMuSqr - fRSqr + atmosphere.bottom_radius * atmosphere.bottom_radius;
    fDistToEarthNear = -fRMu - SafeSqrt(fIntersectEarthDiscriminant);

    bIsMarchToAtmosphere = false;
    if (fRayEndCamDepth > camera.fFarZ)
    {
        fRayMarchLen = fDistToAtmosphereFar;
        bIsMarchToAtmosphere = true;
    }
    bIsMarchToEarth = false;
    bIsIntersectEarth = fIntersectEarthDiscriminant >= 0 && f4RMuMuSNu.y < 0;
    if (bIsIntersectEarth)
    {
        fRayMarchLen = min(fRayMarchLen, fDistToEarthNear);
        bIsMarchToAtmosphere = false;
        bIsMarchToEarth = (fViewRayLenInWorldSpace == fDistToEarthNear);
    }

    return fRayMarchLen;

}