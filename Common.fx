static const float PI = 3.141592654f;

static const int TRANSMITTANCE_TEXTURE_WIDTH = 256; //mu
static const int TRANSMITTANCE_TEXTURE_HEIGHT = 64; //r

static const int SCATTERING_TEXTURE_R_SIZE = 32;
static const int SCATTERING_TEXTURE_MU_SIZE = 128;
static const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
static const int SCATTERING_TEXTURE_NU_SIZE = 8;

static const int SCATTERING_TEXTURE_WIDTH = SCATTERING_TEXTURE_R_SIZE;
static const int SCATTERING_TEXTURE_HEIGHT = SCATTERING_TEXTURE_MU_SIZE;
static const int SCATTERING_TEXTURE_DEPTH = SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE;

static const int IRRADIANCE_TEXTURE_WIDTH = 64;
static const int IRRADIANCE_TEXTURE_HEIGHT = 16;

// The conversion factor between watts and lumens.
static const float MAX_LUMINOUS_EFFICACY = 683.0;

#define USE_LUT_PARAMETERIZATION 1
#define USE_INTEGRAL_OPTIMIZATION 1

SamplerState samLinearClamp
{
    Filter = MIN_MAG_MIP_LINEAR;
    AddressU = Clamp;
    AddressV = Clamp;
    AddressW = Clamp;
};

// Depth stencil state disabling depth test
DepthStencilState DSS_NoDepthTest
{
    DepthEnable = false;
    DepthWriteMask = ZERO;
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
    float padding[2];

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
    int scatter_order;
    float exposure;

    float3 f3CameraPos;
    float nu_power;
    float3 f3EarthCenter;
    float padding2;
    float3 f3SunDir;
    float padding3;
    float3 f3CameraDir;
    float padding4;
};

cbuffer cbMiscDynamicParams
{
    MiscDynamicParams misc;
};

cbuffer cbMatrix
{
    float4x4 InvViewProj;
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

Texture2D<float3> g_tex2DTransmittanceLUT;
Texture2D<float3> g_tex2DDirectIrradianceLUT;
Texture2D<float3> g_tex2DIndirectIrradianceLUT;

Texture3D<float3> g_tex3DSingleScatteringLUT;
Texture3D<float3> g_tex3DMultiScatteringLUT;

Texture3D<float3> g_tex3DSingleMieScatteringLUT;
Texture3D<float4> g_tex3DSingleScatteringCombinedLUT;
Texture3D<float4> g_tex3DMultiScatteringCombinedLUT;

Texture2D<float4> g_tex2DEarthGround;

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