static const float Pi = 3.141592654f;

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

SamplerState samLinearClamp
{
    Filter = MIN_MAG_MIP_LINEAR;
    AddressU = Clamp;
    AddressV = Clamp;
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

//cbuffer cbAtmosphere
//{
//    float3 solar_irradiance;
//    float bottom_radius;
    
//    float3 rayleigh_scattering;
//    float top_radius;
    
//    float3 mie_scattering;
//    float mie_g;
    
//    float3 mie_extinction;
//    float ground_albedo;
    
//    float3 absorption_extinction;
//    float ozone_width;

//    DensityProfileLayer rayleigh_density;
//    DensityProfileLayer mie_density;
//    DensityProfileLayer ozone_density[2];
//};

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
    int IsMultiScatter;
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

float RayleighPhaseFunction(float w)
{
    float k = 3.0 / (16.0 * Pi);
    return k * (1.0 + w * w);
}

float MiePhaseFunction(float g, float w)
{
    float k = 3.0 / (8.0 * Pi) * (1.0 - g * g) / (2.0 + g * g);
    return k * (1.0 + w * w) / pow(max(1.0 + g * g - 2.0 * g * w, 1e-20), 1.5);
}

float GetLayerDesity(DensityProfileLayer layer, float altitude)
{
    float desity = layer.exp_term * exp(layer.exp_scale * altitude) + layer.linear_term * altitude + layer.const_term;
    return clamp(desity, 0.f, 1.f);
}

float DistanceToTopAtmosphereBoundary(float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1) + atmosphere.top_radius * atmosphere.top_radius;
    return max(-r * mu + SafeSqrt(discriminant), 1e-20);
}

float DistanceToBottomAtmosphereBoundary(float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1) + atmosphere.bottom_radius * atmosphere.bottom_radius;
    return max(-r * mu - SafeSqrt(discriminant), 1e-20);
}

float3 ComputeOpticalLengthToTopAtmosphereBoundary(float r, float mu)
{
    const int SAMPLE_COUNT = 500;
    float dx = DistanceToTopAtmosphereBoundary(r, mu) / float(SAMPLE_COUNT);

    float3 result = 0.f;
    for (int i = 0; i < SAMPLE_COUNT; i++)
    {
        float d = dx * i;
        float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
		// float mu_d = (d + r * mu) / r_d;
        
        float altitude = r_d - atmosphere.bottom_radius;
        float rayleigh = GetLayerDesity(atmosphere.rayleigh_density, altitude);
        float mie = GetLayerDesity(atmosphere.mie_density, altitude);
        float ozone = altitude < atmosphere.ozone_width ?
                                GetLayerDesity(atmosphere.ozone_density[0], altitude) :
                                  GetLayerDesity(atmosphere.ozone_density[1], altitude);

        float weight = (i == 0 || i == SAMPLE_COUNT - 1) ? 0.5 : 1.0;

        result += weight * float3(rayleigh, mie, ozone);
    }
    result *= dx;
    return result;
}

float3 ComputeTransmittanceToTopAtmosphereBoundary(float r, float mu)
{
    float3 optical_length = ComputeOpticalLengthToTopAtmosphereBoundary(r, mu);
    return exp(-(atmosphere.rayleigh_scattering * optical_length.x +
                 atmosphere.mie_extinction * optical_length.y +
                 atmosphere.absorption_extinction * optical_length.z));
}

float2 GetRMuFromTransmittanceUV(float2 uv)
{
    //float2 TransmittanceTextureSize = float2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
    //uv /= TransmittanceTextureSize;
    float mu_x = GetUnitRangeFromTextureCoord(uv.x, TRANSMITTANCE_TEXTURE_WIDTH);
    float r_x = GetUnitRangeFromTextureCoord(uv.y, TRANSMITTANCE_TEXTURE_HEIGHT);

    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = H * r_x;
    float r = sqrt(rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
    float d_min = atmosphere.top_radius - r;
    float d_max = rho + H;
    float d = d_min + mu_x * (d_max - d_min);
    float mu = d == 0 ? 1.f : (H * H - rho * rho - d * d) / (2 * r * d);
    
    return float2(r, mu);
}

float2 GetTransmittanceUVFromRMu(float r, float mu)
{
    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = SafeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);

    float r_x = rho / H;
    float d = DistanceToTopAtmosphereBoundary(r, mu);
    float d_min = atmosphere.top_radius - r;
    float d_max = rho + H;
    float mu_x = (d - d_min) / (d_max - d_min);
    
    return float2(GetTextureCoordFromUnitRange(mu_x, TRANSMITTANCE_TEXTURE_WIDTH),
                    GetTextureCoordFromUnitRange(r_x, TRANSMITTANCE_TEXTURE_HEIGHT));

}

float3 ComputeTransmittanceToTopAtmosphereBoundaryTexture(QuadVertexOut In) : SV_Target
{
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float2 RMu = GetRMuFromTransmittanceUV(f2UV);
    
    return ComputeTransmittanceToTopAtmosphereBoundary(RMu.x, RMu.y);
}

technique11 ComputeTransmittanceTex2DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeTransmittanceToTopAtmosphereBoundaryTexture()));
    }
}

float3 GetTransmittanceToTopAtmosphereBoundary(float r, float mu)
{
    float2 uv = GetTransmittanceUVFromRMu(r, mu);
    return g_tex2DTransmittanceLUT.Sample(samLinearClamp, uv);
}

float3 GetTransmittance(float r, float mu, float r_d, float mu_d, bool ray_r_mu_intersects_ground)
{
    if (ray_r_mu_intersects_ground)
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(r_d, -mu_d) /
                max(GetTransmittanceToTopAtmosphereBoundary(r, -mu), 1e-20),
                float3(1.f, 1.f, 1.f));
    }
    else
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(r, mu) /
                max(GetTransmittanceToTopAtmosphereBoundary(r_d, mu_d), 1e-20),
                float3(1.f, 1.f, 1.f));
    }
}

float3 GetTransmittanceToSun(float r, float mu_s)
{
    float sin = atmosphere.bottom_radius / r;
    float cos = -SafeSqrt(1.0 - sin * sin);
    return GetTransmittanceToTopAtmosphereBoundary(r, mu_s) *
            smoothstep(-sin * atmosphere.sun_angular_radius,
                        sin * atmosphere.sun_angular_radius,
                        mu_s - cos);
}

void GetRMuMuSNuFromUVWQ(in float4 uvwq, out float r, out float mu,
                                out float mu_s, out float nu, out bool ray_r_mu_intersects_ground)
{
    float u_r = GetUnitRangeFromTextureCoord(uvwq.x, SCATTERING_TEXTURE_R_SIZE);
    float u_mu = GetUnitRangeFromTextureCoord(1 - 2 * uvwq.y, SCATTERING_TEXTURE_MU_SIZE / 2);
    float u_mu_s = GetUnitRangeFromTextureCoord(uvwq.z, SCATTERING_TEXTURE_MU_S_SIZE);
    float u_nu = GetUnitRangeFromTextureCoord(uvwq.w, SCATTERING_TEXTURE_NU_SIZE);
    //float u_nu = uvwq.w;

    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
                   atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = u_r * H;

    r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
    
    if (uvwq.y < 0.5) // [-1,mu_horizon]
    {
        float d_min = r - atmosphere.bottom_radius;
        float d_max = rho;
        float d = u_mu * (d_max - d_min) + d_min;
        mu = d == 0.f ? -1.f :
                -(d * d + rho * rho) / (2.f * r * d);
        ray_r_mu_intersects_ground = true;
    }
    else
    {
        float d_min = atmosphere.top_radius;
        float d_max = rho + H;
        float d = u_mu * (d_max - d_min) + d_min;
        mu = d == 0.f ? 1.f :
                -(d * d + H * H - rho * rho) / (2.f * r * d);
        ray_r_mu_intersects_ground = false;
    }
    mu = clamp(mu, -1.f, 1.f);

    //float d_min = atmosphere.top_radius - atmosphere.bottom_radius;
    //float d_max = H;
    //float A = -2.f * atmosphere.mu_s_min * atmosphere.bottom_radius / (d_max - d_min);
    //float a = (A - u_mu_s * A) / max(u_mu_s * A + 1.f, 1e-20);
    //float d = max(a * (d_max - d_min) + d_min, 1e-20);
    //mu_s = clamp(-(d * d - H * H) / (2.f * d * atmosphere.bottom_radius), -1.f, 1.f);

    //float a = 1 - exp(-2.8 * mu_s - 0.8);
    float A = 1 - exp(-3.6);
    float a = A * u_mu_s;
    mu_s = clamp(-(log(1 - a) + 0.8) / 2.8, -1.f, 1.f);
    
    nu = clamp((u_nu - 0.5f) * 2.f, -1.f, 1.f);

}

float4 GetUVWQFromrRMuMuSNu(in float r, in float mu,
                                in float mu_s, in float nu)
{
    float u_r, u_mu, u_mu_s, u_nu;
    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = SafeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);

    u_r = GetTextureCoordFromUnitRange(rho / H, SCATTERING_TEXTURE_R_SIZE);

    float r_mu = r * mu;
    float discriminant = r_mu * r_mu - r * r + atmosphere.bottom_radius * atmosphere.bottom_radius;

    if (discriminant > 0.f && r_mu < 0) // [-1,mu_horizon]
    {
        float d = -r_mu - SafeSqrt(discriminant);
        float d_min = r - atmosphere.bottom_radius;
        float d_max = rho;
        u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(
                            d_max == d_min ? 0.f : (d - d_min) / (d_max - d_min),
                            SCATTERING_TEXTURE_MU_SIZE / 2);
    }
    else // [mu_horizon,1]
    {
        float d = -r_mu + SafeSqrt(discriminant + H * H); //distance to top atmosphere
        float d_min = atmosphere.top_radius - r;
        float d_max = rho + H;
        u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(
                            (d - d_min) / (d_max - d_min),
                            SCATTERING_TEXTURE_MU_SIZE / 2);
    }
    float a = 1 - exp(-2.8 * mu_s - 0.8);
    float A = 1 - exp(-3.6);
    u_mu_s = GetTextureCoordFromUnitRange(max(a / A, 1e-20), SCATTERING_TEXTURE_MU_S_SIZE);

    //float d = DistanceToTopAtmosphereBoundary(atmosphere.bottom_radius, mu_s);
    //float d_min = atmosphere.top_radius - atmosphere.bottom_radius;
    //float d_max = H;
    //float a = (d - d_min) / (d_max - d_min);
    //float A = -2.f * atmosphere.mu_s_min * atmosphere.bottom_radius / (d_max - d_min);
    //u_mu_s = GetTextureCoordFromUnitRange(
    //            max(1.f - a / A, 1e-20) / (1.f + a),
    //            SCATTERING_TEXTURE_MU_S_SIZE);

    u_nu = GetTextureCoordFromUnitRange((1.f + nu) / 2.f, SCATTERING_TEXTURE_NU_SIZE);
    //u_nu = (1.f + nu) / 2.f;

    return float4(u_r, u_mu, u_mu_s, u_nu);
}

struct SingleScatterOutput
{
    float4 single_scatter : SV_Target0;
    float4 scatter_combined : SV_Target1;
    float4 scatter_mie : SV_Target2;
};

SingleScatterOutput ComputeSingleScatteringTexture(QuadVertexOut In) : SV_Target
{
    SingleScatterOutput res;
    float r, mu, mu_s, nu;
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    bool ray_r_mu_intersects_ground;

    GetRMuMuSNuFromUVWQ(float4(f2UV, misc.f2WQ), r, mu, mu_s, nu, ray_r_mu_intersects_ground);

    const int SAMPLE_COUNT = 50;

    float dx = ray_r_mu_intersects_ground ? DistanceToBottomAtmosphereBoundary(r, mu) :
                                                DistanceToTopAtmosphereBoundary(r, mu);
    dx /= SAMPLE_COUNT;

    float3 rayleigh = 0.f;
    float3 mie = 0.f;
    for (int i = 0; i < SAMPLE_COUNT; ++i)
    {
        float d = dx * i;
        float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
        float mu_d = CosClamp((d + r * mu) / r_d);
        float mu_s_d = CosClamp((r * mu_s + d * nu) / r_d);
        float altitude_d = r_d - atmosphere.bottom_radius;

        float3 transmittance = GetTransmittance(r, mu, r_d, mu_d, ray_r_mu_intersects_ground) *
                                GetTransmittanceToSun(r_d, mu_s_d);
        float weight = (i == 0 || i == SAMPLE_COUNT) ? 0.5f : 1.0f;
        
        rayleigh += weight * transmittance * GetLayerDesity(atmosphere.rayleigh_density, altitude_d);
        mie += weight * transmittance * GetLayerDesity(atmosphere.mie_density, altitude_d);
    }
    rayleigh *= (dx * atmosphere.solar_irradiance * atmosphere.rayleigh_scattering);
    mie *= (dx * atmosphere.solar_irradiance * atmosphere.mie_scattering);

    res.single_scatter = float4(rayleigh * RayleighPhaseFunction(nu) + mie * MiePhaseFunction(nu, atmosphere.mie_g), 1);
    res.scatter_combined = float4(rayleigh, mie.r);
    res.scatter_mie = float4(mie, 1);
    return res;
}

technique11 ComputeSingleScatterTex3DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeSingleScatteringTexture()));
    }
}

void GetScatteringUVW(float r, float mu, float mu_s, float nu, out float3 f3UVW0, out float3 f3UVW1, out float fQWeight)
{
    float4 f4UVWQ = GetUVWQFromrRMuMuSNu(r, mu, mu_s, nu);
    f3UVW0 = f4UVWQ.xyz;
    float fQ0Slice = floor(f4UVWQ.w * SCATTERING_TEXTURE_NU_SIZE - 0.5);
    fQ0Slice = clamp(fQ0Slice, 0, SCATTERING_TEXTURE_NU_SIZE - 1);
    fQWeight = (f4UVWQ.w * SCATTERING_TEXTURE_NU_SIZE - 0.5) - fQ0Slice;
    fQWeight = max(fQWeight, 0.f);
    f3UVW0.z = (fQ0Slice + f3UVW0.z) / SCATTERING_TEXTURE_NU_SIZE;

    float fQ1Slice = min(fQ0Slice + 1, SCATTERING_TEXTURE_NU_SIZE - 1);
    float fSliceOffset = (fQ1Slice - fQ0Slice) / SCATTERING_TEXTURE_NU_SIZE;
    f3UVW1 = f3UVW0 + float3(0, 0, fSliceOffset);
}

float3 GetScattering(float r, float mu, float mu_s, float nu, int scatter_order)
{
    float3 f3UVW0, f3UVW1;
    float fQWeight;
    GetScatteringUVW(r, mu, mu_s, nu, f3UVW0, f3UVW1, fQWeight);
    
    float3 f3Insctr0, f3Insctr1;
    if (1 == scatter_order)
    {
        f3Insctr0 = g_tex3DSingleScatteringLUT.SampleLevel(samLinearClamp, f3UVW0, 0);
        f3Insctr1 = g_tex3DSingleScatteringLUT.SampleLevel(samLinearClamp, f3UVW1, 0);
    }
    else
    {
        f3Insctr0 = g_tex3DMultiScatteringLUT.SampleLevel(samLinearClamp, f3UVW0, 0);
        f3Insctr1 = g_tex3DMultiScatteringLUT.SampleLevel(samLinearClamp, f3UVW1, 0);
    }
    float3 f3Inscattering = lerp(f3Insctr0, f3Insctr1, fQWeight);

    return f3Inscattering;
}

float4 GetScatteringCombined(float r, float mu, float mu_s, float nu, int scatter_order)
{
    float3 f3UVW0, f3UVW1;
    float fQWeight;
    GetScatteringUVW(r, mu, mu_s, nu, f3UVW0, f3UVW1, fQWeight);
    
    float4 f4Insctr0, f4Insctr1;
    if (1 == scatter_order)
    {
        f4Insctr0 = g_tex3DSingleScatteringCombinedLUT.SampleLevel(samLinearClamp, f3UVW0, 0);
        f4Insctr1 = g_tex3DSingleScatteringCombinedLUT.SampleLevel(samLinearClamp, f3UVW1, 0);
    }
    else
    {
        f4Insctr0 = g_tex3DMultiScatteringCombinedLUT.SampleLevel(samLinearClamp, f3UVW0, 0);
        f4Insctr1 = g_tex3DMultiScatteringCombinedLUT.SampleLevel(samLinearClamp, f3UVW1, 0);
    }

    float4 f4Inscattering = lerp(f4Insctr0, f4Insctr1, fQWeight);

    return f4Inscattering;
}

float3 GetScatteringMie(float r, float mu, float mu_s, float nu)
{
    float3 f3UVW0, f3UVW1;
    float fQWeight;
    GetScatteringUVW(r, mu, mu_s, nu, f3UVW0, f3UVW1, fQWeight);
    
    float3 f3Insctr0, f3Insctr1;

    f3Insctr0 = g_tex3DSingleMieScatteringLUT.SampleLevel(samLinearClamp, f3UVW0, 0);
    f3Insctr1 = g_tex3DSingleMieScatteringLUT.SampleLevel(samLinearClamp, f3UVW1, 0);
    
    float3 f4Inscattering = lerp(f3Insctr0, f3Insctr1, fQWeight);

    return f4Inscattering;
}

float2 GetIrradianceUVFromRMuS(float r, float mu_s)
{
    float u_x = (1.0 + mu_s) / 2.0;
    float v_x = (r - atmosphere.bottom_radius) / (atmosphere.top_radius - atmosphere.bottom_radius);
    return float2(GetTextureCoordFromUnitRange(u_x, IRRADIANCE_TEXTURE_WIDTH),
                    GetTextureCoordFromUnitRange(v_x, IRRADIANCE_TEXTURE_HEIGHT));

}

float3 ComputeScatteringDensity(float r, float mu, float mu_s, float nu, int scatter_order)
{
    const int SAMPLE_COUNT = 16;
    const float dtheta = Pi / SAMPLE_COUNT;
    const float dphi = Pi / SAMPLE_COUNT;

    float3 zenith = float3(0.f, 1.f, 0.f);
    float3 omega = float3(sqrt(1.f - mu * mu), mu, 0.f);
    float sun_dir_x = omega.x == 0.f ? 0.f : (nu - mu * mu_s) / omega.x;
    float sun_dir_z = sqrt(1 - sun_dir_x * sun_dir_x * mu_s * mu_s);
    float3 omega_s = float3(sun_dir_x, mu_s, sun_dir_z);

    float3 rayleigh = 0.f;
    float3 mie = 0.f;

    float altitude = r - atmosphere.bottom_radius;
    float rayleigh_density = GetLayerDesity(atmosphere.rayleigh_density, altitude);
    float mie_density = GetLayerDesity(atmosphere.mie_density, altitude);

    for (int m = 0; m < SAMPLE_COUNT; ++m)
    {
        float theta = (m + 0.5) * dtheta;
        float cos_theta = cos(theta);
        float sin_theta = sin(theta);
        float dSolidAngle = sin_theta * dphi * dtheta;

        float distance_to_ground = 0.f;
        float mu_ground = 0.f;
        float3 transmittance_to_ground = 0.f;
        float ground_albedo = 0.f;
        float discriminant = r * r * (mu * mu - 1) + atmosphere.bottom_radius * atmosphere.bottom_radius;
        bool ray_r_theta_intersects_ground = discriminant >= 0.f && mu < 0;
        if (ray_r_theta_intersects_ground)
        {
            distance_to_ground = -r * mu - sqrt(discriminant);
            mu_ground = (r * mu + distance_to_ground) / atmosphere.bottom_radius;
            transmittance_to_ground = GetTransmittance(r, mu, atmosphere.bottom_radius, mu_ground, true);
            ground_albedo = atmosphere.ground_albedo;
        }

        for (int n = 0; n < SAMPLE_COUNT * 2; ++n)
        {
            float phi = (n + 0.5) * dphi;
            
            float3 omega_i = float3(sin_theta * cos(phi), cos_theta, sin_theta * sin(phi));
            float nu_i = dot(omega_s, omega_i);
            float3 incident_radiance = GetScattering(r, omega_i.y, mu_s, nu_i, scatter_order - 1);
           
            float3 ground_normal = normalize(zenith * r + omega_i * distance_to_ground);
            float w = dot(omega, omega_i);
            float2 irradiance_uv = GetIrradianceUVFromRMuS(atmosphere.bottom_radius, nu_i);
            float3 ground_irradiance = g_tex2DIndirectIrradianceLUT.Sample(samLinearClamp, irradiance_uv) * dot(ground_normal, omega_s);
            incident_radiance += transmittance_to_ground * ground_albedo * (1.0 / Pi) * ground_irradiance;

            rayleigh += incident_radiance * RayleighPhaseFunction(w) * dSolidAngle;
            mie += incident_radiance * MiePhaseFunction(atmosphere.mie_g, w) * dSolidAngle;
        }
    }
    rayleigh *= rayleigh_density * atmosphere.rayleigh_scattering;
    mie *= mie_density * atmosphere.mie_scattering;

    return rayleigh + mie;
}

struct MultiScatterOutput
{
    float4 multi_scatter : SV_Target0;
    float4 scatter_combined : SV_Target1;
};

MultiScatterOutput ComputeMultiScatteringTexture(QuadVertexOut In) : SV_Target
{
    MultiScatterOutput multi_scatter;
    float3 rayleigh_mie = 0.f;
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float r, mu, mu_s, nu;
    bool ray_r_mu_intersects_ground;
    GetRMuMuSNuFromUVWQ(float4(f2UV, misc.f2WQ), r, mu, mu_s, nu, ray_r_mu_intersects_ground);

    const int SAMPLE_COUNT = 50;

    float dx = ray_r_mu_intersects_ground ? DistanceToBottomAtmosphereBoundary(r, mu) :
                                                DistanceToTopAtmosphereBoundary(r, mu);
    dx /= SAMPLE_COUNT;
    
    for (int i = 0; i < SAMPLE_COUNT; ++i)
    {
        float d = dx * i;
        float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
        float mu_d = CosClamp((r * mu + d) / r_d);
        float mu_s_d = CosClamp((r * mu_s + d * nu) / r_d);

        float altitude = r_d - atmosphere.bottom_radius;
        float3 transmittance = GetTransmittance(r, mu, r_d, mu_d, ray_r_mu_intersects_ground);
        float3 scatter_density = ComputeScatteringDensity(r_d, mu_d, mu_s_d, nu, misc.scatter_order);
        
        float weight = (i == 0 || i == SAMPLE_COUNT - 1) ? 0.5 : 1.0;

        float3 res = weight * transmittance * scatter_density;

        rayleigh_mie += res;
    }
    rayleigh_mie *= dx;
    float w = RayleighPhaseFunction(nu);
    multi_scatter.multi_scatter = float4(rayleigh_mie, 1);
    multi_scatter.scatter_combined = float4(w == 0 ?  float3(0,0,0): rayleigh_mie / w,0) +
                                     GetScatteringCombined(r, mu, mu_s, nu, 1);
    return multi_scatter;
}

technique11 ComputeMultiScatterTex3DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeMultiScatteringTexture()));
    }
}

void GetRMuSFromIrradianceUV(float2 uv, out float r, out float mu_s)
{
    float mu_s_x = GetUnitRangeFromTextureCoord(uv.x, IRRADIANCE_TEXTURE_WIDTH);
    float r_x = GetUnitRangeFromTextureCoord(uv.y, IRRADIANCE_TEXTURE_HEIGHT);

    r = r_x * (atmosphere.top_radius - atmosphere.bottom_radius) + atmosphere.bottom_radius;
    mu_s = mu_s_x * 2.0 - 1.0;
}

float3 ComputeDirectIrradianceTexture(QuadVertexOut In) : SV_Target
{
    float r, mu_s;
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    GetRMuSFromIrradianceUV(f2UV, r, mu_s);

    float alpha_s = atmosphere.sun_angular_radius;
    // Approximate average of the cosine factor mu_s over the visible fraction of
    // the Sun disc.
    float average_cosine_factor =
                mu_s < -alpha_s ? 0.0 : (mu_s > alpha_s ? mu_s :
                    (mu_s + alpha_s) * (mu_s + alpha_s) / (4.0 * alpha_s));

    return atmosphere.solar_irradiance *
                GetTransmittanceToTopAtmosphereBoundary(r, mu_s) * average_cosine_factor;

}

technique11 ComputeDirectIrradiance2DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeDirectIrradianceTexture()));
    }
}

float3 ComputeIndirectIrradianceTexture(QuadVertexOut In) : SV_Target
{
    float3 result = 0.f;
    float r, mu_s;
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    GetRMuSFromIrradianceUV(f2UV, r, mu_s);

    const int SAMPLE_COUNT = 32;
    float dphi = Pi / SAMPLE_COUNT;
    float dtheta = Pi / SAMPLE_COUNT;

    float3 omega_s = float3(sqrt(1 - mu_s * mu_s), mu_s, 0.f);
    for (int m = 0; m < SAMPLE_COUNT / 2; ++m)
    {
        float theta = (m + 0.5) * dtheta;
        float theta_sin = sin(theta);
        float theta_cos = cos(theta);
        float dSolidAngle = theta_sin * dtheta * dphi;

        for (int n = 0; n < SAMPLE_COUNT * 2; ++n)
        {
            float phi = (n + 0.5) * dphi;
            float3 omega_i = float3(theta_sin * cos(phi), theta_cos, theta_sin * sin(phi));
            float nu = dot(omega_i, omega_s);
            float3 scatter = GetScattering(r, omega_i.z, mu_s, nu, misc.scatter_order);
            result += scatter * dSolidAngle * omega_i.z;
        }
    }
    return result;
}


technique11 ComputeIndirectIrradiance2DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeIndirectIrradianceTexture()));
    }
}

//struct QuadViewRayOut
//{
//    float4 m_f4Pos : SV_Position;
//    float2 m_f2PosPS : PosPS; // Position in projection space [-1,1]x[-1,1]
//    float m_fInstID : InstanceID;
//    float3 m_f4Ray : ViewRay;
//};

//QuadViewRayOut GenerateViewRayVS(in uint VertexId : SV_VertexID,
//                                                 in uint InstID : SV_InstanceID)
//{
//    float4 MinMaxUV = float4(-1, -1, 1, 1);
    
//    QuadViewRayOut ViewRayOut[4] =
//    {
//        { float4(MinMaxUV.xy, 1.0, 1.0), MinMaxUV.xy, InstID, float3(0.f, 0.f, 0.f) },
//        { float4(MinMaxUV.xw, 1.0, 1.0), MinMaxUV.xw, InstID, float3(0.f, 0.f, 0.f) },
//        { float4(MinMaxUV.zy, 1.0, 1.0), MinMaxUV.zy, InstID, float3(0.f, 0.f, 0.f) },
//        { float4(MinMaxUV.zw, 1.0, 1.0), MinMaxUV.zw, InstID, float3(0.f, 0.f, 0.f) }
//    };
//    ViewRayOut[VertexId].m_f4Ray = mul(ViewRayOut[VertexId].m_f4Pos, InvViewProj).xyz;
//    return ViewRayOut[VertexId];

//}

float3 GetIrradianceDirectFromSun(float r, float mu_s, float nu)
{
    float3 irradiance_from_sun = 0.f;
    float3 transmmitance_to_sun;
    if (r > atmosphere.top_radius)
    {
        transmmitance_to_sun = 1.f;
        irradiance_from_sun = transmmitance_to_sun * atmosphere.solar_irradiance * max(nu, 0.f);
    }
    else
    {
        float2 irradiance_uv = GetIrradianceUVFromRMuS(r, mu_s);
        transmmitance_to_sun = GetTransmittanceToSun(r, mu_s);
        float alpha_s = atmosphere.sun_angular_radius;
        //Approximate average of the cosine factor mu_s over the visible fraction of
        //the Sun disc.
        float average_cosine_factor =
                mu_s < -alpha_s ? 0.0 : (mu_s > alpha_s ? mu_s :
                    (mu_s + alpha_s) * (mu_s + alpha_s) / (4.0 * alpha_s));
        irradiance_from_sun = atmosphere.solar_irradiance *
                GetTransmittanceToTopAtmosphereBoundary(r, mu_s) * average_cosine_factor * max(nu, 0.f);
        //irradiance_from_sun = g_tex2DDirectIrradianceLUT.Sample(samLinearClamp, irradiance_uv) * max(nu, 0.f);
    }
    if (nu > cos(atmosphere.sun_angular_radius))
        irradiance_from_sun += transmmitance_to_sun * atmosphere.solar_irradiance /
                                (Pi * atmosphere.sun_angular_radius * atmosphere.sun_angular_radius);
    return irradiance_from_sun;
}

float3 GetSkyMultiScatter(float r, float mu, float mu_s, float nu)
{
    float4 scatter_combined = GetScatteringCombined(r, mu, mu_s, nu, misc.IsMultiScatter);
    
    float3 mie_scatter = scatter_combined.r == 0 ? float3(0,0,0):
        scatter_combined.rgb * scatter_combined.a / scatter_combined.r *
	    (atmosphere.rayleigh_scattering.r / atmosphere.mie_scattering.r) *
	    (atmosphere.mie_scattering / atmosphere.rayleigh_scattering);
    
    return (mie_scatter * MiePhaseFunction(atmosphere.mie_g, nu) + scatter_combined.rgb * RayleighPhaseFunction(nu));
}

float3 GetSkyMultiScatterToGround(float r, float mu, float mu_s,float r_d, float mu_d, float mu_s_d, 
                                    float nu, float shadow_length, float d, float transmittance)
{
    // S[L]x - T(x,xs)S[L]xs=x0-lv
    float3 scatter = GetSkyMultiScatter(r, mu, mu_s, nu);
    float3 scatter_d;

    float shadow_transmittance;
    if (0.f == shadow_length)
    {
        scatter_d = GetSkyMultiScatter(r_d, mu_d, mu_s_d, nu);
        shadow_transmittance = transmittance;
    }
    else
    {
        float d = max(d - shadow_length, 0.f);
        float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
        float mu_d = (r * mu + d) / r_d;
        float mu_s_d = (r * mu_s + d * nu) / r_d;
        scatter_d = GetSkyMultiScatter(r_d, mu_d, mu_s_d, nu);
        shadow_transmittance = GetTransmittance(r, mu, r_d, mu_d, true);
    }
    return scatter - shadow_length * scatter_d;
}

float3 GetSkyMultiScatterToAtmosphere(float r, float mu, float mu_s, float nu, float shadow_length, float d)
{
    // T(x,xs)S[L]xs=x+lv
    float3 scatter_d;
    float shadow_transmittance;
    if (0.f == shadow_length)
    {
        return GetSkyMultiScatter(r, mu, mu_s, nu);
    }
    else
    {
        float d = shadow_length;
        float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
        float mu_d = (r * mu + d) / r_d;
        float mu_s_d = (r * mu_s + d * nu) / r_d;
        scatter_d = GetSkyMultiScatter(r_d, mu_d, mu_s_d, nu);
        shadow_transmittance = GetTransmittance(r, mu, r_d, mu_d, false);
        return shadow_transmittance * scatter_d;
    }
}

float3 HDR(float3 L)
{
    L = L * misc.exposure;
    L.r = L.r < 1.413 ? pow(L.r * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.r);
    L.g = L.g < 1.413 ? pow(L.g * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.g);
    L.b = L.b < 1.413 ? pow(L.b * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.b);
    return L;
}

float4 DrawGroundAndSky(QuadVertexOut In) : SV_Target
{
    float3 v = mul(float4(In.m_f2PosPS, 1.0f, 1.0f), InvViewProj).xyz;
    float v_length = length(v);
    v = v / v_length;
    
    //float fragment_angular_size =
    //  length(ddx(In.m_f4Ray) + ddy(In.m_f4Ray)) / length(In.m_f4Ray);
    
    float3 camera_pos = misc.f3CameraPos - misc.f3EarthCenter;
    float3 sun_dir = misc.f3SunDir;
    float r = length(camera_pos);
    float mu = dot(camera_pos, v) / r;

    float nu = dot(v, sun_dir);
    float mu_s = dot(camera_pos, sun_dir) / r;

    float3 ground_radiance = 0.f;
    float3 sky_radiance = 0.f;
    float3 irradiance_from_sun = GetIrradianceDirectFromSun(r, mu_s, nu); //L0

    float intersect_atmosphere_discriminant = r * r * (mu * mu - 1) + atmosphere.top_radius * atmosphere.top_radius;
    if (r > atmosphere.top_radius && 
        (intersect_atmosphere_discriminant < 0 || intersect_atmosphere_discriminant >= 0 && mu > 0))
    {
        return float4(HDR(irradiance_from_sun), 1);
    }

    float distance_to_atmosphere_near = -r * mu - sqrt(intersect_atmosphere_discriminant);
    if (distance_to_atmosphere_near > 0.f)
    {
        camera_pos += distance_to_atmosphere_near * v;
        r = atmosphere.top_radius;
        mu = dot(camera_pos, v) / r;
        mu_s = dot(camera_pos, sun_dir) / r;
    }
    
    float intersect_ground_discriminant = r * r * (mu * mu - 1) + misc.f3EarthCenter.y * misc.f3EarthCenter.y;
    
    if (intersect_ground_discriminant >= 0 && mu < 0)
    {
        float distance_to_ground = -r * mu - sqrt(intersect_ground_discriminant);
        float3 pos_d = camera_pos + distance_to_ground * v;
        float r_d = length(pos_d);
        float3 normal_d = pos_d / r_d;
        float mu_d = CosClamp((distance_to_ground + r * mu) / r_d);
        float mu_s_d = CosClamp((distance_to_ground * nu + r * mu_s) / r_d);
        
        float2 coords = float2(normal_d.x == 0 ? 0 : atan2(normal_d.z,normal_d.x),acos(normal_d.y)) * float2(0.5, 1.0) / Pi + float2(0.5, 0.0);
        float4 reflectance = g_tex2DEarthGround.SampleLevel(samLinearClamp, coords,0); //* float4(0.2, 0.2, 0.2, 1.0);
       
        float transmittance = GetTransmittance(r, mu, r_d, mu_d, true);
        
        float2 irradiance_uv = GetIrradianceUVFromRMuS(r_d, mu_s_d);

        // R[L0] + R[L*]
        float3 sunColor = g_tex2DDirectIrradianceLUT.Sample(samLinearClamp, irradiance_uv) * max(mu_s_d, 0.f);
        //float3 skyColor = g_tex2DIndirectIrradianceLUT.Sample(samLinearClamp, irradiance_uv) * (1.0 + dot(normal, pos) / r) * 0.5;
        float3 skyColor = g_tex2DIndirectIrradianceLUT.Sample(samLinearClamp, irradiance_uv);
        
        float shadow_length = 0;
        // S[L]x - T(x,xs)S[L]xs=x0-lv
        float3 scatter = GetSkyMultiScatterToGround(r, mu, mu_s, r_d, mu_d, mu_s_d, nu, shadow_length, distance_to_ground, transmittance);
        sky_radiance = scatter;

        ground_radiance = reflectance.rgb * atmosphere.ground_albedo * (1.0 / Pi) * (sunColor + skyColor) * transmittance;
        
        return float4(HDR(ground_radiance + sky_radiance), 1.0);
    }
    else
    {
        float distance_to_atmosphere = -r * mu + sqrt(intersect_atmosphere_discriminant);
        float3 pos_d = camera_pos + distance_to_atmosphere * v;
        float r_d = length(pos_d);
        float mu_d = CosClamp((distance_to_atmosphere + r * mu) / r_d);
        float mu_s_d = CosClamp((distance_to_atmosphere * nu + r * mu_s) / r_d);
    
        float transmittance = GetTransmittance(r, mu, r_d, mu_d, false);
        
        // T(x,xs)S[L]xs=x+lv
        float3 scatter = GetSkyMultiScatterToAtmosphere(r, mu, mu_s, nu, 0, distance_to_atmosphere);
        sky_radiance = scatter;
        
        return float4(HDR(sky_radiance + irradiance_from_sun), 1);
    }
    return float4(HDR(irradiance_from_sun + ground_radiance + sky_radiance), 1.0f);
}

technique11 DrawGroundAndSkyTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, DrawGroundAndSky()));
    }
}

float4 TestTexture(QuadVertexOut In) : SV_Target
{
    float3 v = mul(float4(In.m_f2PosPS, 1.0f, 1.0f), InvViewProj).xyz;
    float v_length = length(v);
    v = v / v_length;
    
    //float fragment_angular_size =
    //  length(ddx(In.m_f4Ray) + ddy(In.m_f4Ray)) / length(In.m_f4Ray);
    
    float3 camera_pos = misc.f3CameraPos - misc.f3EarthCenter;
    float3 sun_dir = misc.f3SunDir;
    float r = length(camera_pos);
    float mu = dot(camera_pos, v) / r;

    float nu = dot(v, sun_dir);
    float mu_s = dot(camera_pos, sun_dir) / r;

    float3 scatter = GetSkyMultiScatter(r, mu, mu_s, nu);

    return float4(scatter, 1.0);
}

technique11 TestTextureTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, TestTexture()));
    }
}