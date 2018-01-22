static const int TRANSMITTANCE_TEXTURE_WIDTH = 256;    //mu
static const int TRANSMITTANCE_TEXTURE_HEIGHT = 64; //r

static const int SCATTERING_TEXTURE_R_SIZE = 32;
static const int SCATTERING_TEXTURE_MU_SIZE = 128;
static const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
static const int SCATTERING_TEXTURE_NU_SIZE = 8;

static const  int SCATTERING_TEXTURE_WIDTH = SCATTERING_TEXTURE_R_SIZE;
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
};

cbuffer cbMiscDynamicParams
{
    MiscDynamicParams misc;
};

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
	float2 Tex		: TEXCOORD0;
	float3 PosW		: POSITION1;
};

struct QuadVertexOut
{
    float4 m_f4Pos : SV_Position;
    float2 m_f2PosPS : PosPS; // Position in projection space [-1,1]x[-1,1]
    float m_fInstID : InstanceID;
};

Texture2D<float3> g_tex2DTransmittanceLUT;
Texture3D<float4> g_tex3DSingleScatteringLUT;

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

float GetLayerDesity(DensityProfileLayer layer, float altitude)
{
	float desity = layer.exp_term * exp(layer.exp_scale * altitude) + layer.linear_term * altitude + layer.const_term;
	return clamp(desity, 0.f, 1.f);
}

float DistanceToTopAtmosphereBoundary(float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1) + atmosphere.top_radius * atmosphere.top_radius;
    return max(-r * mu + SafeSqrt(discriminant),1e-20);
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

        result += weight * float3(rayleigh,mie,ozone);
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
    float rho = SafeSqrt( r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);

    float r_x = rho / H;
    float d = DistanceToTopAtmosphereBoundary(r, mu);
    float d_min = atmosphere.top_radius - r;
    float d_max = rho + H;
    float mu_x = (d - d_min) / (d_max - d_min);
    
    return float2(GetTextureCoordFromUnitRange(mu_x, TRANSMITTANCE_TEXTURE_WIDTH),
                    GetTextureCoordFromUnitRange(r_x, TRANSMITTANCE_TEXTURE_HEIGHT));

}

float3 ComputeTransmittanceToTopAtmosphereBoundaryTexture(QuadVertexOut In):SV_Target
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
    return g_tex2DTransmittanceLUT.Sample(samLinearClamp,uv);
}

float3 GetTransmittance(float r, float mu, float r_d,float mu_d, bool ray_r_mu_intersects_ground)
{
    if(ray_r_mu_intersects_ground)
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(r_d, -mu_d) /
                GetTransmittanceToTopAtmosphereBoundary(r, mu),
                float3(1.f,1.f,1.f));
    }
    else
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(r, mu) /
                GetTransmittanceToTopAtmosphereBoundary(r_d, mu_d),
                float3(1.f, 1.f, 1.f));
    }
}

float3 GetTransmittanceToSun(float r,float mu_s)
{
    float sin = atmosphere.bottom_radius / r;
    float cos = -SafeSqrt(1.0 - sin * sin);
    return GetTransmittanceToTopAtmosphereBoundary(r, mu_s) *
            smoothstep(-sin * atmosphere.sun_angular_radius,
                        sin * atmosphere.sun_angular_radius,
                        mu_s - cos);          
}

void GetRMuMuSNuFromUVWQ(in float4 uvwq,out float r,out float mu,
                                out float mu_s,out float nu,out bool ray_r_mu_intersects_ground)
{
    float u_r = GetUnitRangeFromTextureCoord(uvwq.x, SCATTERING_TEXTURE_R_SIZE);
    float u_mu = GetUnitRangeFromTextureCoord(1-2*uvwq.y, SCATTERING_TEXTURE_MU_SIZE/2);
    float u_mu_s = GetUnitRangeFromTextureCoord(uvwq.z, SCATTERING_TEXTURE_MU_S_SIZE);
    float u_nu = GetUnitRangeFromTextureCoord(uvwq.w, SCATTERING_TEXTURE_NU_SIZE);

    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
                   atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = u_r * H;

    r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
    
    if (uvwq.y<0.5) // [-1,mu_horizon]
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

    float d_min = atmosphere.top_radius - atmosphere.bottom_radius;
    float d_max = H;
    float A = -2.f * atmosphere.mu_s_min * atmosphere.bottom_radius / (d_max - d_min);
    float a = (A - u_mu_s * A) / max(u_mu_s * A + 1.f, 1e-20);
    float d = max(a * (d_max - d_min) + d_min, 1e-20);
    mu_s = clamp(-(d * d - H * H) / (2.f * d * atmosphere.bottom_radius), -1.f, 1.f);
    
    nu = clamp((u_nu - 0.5f) * 2.f, -1.f, 1.f);

}

float4 GetUVWQFromrRMuMuSNu(in float r,in float mu,
                                in float mu_s,in float nu)
{
    float u_r, u_mu, u_mu_s, u_nu;
    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = SafeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);

    u_r = GetTextureCoordFromUnitRange(rho / H, SCATTERING_TEXTURE_R_SIZE);

    float r_mu = r * mu;
    float discriminant = r_mu * r_mu - r * r + atmosphere.bottom_radius * atmosphere.bottom_radius;

    if(discriminant>0.f) // [-1,mu_horizon]
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
        float d = -r_mu + SafeSqrt(discriminant + H * H);//distance to top atmosphere
        float d_min = atmosphere.top_radius - r;
        float d_max = rho + H;
        u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(
                            (d - d_min) / (d_max - d_min),
                            SCATTERING_TEXTURE_MU_SIZE / 2);
    }

    float d = DistanceToTopAtmosphereBoundary(atmosphere.bottom_radius, mu_s);
    float d_min = atmosphere.top_radius - atmosphere.bottom_radius;
    float d_max = H;
    float a = (d - d_min) / (d_max - d_min);
    float A = -2.f * atmosphere.mu_s_min * atmosphere.bottom_radius / (d_max - d_min);
    u_mu_s = GetTextureCoordFromUnitRange(
                max(1.f - a / A, 1e-20) / (1.f + a),
                SCATTERING_TEXTURE_MU_S_SIZE);

    u_nu = (1.f + nu) / 2.f;

    return float4(u_r, u_mu, u_mu_s, u_nu);
}

struct SingleScatterTex
{
    float3 rayleigh;
    float3 mie;
    float3 single_scatter;
};

SingleScatterTex ComputeSingleScatteringTexture(QuadVertexOut In) : SV_Target
{
    SingleScatterTex res;
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
    for (int i = 0; i < SAMPLE_COUNT;++i)
    {
        float d = dx * i;
        float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
        float mu_d = (d + r * mu) / r_d;
        float mu_s_d = (r * mu_s + d * nu) / r_d;
        float altitude_d = r_d - atmosphere.bottom_radius;

        float transmittance = GetTransmittance(r, mu,r_d,mu_d,ray_r_mu_intersects_ground) *
                                GetTransmittanceToSun(r_d, mu_s_d);
        float weight = (i == 0 || i == SAMPLE_COUNT) ? 0.5f : 1.0f;
        
        rayleigh += weight * transmittance * GetLayerDesity(atmosphere.rayleigh_density, altitude_d);
        mie += weight * transmittance * GetLayerDesity(atmosphere.mie_density, altitude_d);
    }
    rayleigh *= dx * atmosphere.solar_irradiance * atmosphere.rayleigh_scattering;
    mie *= dx * atmosphere.solar_irradiance * atmosphere.mie_scattering;
    
    res.rayleigh = rayleigh;
    res.mie = mie;
    res.single_scatter = rayleigh + mie;
    return res;
}

technique11 ComputeSingleScaterTex3DTech
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