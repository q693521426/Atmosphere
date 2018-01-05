const int TRANSMITTANCE_TEXTURE_WIDTH = 256;    //mu
const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;    //r

const int SCATTERING_TEXTURE_R_SIZE = 32;
const int SCATTERING_TEXTURE_MU_SIZE = 128;
const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
const int SCATTERING_TEXTURE_NU_SIZE = 8;

const int SCATTERING_TEXTURE_WIDTH = 
            SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE;
const int SCATTERING_TEXTURE_HEIGHT = SCATTERING_TEXTURE_MU_SIZE;
const int SCATTERING_TEXTURE_DEPTH = SCATTERING_TEXTURE_R_SIZE;

const int IRRADIANCE_TEXTURE_WIDTH = 64;
const int IRRADIANCE_TEXTURE_HEIGHT = 16;

// The conversion factor between watts and lumens.
const double MAX_LUMINOUS_EFFICACY = 683.0;

SamplerState samLinearClamp
{
    Filter = MIN_MAG_MIP_LINEAR;
    AddressU = Clamp;
    AddressV = Clamp;
};

struct DensityProfileLayer
{
	float width;
	float exp_term;
	float exp_scale;
	float linear_term;
	float const_term;
};

struct DensityProfile
{
	DensityProfileLayer layer[2];
};

struct AtmosphereParameters
{
	float bottom_radius;
	float top_radius;

	DensityProfile rayleigh_density;
	float3 rayleigh_scattering;

	DensityProfile mie_density;
	float3 mie_scattering;


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

cbuffer cbAtmosphereParameters
{
	AtmosphereParameters atmosphere;
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

float GetProfileDesity(DensityProfile profile, float altitude)
{
	return altitude < profile.layer[0].width ?
		GetLayerDesity(profile.layer[0], altitude) : GetLayerDesity(profile.layer[1], altitude);
}

float DistanceToTopAtmosphereBoundary(float r, float mu)
{
	float discriminant = r * r * (mu * mu - 1) + atmosphere.top_radius * atmosphere.top_radius;
	return max(-r * mu + sqrt(max(discriminant, 0.f)), 0.f);
}

float DistanceToBottomAtmosphereBoundary(float r, float mu)
{
	float discriminant = r * r * (mu * mu - 1) + atmosphere.bottom_radius * atmosphere.bottom_radius;
	return max(-r * mu - sqrt(max(discriminant, 0.f)), 0.f);
}

float ComputeOpticalLengthToTopAtmosphereBoundary(DensityProfile profile,float r, float mu)
{
	const int SAMPLE_COUNT = 500;
	float dx = DistanceToTopAtmosphereBoundary(r, mu);

	float result = 0.f;
	for (int i = 0; i < SAMPLE_COUNT; i++)
	{
		float d = dx * i;
		float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
		// float mu_d = (d + r * mu) / r_d;
		float desity = GetProfileDesity(profile, r_d - atmosphere.bottom_radius);

		float weight = (i == 0 || i == SAMPLE_COUNT - 1) ? 0.5 : 1.0;

		result += weight * dx * desity;
	}
	return result;
}

float3 ComputeTransmittanceToTopAtmosphereBoundary(float r, float mu)
{
	//return exp(-(
	//	atmosphere.rayleigh_scattering *
	//	ComputeOpticalLengthToTopAtmosphereBoundary(
	//		atmosphere, atmosphere.rayleigh_density, r, mu) +
	//	atmosphere.mie_extinction *
	//	ComputeOpticalLengthToTopAtmosphereBoundary(
	//		atmosphere, atmosphere.mie_density, r, mu) +
	//	atmosphere.absorption_extinction *
	//	ComputeOpticalLengthToTopAtmosphereBoundary(
	//		atmosphere, atmosphere.absorption_density, r, mu)));
	return exp(-(
		atmosphere.rayleigh_scattering * 
				ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere.rayleigh_density, r, mu) +
		atmosphere.mie_scattering * 
				ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere.mie_density, r, mu)));
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
    float rho = sqrt(max(0.f, r * r - atmosphere.bottom_radius * atmosphere.bottom_radius));

    float r_x = rho / H;
    float d = DistanceToTopAtmosphereBoundary(r, mu);
    float d_min = atmosphere.top_radius - r;
    float d_max = rho + H;
    float mu_x = (d - d_min) / (d_max - d_min);
    
    return float2(GetTextureCoordFromUnitRange(mu_x, TRANSMITTANCE_TEXTURE_WIDTH),
                    GetTextureCoordFromUnitRange(r_x, TRANSMITTANCE_TEXTURE_HEIGHT));

}

float3 ComputeTransmittanceToTopAtmosphereBoundaryTexture(QuadVertexOut In)
{
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float2 RMu = GetRMuFromTransmittanceUV(f2UV);

    return ComputeTransmittanceToTopAtmosphereBoundary(RMu.x, RMu.y);
}

float3 GetTransmittanceToTopAtmosphereBoundary(float r, float mu)
{
    float2 uv = GetTransmittanceUVFromRMu(r, mu);
    return g_tex2DTransmittanceLUT.Sample(samLinearClamp,uv);
}

float3 GetTransmittance(float r, float mu, float d, bool ray_r_mu_intersects_ground)
{
    float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
    float mu_d = (d + r * mu) / r_d;
    
    if(ray_r_mu_intersects_ground)
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(r_d, -mu_d) /
                GetTransmittanceToTopAtmosphereBoundary(r, mu),
                float3(1.f));
    }
    else
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(r, mu) /
                GetTransmittanceToTopAtmosphereBoundary(r_d, mu_d),
                float3(1.f));
    }
}


technique11 PreComputeTransimittanceTextureTech
{
    pass P0
    {
        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeTransmittanceToTopAtmosphereBoundaryTexture()));
    }
}