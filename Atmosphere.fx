#include "Common.fx"

float RayleighPhaseFunction(float w)
{
    float k = 3.0 / (16.0 * PI);
    return k * (1.0 + w * w);
}

float MiePhaseFunction(float g, float w)
{
    float g2 = g * g;
    float k = 3.0 / (8.0 * PI) * (1.0 - g2) / (2.0 + g2);
    return k * (1.0 + w * w) / pow(max(1.0 + g2 - 2.0 * g * w, 1e-20), 1.5);
}

float GetLayerDensity(DensityProfileLayer layer, float altitude)
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

float3 GetTransmittanceIntegralAnalytic(float r, float mu, float d)
{
    float2 f2A = sqrt((0.5 / float2(8,1.2)) * r);
    float4 f4A01 = f2A.xxyy * float2(mu, mu + d / r).xyxy;
    float4 f4A01s = sign(f4A01);
    float4 f4A01sq = f4A01 * f4A01;
    
    float2 f2X;
    f2X.x = f4A01s.y > f4A01s.x ? exp(f4A01sq.x) : 0.0;
    f2X.y = f4A01s.w > f4A01s.z ? exp(f4A01sq.z) : 0.0;
    
    float4 f4Y = f4A01s / (2.3193 * abs(f4A01) + sqrt(1.52 * f4A01sq + 4.0)) * float3(1.0, exp(-d / float2(8, 1.2) * (d / (2.0 * r) + mu))).xyxz;

    float2 optical_length = sqrt((6.2831 * float2(8, 1.2)) * r) * exp((atmosphere.bottom_radius - r) / float2(8, 1.2)) * (f2X + float2(dot(f4Y.xy, float2(1.0, -1.0)), dot(f4Y.zw, float2(1.0, -1.0))));

    return exp(-(atmosphere.rayleigh_scattering * optical_length.x +
                 atmosphere.mie_extinction * optical_length.y));
}

float3 ComputeOpticalLengthToTopAtmosphereBoundary(float r, float mu)
{
    const int SAMPLE_COUNT = 128;
    float dx = DistanceToTopAtmosphereBoundary(r, mu) / float(SAMPLE_COUNT);
    float3 result = 0.f;
    float3 pre = 0.f;
    for (int i = 0; i <= SAMPLE_COUNT; i++)
    {
        float d = dx * i;
        float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
		// float mu_d = (d + r * mu) / r_d;
        
        float altitude = r_d - atmosphere.bottom_radius;
        float rayleigh = GetLayerDensity(atmosphere.rayleigh_density, altitude);
        float mie = GetLayerDensity(atmosphere.mie_density, altitude);
        float ozone = altitude < atmosphere.ozone_width ?
                                GetLayerDensity(atmosphere.ozone_density[0], altitude) :
                                  GetLayerDensity(atmosphere.ozone_density[1], altitude);

        //float weight = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;

        //result += weight * float3(rayleigh, mie, ozone)*dx;
        result += (float3(rayleigh, mie, ozone) + pre) * dx / 2;
        pre = float3(rayleigh, mie, ozone);
    }
    return result;
}

float3 ComputeTransmittanceToTopAtmosphereBoundary(float r, float mu)
{
    float3 optical_length = ComputeOpticalLengthToTopAtmosphereBoundary(r, mu);
#if USE_OZONE_DENSITY
    return exp(-(atmosphere.rayleigh_scattering * optical_length.x +
                 atmosphere.mie_extinction * optical_length.y +
                 atmosphere.absorption_extinction * optical_length.z));
#else
   return exp(-(atmosphere.rayleigh_scattering * optical_length.x +
                atmosphere.mie_extinction * optical_length.y));
#endif
}

float2 GetRMuFromTransmittanceUV(float2 uv)
{
    float mu_x = GetUnitRangeFromTextureCoord(uv.x, TRANSMITTANCE_TEXTURE_WIDTH);
    float r_x = GetUnitRangeFromTextureCoord(uv.y, TRANSMITTANCE_TEXTURE_HEIGHT);

    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = H * r_x;
    float r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
    float d_min = atmosphere.top_radius - r;
    float d_max = rho + H;
    float d = d_min + mu_x * (d_max - d_min);
    float mu = d == 0 ? 1.0 : (H * H - rho * rho - d * d) / (2 * r * d);
    mu = CosClamp(mu);
    //r = clamp(r, atmosphere.bottom_radius + 10, atmosphere.top_radius - 10);
    return float2(r, mu);

}

float2 GetTransmittanceUVFromRMu(float r, float mu)
{
    //r = clamp(r, atmosphere.bottom_radius + 10, atmosphere.top_radius - 10);
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

float3 ComputeOpticalToTopAtmosphereBoundaryTexture(QuadVertexOut In) : SV_Target
{
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float2 RMu = GetRMuFromTransmittanceUV(f2UV);
    
    return ComputeOpticalLengthToTopAtmosphereBoundary(RMu.x, RMu.y);
}

float3 ComputeTransmittanceToTopAtmosphereBoundaryTexture(QuadVertexOut In) : SV_Target
{
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float2 RMu = GetRMuFromTransmittanceUV(f2UV);
    
    return ComputeTransmittanceToTopAtmosphereBoundary(RMu.x, RMu.y);
}

technique11 ComputeOpticalLengthTex2DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeOpticalToTopAtmosphereBoundaryTexture()));
    }
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
//#if !USE_OPTICAL_LUT
//    return g_tex2DTransmittanceLUT.SampleLevel(samLinearClamp, uv,0);
//#endif
    float3 optical_length = g_tex2DOpticalLengthLUT.SampleLevel(samLinearClamp, uv, 0);
#if USE_OZONE_DENSITY
    return exp(-(atmosphere.rayleigh_scattering * optical_length.x +
                 atmosphere.mie_extinction * optical_length.y +
                 atmosphere.absorption_extinction * optical_length.z));
#else
   return exp(-(atmosphere.rayleigh_scattering * optical_length.x +
                atmosphere.mie_extinction * optical_length.y));
#endif

}

float3 GetTransmittance(float r, float mu, float r_d, float mu_d, bool ray_r_mu_intersects_ground)
{
    if (ray_r_mu_intersects_ground)
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(r_d, -mu_d) /
                GetTransmittanceToTopAtmosphereBoundary(r, -mu),
                float3(1.f, 1.f, 1.f));
        //return float3(1.f, 1.f, 1.f);
    }
    else
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(r, mu) /
                GetTransmittanceToTopAtmosphereBoundary(r_d, mu_d),
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
    float u_mu_s = GetUnitRangeFromTextureCoord(uvwq.z, SCATTERING_TEXTURE_MU_S_SIZE);
    float u_nu = GetUnitRangeFromTextureCoord(uvwq.w, SCATTERING_TEXTURE_NU_SIZE);
#if USE_LUT_PARAMETERIZATION
    float h = pow(u_r, 2) * (atmosphere.top_radius - atmosphere.bottom_radius);
    r = h + atmosphere.bottom_radius;
    float mu_horizon = -sqrt(h * (2 * atmosphere.bottom_radius + h)) / (h + atmosphere.bottom_radius);
    if (uvwq.y > 0.5)
    {
        float u_mu = GetUnitRangeFromTextureCoord((uvwq.y - 0.5) * 2, SCATTERING_TEXTURE_MU_SIZE/2);
        mu = pow(u_mu, 1 / 0.2) * (1 - mu_horizon) + mu_horizon;
        ray_r_mu_intersects_ground = false;
    }
    else
    {
        float u_mu = GetUnitRangeFromTextureCoord(uvwq.y * 2, SCATTERING_TEXTURE_MU_SIZE/2);
        mu = mu_horizon - pow(u_mu, 1 / 0.2) * (mu_horizon - (-1));
        ray_r_mu_intersects_ground = true;
    }
    mu_s = tan((2.0 * u_mu_s - 1.0 + 0.26) * 1.1) / tan(1.26 * 1.1);

    nu = sign(u_nu - 0.5) * pow(abs((u_nu - 0.5) * 2), 1 / atmosphere.nu_power) / 2 + 0.5;
    nu = cos(nu * PI);
    //nu = clamp((u_nu - 0.5f) * 2.f, -1.f, 1.f);
#else
    //float u_nu = uvwq.w;

    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
                   atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = u_r * H;

    r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
    
    if (uvwq.y < 0.5) // [-1,mu_horizon]
    {
        float d_min = r - atmosphere.bottom_radius;
        float d_max = rho;
        float d = GetUnitRangeFromTextureCoord(1 - 2 * uvwq.y, SCATTERING_TEXTURE_MU_SIZE / 2) * (d_max - d_min) + d_min;
        mu = d == d_min ? -1.f :
                -(d * d + rho * rho) / (2.f * r * d);
        ray_r_mu_intersects_ground = true;
    }
    else
    {
        float d_min = atmosphere.top_radius - r;
        float d_max = rho + H;
        float d = GetUnitRangeFromTextureCoord(2 * uvwq.y -1 , SCATTERING_TEXTURE_MU_SIZE / 2) * (d_max - d_min) + d_min;
        mu = d == 0.f ? 1.f :
                -(d * d - H * H + rho * rho) / (2.f * r * d);
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
#endif
}

float4 GetUVWQFromrRMuMuSNu(in float r, in float mu,
                                in float mu_s, in float nu)
{
#if USE_LUT_PARAMETERIZATION
    float h = max(r - atmosphere.bottom_radius,0);
    float u_h = GetTextureCoordFromUnitRange(pow(h / (atmosphere.top_radius - atmosphere.bottom_radius), 0.5), SCATTERING_TEXTURE_R_SIZE);
    float mu_horizon = -sqrt(h * (2 * atmosphere.bottom_radius + h)) / (h + atmosphere.bottom_radius);
    float u_mu;
    if(mu>mu_horizon)
    {
        u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(pow((mu - mu_horizon) / (1 - mu_horizon), 0.2), 
                                                        SCATTERING_TEXTURE_MU_SIZE / 2);

    }
    else
    {
        u_mu = 0.5 * GetTextureCoordFromUnitRange(pow((mu_horizon - mu) / (mu_horizon - (-1)), 0.2),
                                                        SCATTERING_TEXTURE_MU_SIZE / 2);
    }
    float u_mu_s = GetTextureCoordFromUnitRange((atan(max(mu_s, atmosphere.mu_s_min) * tan(1.26 * 1.1)) / 1.1 + (1.0 - 0.26)) * 0.5,
                                                        SCATTERING_TEXTURE_MU_S_SIZE);
    
    float u_nu = acos(nu) / PI;
    u_nu = sign(u_nu - 0.5) * pow(abs((u_nu - 0.5) / 0.5), atmosphere.nu_power) / 2 + 0.5;
    u_nu = GetTextureCoordFromUnitRange(u_nu, SCATTERING_TEXTURE_NU_SIZE);
    //float u_nu = GetTextureCoordFromUnitRange((1.f + nu) / 2.f, SCATTERING_TEXTURE_NU_SIZE);
    
    //nu = sign(u_nu - 0.5) * pow(abs((u_nu - 0.5) * 2), 1 / misc.nu_power) / 2 + 0.5;
    //nu = cos(nu * PI);

    return float4(u_h, u_mu, u_mu_s, u_nu);
#else
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
                            d_max == d_min ? 0.f : (d - d_min) / (d_max - d_min),
                            SCATTERING_TEXTURE_MU_SIZE / 2);
    }
    float a = 1 - exp(-2.8 * mu_s - 0.8);
    float A = 1 - exp(-3.6);
    u_mu_s = GetTextureCoordFromUnitRange(max(a / A, 1e-20), SCATTERING_TEXTURE_MU_S_SIZE);

    u_nu = GetTextureCoordFromUnitRange((1.f + nu) / 2.f, SCATTERING_TEXTURE_NU_SIZE);
    //u_nu = (1.f + nu) / 2.f;

    return float4(u_r, u_mu, u_mu_s, u_nu);
#endif
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

    float3 rayleigh, mie, pre_rayleigh, pre_mie;
    float3 total_density, pre_density;
    float3 transmittance_cam = 1;
    for (int i = 0; i <= SAMPLE_COUNT;i++)
    {
        float d = dx * i;
        float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
        //r_d = clamp(r_d,atmosphere.bottom_radius, atmosphere.top_radius);
        float mu_d = CosClamp((d + r * mu) / r_d);
        float mu_s_d = CosClamp((r * mu_s + d * nu) / r_d);
        float altitude_d = max(r_d - atmosphere.bottom_radius, 1e-20);

        float rayleigh_density = GetLayerDensity(atmosphere.rayleigh_density, altitude_d);
        float mie_density = GetLayerDensity(atmosphere.mie_density, altitude_d);
        float ozone_density = altitude_d < atmosphere.ozone_width ?
                                GetLayerDensity(atmosphere.ozone_density[0], altitude_d) :
                                  GetLayerDensity(atmosphere.ozone_density[1], altitude_d);
        float3 density = float3(rayleigh_density, mie_density, ozone_density);
        total_density += (pre_density + density) * dx / 2;
        pre_density = density;
#if USE_OPTICAL_LUT       
        transmittance_cam = GetTransmittance(r, mu, r_d, mu_d, ray_r_mu_intersects_ground);
#elif USE_TRANSMITTANCE_ANALYTIC
        if(i==0)
            transmittance_cam = 1;
        else
        {
            float d = dx * (i - 1);
            float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
            //r_d = clamp(r_d,atmosphere.bottom_radius, atmosphere.top_radius);
            float mu_d = CosClamp((d + r * mu) / r_d);
            float3 cur_transmittance = GetTransmittanceIntegralAnalytic(r_d, mu_d, dx);
            transmittance_cam *= cur_transmittance;
        }
#elif USE_OZONE_DENSITY
        transmittance_cam = exp(-(atmosphere.rayleigh_scattering * total_density.x +
                                    atmosphere.mie_extinction * total_density.y +
                                        atmosphere.absorption_extinction * total_density.z));
#else
        transmittance_cam = exp(-(atmosphere.rayleigh_scattering * total_density.x +
                                    atmosphere.mie_extinction * total_density.y));
#endif
        float3 transmittance_sun = GetTransmittanceToTopAtmosphereBoundary(r_d, mu_s_d);
        
        float3 transmittance = transmittance_cam * transmittance_sun;
        
        float3 cur_rayleigh = transmittance * rayleigh_density;
        float3 cur_mie = transmittance * mie_density;

        rayleigh += (pre_rayleigh + cur_rayleigh) * dx / 2;
        pre_rayleigh = cur_rayleigh;

        mie += (pre_mie + cur_mie) * dx / 2;
        pre_mie = cur_mie;
    }
    rayleigh *= (atmosphere.solar_irradiance * atmosphere.rayleigh_scattering);
    mie *= (atmosphere.solar_irradiance * atmosphere.mie_scattering);

    res.single_scatter = float4(rayleigh * RayleighPhaseFunction(nu) + mie * MiePhaseFunction(atmosphere.mie_g, nu), 1);
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

    float2 f2SliceMinMaxZ = float2(fQ0Slice, fQ0Slice + 1) / SCATTERING_TEXTURE_NU_SIZE + float2(0.5, -0.5) / (SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE);
    f3UVW0.z = clamp(f3UVW0.z, f2SliceMinMaxZ.x, f2SliceMinMaxZ.y);

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
    const float dtheta = PI / SAMPLE_COUNT;
    const float dphi = PI / SAMPLE_COUNT;

    float3 zenith = float3(0.f, 1.f, 0.f);
    float3 omega = float3(sqrt(1.f - mu * mu), mu, 0.f);
    float sun_dir_x = omega.x == 0.f ? 0.f : (nu - mu * mu_s) / omega.x;
    float sun_dir_z = sqrt(1 - sun_dir_x * sun_dir_x - mu_s * mu_s);
    float3 omega_s = float3(sun_dir_x, mu_s, sun_dir_z);

    float3 rayleigh = 0.f;
    float3 mie = 0.f;

    float altitude = r - atmosphere.bottom_radius;
    float rayleigh_density = GetLayerDensity(atmosphere.rayleigh_density, altitude);
    float mie_density = GetLayerDensity(atmosphere.mie_density, altitude);

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
            transmittance_to_ground = float3(1, 1, 1);
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
            incident_radiance += transmittance_to_ground * ground_albedo * (1.0 / PI) * ground_irradiance;

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
    
    float3 transmittance = 1;
    float3 pre_density, total_density;
    float3 transmittance_cam;
    for (int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        float d = dx * i;
        float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
        float mu_d = CosClamp((r * mu + d) / r_d);
        float mu_s_d = CosClamp((r * mu_s + d * nu) / r_d);

        float altitude_d = r_d - atmosphere.bottom_radius;

        float3 scatter_density = ComputeScatteringDensity(r_d, mu_d, mu_s_d, nu, misc.scatter_order);
        
        float rayleigh_density = GetLayerDensity(atmosphere.rayleigh_density, altitude_d);
        float mie_density = GetLayerDensity(atmosphere.mie_density, altitude_d);
        float ozone_density = altitude_d < atmosphere.ozone_width ?
                                GetLayerDensity(atmosphere.ozone_density[0], altitude_d) :
                                  GetLayerDensity(atmosphere.ozone_density[1], altitude_d);
        float3 density = float3(rayleigh_density, mie_density, ozone_density);
        total_density += (pre_density + density) * dx / 2;
        pre_density = density;
#if USE_OPTICAL_LUT       
        transmittance_cam = GetTransmittance(r, mu, r_d, mu_d, ray_r_mu_intersects_ground);
#elif USE_TRANSMITTANCE_ANALYTIC
        if (i == 0)
            transmittance_cam = 1;
        else
        {
            float d = dx * (i - 1);
            float r_d = sqrt(r * r + d * d + 2 * r * d * mu);
            //r_d = clamp(r_d,atmosphere.bottom_radius, atmosphere.top_radius);
            float mu_d = CosClamp((d + r * mu) / r_d);
            float3 cur_transmittance = GetTransmittanceIntegralAnalytic(r_d, mu_d, dx);
            transmittance_cam *= cur_transmittance;
        }
#elif USE_OZONE_DENSITY
        transmittance_cam = exp(-(atmosphere.rayleigh_scattering * total_density.x +
                                    atmosphere.mie_extinction * total_density.y +
                                        atmosphere.absorption_extinction * total_density.z));
#else
        transmittance_cam = exp(-(atmosphere.rayleigh_scattering * total_density.x +
                                    atmosphere.mie_extinction * total_density.y));
#endif
        
        float3 res = transmittance_cam * scatter_density;
        rayleigh_mie += res;
    }
    float w = RayleighPhaseFunction(nu);
    multi_scatter.multi_scatter = float4(rayleigh_mie, 1);
    multi_scatter.scatter_combined = float4(w == 0 ?  float3(0,0,0): rayleigh_mie / w,0) +
                                     GetScatteringCombined(r, mu, mu_s, nu,1);
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
    float dphi = PI / SAMPLE_COUNT;
    float dtheta = PI / SAMPLE_COUNT;

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
            float3 scatter = GetScattering(r, omega_i.y, mu_s, nu, misc.scatter_order);
            result += scatter * dSolidAngle * omega_i.y;
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
        //float alpha_s = atmosphere.sun_angular_radius;
        ////Approximate average of the cosine factor mu_s over the visible fraction of
        ////the Sun disc.
        //float average_cosine_factor =
        //        mu_s < -alpha_s ? 0.0 : (mu_s > alpha_s ? mu_s :
        //            (mu_s + alpha_s) * (mu_s + alpha_s) / (4.0 * alpha_s));
        //irradiance_from_sun = atmosphere.solar_irradiance * transmmitance_to_sun * 
        //                          average_cosine_factor * max(nu, 0.f);
        irradiance_from_sun = g_tex2DDirectIrradianceLUT.Sample(samLinearClamp, irradiance_uv) * max(nu, 0.f);
    }
    if (nu > cos(atmosphere.sun_angular_radius))
        irradiance_from_sun += transmmitance_to_sun * atmosphere.solar_irradiance /
                                (PI * atmosphere.sun_angular_radius * atmosphere.sun_angular_radius);
    else
        irradiance_from_sun = float3(0, 0, 0);
    return irradiance_from_sun;
}

float3 GetSkyMultiScatter(float r, float mu, float mu_s, float nu)
{
#if USE_SCATTER_COMBINED
    float4 scatter_combined = GetScatteringCombined(r, mu, mu_s, nu, misc.scatter_order);
   
    float3 mie_scatter = scatter_combined.r == 0 ? float3(0,0,0):
        scatter_combined.rgb * scatter_combined.a / scatter_combined.r *
	    (atmosphere.rayleigh_scattering.r / atmosphere.mie_scattering.r) *
	    (atmosphere.mie_scattering / atmosphere.rayleigh_scattering);
    
    return (mie_scatter * MiePhaseFunction(atmosphere.mie_g, nu) + scatter_combined.rgb * RayleighPhaseFunction(nu));
#else
    if(misc.scatter_order==1)
        return GetScattering(r,mu,mu_s,nu,1);
    else
        return GetScattering(r,mu,mu_s,nu,2);
#endif
}

float3 GetSkyMultiScatterToGround(float r, float mu, float mu_s,float r_d, float mu_d, float mu_s_d, 
                                    float nu, float shadow_length, float d, float3 transmittance)
{
    // S[L]x - T(x,xs)S[L]xs=x0-lv
    float3 scatter = GetSkyMultiScatter(r, mu, mu_s, nu);
    float3 scatter_d;

    float3 shadow_transmittance;
    if (0.f == shadow_length)
    {
        scatter_d = GetSkyMultiScatter(r_d, mu_d, mu_s_d, nu);
        shadow_transmittance = transmittance;
    }
    else
    {
        d = max(d - shadow_length, 0.f);
        r_d = sqrt(r * r + d * d + 2 * r * d * mu);
        mu_d = (r * mu + d) / r_d;
        mu_s_d = (r * mu_s + d * nu) / r_d;
        scatter_d = GetSkyMultiScatter(r_d, mu_d, mu_s_d, nu);
        shadow_transmittance = GetTransmittance(r, mu,r_d,mu_d,true);
        //shadow_transmittance = GetTransmittanceIntegralAnalytic(r, mu, d);
    }
    return scatter - shadow_transmittance * scatter_d;
}

float3 GetSkyMultiScatterToAtmosphere(float r, float mu, float mu_s, float nu, float shadow_length)
{
    // T(x,xs)S[L]xs=x+lv
    float3 scatter_d;
    float3 shadow_transmittance;
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
        //shadow_transmittance = GetTransmittanceIntegralAnalytic(r, mu, d);
        return shadow_transmittance * scatter_d;
    }
}

float3 HDR(float3 L)
{
    L = L * 10.f;
    L.r = L.r < 1.413 ? pow(L.r * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.r);
    L.g = L.g < 1.413 ? pow(L.g * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.g);
    L.b = L.b < 1.413 ? pow(L.b * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.b);
    return L;
}

float4 DrawGroundAndSky(QuadVertexOut In) : SV_Target
{
    float4 f4RayEndWorldPos = mul(float4(In.m_f2PosPS, 1.0f, 1.0f), camera.InvViewProj);
    float3 v = f4RayEndWorldPos.xyz / f4RayEndWorldPos.w - camera.f3CameraPos;
    float v_length = length(v);
    v = v / v_length;
    
    float fDepth = g_tex2DSpaceDepth.Load(uint3(In.m_f4Pos.xy, 0));
    if(fDepth != 1)
    {
        return float4(g_tex2DColorBuffer.Load(uint3(In.m_f4Pos.xy, 0)), 1);
    }

    //float fragment_angular_size =
    //  length(ddx(In.m_f4Ray) + ddy(In.m_f4Ray)) / length(In.m_f4Ray);
    float3 f3EarthCenter = float3(0, -atmosphere.bottom_radius, 0);
    float3 camera_pos = camera.f3CameraPos - f3EarthCenter;
    float3 sun_dir = light.f3LightDir;
    float r = length(camera_pos);
    float mu = dot(camera_pos, v) / r; 

    float nu = dot(v, sun_dir);
    float mu_s = dot(camera_pos, sun_dir) / r;

    float3 ground_radiance = 0.f;
    float3 sky_radiance = 0.f;

    float intersect_ground_discriminant = r * r * (mu * mu - 1) + atmosphere.bottom_radius * atmosphere.bottom_radius;
    bool ray_r_mu_intersects_ground = intersect_ground_discriminant >= 0 && mu < 0;
    float3 irradiance_from_sun = GetIrradianceDirectFromSun(r, mu_s, nu); //L0
    //return float4(irradiance_from_sun, 1);
    float intersect_atmosphere_discriminant = r * r * (mu * mu - 1) + atmosphere.top_radius * atmosphere.top_radius;
    if (r > atmosphere.top_radius &&
        (intersect_atmosphere_discriminant < 0 || intersect_atmosphere_discriminant >= 0 && mu > 0))
    {
        return float4(0, 0, 0, 1);
    }

    float distance_to_atmosphere_near = -r * mu - sqrt(intersect_atmosphere_discriminant);
    if (distance_to_atmosphere_near > 0.f)
    {
        camera_pos += distance_to_atmosphere_near * v;
        r = atmosphere.top_radius;
        mu = dot(camera_pos, v) / r;
        mu_s = dot(camera_pos, sun_dir) / r;
        intersect_ground_discriminant = r * r * (mu * mu - 1) + atmosphere.bottom_radius * atmosphere.bottom_radius;
    }
    

    if (ray_r_mu_intersects_ground)
    {
        float distance_to_ground = -r * mu - sqrt(intersect_ground_discriminant);
        float3 pos_d = camera_pos + distance_to_ground * v;
        float r_d = length(pos_d);
        float3 normal_d = pos_d / r_d;
        float mu_d = CosClamp((distance_to_ground + r * mu) / r_d);
        float mu_s_d = CosClamp((distance_to_ground * nu + r * mu_s) / r_d);
        //if (mu_s_d <= dot(pos_d, sun_dir) / r_d + 0.01 && mu_s_d >= dot(pos_d, sun_dir) / r_d - 0.01)
        //    return float4(1, 1, 1, 1);
        
        float2 coords = float2(normal_d.x == 0 ? 0 : atan2(normal_d.z, normal_d.x), acos(normal_d.y)) * float2(0.5, 1.0) / PI + float2(0.5, 0.0);
        float4 reflectance = g_tex2DEarthGround.SampleLevel(samLinearClamp, coords, 0); //* float4(0.2, 0.2, 0.2, 1.0);
       
        //float3 transmittance = GetTransmittanceIntegralAnalytic(r, mu, distance_to_ground);
        float3 transmittance = GetTransmittance(r, mu, r_d, mu_d, true);

        float2 irradiance_uv = GetIrradianceUVFromRMuS(r_d, mu_s_d);

        // R[L0] + R[L*]
        float3 sunColor = g_tex2DDirectIrradianceLUT.Sample(samLinearClamp, irradiance_uv) * max(mu_s_d, 0);
        float3 skyColor = g_tex2DIndirectIrradianceLUT.Sample(samLinearClamp, irradiance_uv) * (1.0 + dot(normal_d, camera_pos) / r) * 0.5;
        //float3 skyColor = g_tex2DIndirectIrradianceLUT.Sample(samLinearClamp, irradiance_uv);
        
        float shadow_length = 0;
        // S[L]x - T(x,xs)S[L]xs=x0-lv
        float3 scatter = GetSkyMultiScatterToGround(r, mu, mu_s, r_d, mu_d, mu_s_d, nu, shadow_length, distance_to_ground, transmittance);
        sky_radiance = scatter;

        ground_radiance = reflectance.rgb * atmosphere.ground_albedo * (1.0 / PI) * (sunColor + skyColor) * transmittance;
        
        return float4(ToneMap(ground_radiance + sky_radiance), 1.0);
    }
    else
    {
        float distance_to_atmosphere = -r * mu + sqrt(intersect_atmosphere_discriminant);
        float3 pos_d = camera_pos + distance_to_atmosphere * v;
        float r_d = length(pos_d);
        float mu_d = CosClamp((distance_to_atmosphere + r * mu) / r_d);
        float mu_s_d = CosClamp((distance_to_atmosphere * nu + r * mu_s) / r_d);
    
        //float3 transmittance = GetTransmittanceIntegralAnalytic(r, mu, distance_to_atmosphere);
        float3 transmittance = GetTransmittance(r, mu, r_d, mu_d, false);
        
        // T(x,xs)S[L]xs=x+lv
        float3 scatter = GetSkyMultiScatterToAtmosphere(r, mu, mu_s, nu, 0);
        sky_radiance = scatter;
        
        return float4(ToneMap(sky_radiance + irradiance_from_sun), 1);
    }
    float3 white_point = 1;
    float3 radiance = irradiance_from_sun + ground_radiance + sky_radiance;
    return float4(ToneMap(radiance), 1);
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

float ComputeSpaceLinearDepthTex2D(QuadVertexOut In) : SV_Target
{
    float depth = g_tex2DSpaceDepth.Load(int3(In.m_f4Pos.xy, 0));
    return camera.Proj[3][2] / (depth - camera.Proj[2][2]);
}

technique11 ComputeSpaceLinearDepthTex2DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeSpaceLinearDepthTex2D()));
    }
}

float4 INVALID_EPIPOLAR_LINE = float4(-1000, -1000, -100, -100);
float4 ComputeSliceEndTex2D(QuadVertexOut In) : SV_Target
{
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float fEpipolarSlice = GetUnitRangeFromTextureCoord(f2UV.x, EPIPOLAR_SLICE_NUM);
    //float fEpipolarSlice = f2UV.x - 0.5 / (float)EPIPOLAR_SLICE_NUM;

    float4 f4ScreenPixelCoord = float4(-1, -1, 1, 1) + float4(1, 1, -1, -1) / float2(SCREEN_WIDTH, SCREEN_HEIGHT).xyxy;

    float fBoundary = fEpipolarSlice * 4;
    uint uiBoundary = floor(fBoundary);
    uiBoundary = uiBoundary & 3;

    float fPosOnBoundary = frac(fBoundary);
    uint4 ui4BoundaryFlags = uint4(uiBoundary.xxxx == uint4(0, 1, 2, 3));

    float4 f4LightXYDistToBoundary = (light.f4LightScreenPos.xyxy - f4ScreenPixelCoord) * float4(1, 1, -1, -1);
    uint4 ui4IsInvalidBoundary = uint4(f4LightXYDistToBoundary <= 0);
    if (dot(ui4BoundaryFlags, ui4IsInvalidBoundary))
        return INVALID_EPIPOLAR_LINE;

    float4 f4BoundaryX = float4(0, fPosOnBoundary, 1, 1 - fPosOnBoundary);
    float4 f4BoundaryY = float4(1 - fPosOnBoundary, 0, fPosOnBoundary, 1);
    float2 f2SliceEndPos = float2(dot(f4BoundaryX, ui4BoundaryFlags), dot(f4BoundaryY, ui4BoundaryFlags));
    f2SliceEndPos = lerp(f4ScreenPixelCoord.xy, f4ScreenPixelCoord.zw, f2SliceEndPos);

    float2 f2SliceStartPos;
    if (all(f4LightXYDistToBoundary >= 0))
    {
        f2SliceStartPos = light.f4LightScreenPos.xy;
    }
    else
    {
        float2 f2RayDir = f2SliceEndPos - light.f4LightScreenPos.xy;
        float fRayLength = length(f2RayDir);
        f2RayDir /= fRayLength;
        
        bool4 b4IsCorrectIntersectionFlag = abs(f2RayDir.xyxy) > 1e-5;
        float4 f4DistToBoundary = (f4ScreenPixelCoord - light.f4LightScreenPos.xyxy) / (f2RayDir.xyxy + !b4IsCorrectIntersectionFlag);
        b4IsCorrectIntersectionFlag = b4IsCorrectIntersectionFlag && (f4DistToBoundary < (fRayLength * (1 - 1e-5)));
        f4DistToBoundary = b4IsCorrectIntersectionFlag * f4DistToBoundary +
                            !b4IsCorrectIntersectionFlag * float4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX);

        float2 f2FirstIntersecDist = max(f4DistToBoundary.xy, f4DistToBoundary.zw);
        float fFirstIntersecDist = max(f2FirstIntersecDist.x, f2FirstIntersecDist.y);

        f2SliceStartPos = light.f4LightScreenPos.xy + f2RayDir * fFirstIntersecDist;
    }

    if (IsValidScreenLocation(f2SliceStartPos))
    {
        // Compute length of the epipolar line in screen pixels:
        float fEpipolarSliceScreenLen = length((f2SliceEndPos - f2SliceStartPos) * float2(SCREEN_WIDTH, SCREEN_HEIGHT) / 2);
        // If epipolar line is too short, update epipolar line exit point to provide 1:1 texel to screen pixel correspondence:
        f2SliceEndPos = f2SliceStartPos + (f2SliceEndPos - f2SliceStartPos) * max((float) EPIPOLAR_SAMPLE_NUM / fEpipolarSliceScreenLen, 1);
    }

    return float4(f2SliceStartPos, f2SliceEndPos);
}

technique11 ComputeSliceEndTex2DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeSliceEndTex2D()));
    }
}

void ComputeEpipolarCoordTex2D(QuadVertexOut In, out float2 f2XY : SV_Target0, out float fCameraZ : SV_Target1) 
{
    float4 f4SliceEnd = g_tex2DSliceEnd.Load(uint3(In.m_f4Pos.y, 0,0));
    if (!IsValidScreenLocation(f4SliceEnd.xy))
    {
        discard;
    }
    //float fEpipolarSample = (In.m_f4Pos.x - 0.5) / ((float)EPIPOLAR_SAMPLE_NUM - 1);
    float fEpipolarSample = GetUnitRangeFromTextureCoord(ProjToUV(In.m_f2PosPS).x, EPIPOLAR_SAMPLE_NUM);
    f2XY = lerp(f4SliceEnd.xy, f4SliceEnd.zw, fEpipolarSample);
    if (!IsValidScreenLocation(f2XY))
    {
        discard;
    }
    fCameraZ = g_tex2DSpaceLinearDepth.SampleLevel(samLinearClamp, ProjToUV(f2XY), 0);
}

static const float3 m_f4ExtraLight = float3(5.f, 5.f, 5.f);

float4 ComputeUnshadowedSampleScatter(float2 f2ScreenXY, float fCamDepth, float fDepth)
{
    float4 f4RayEndInWorldSpace = mul(float4(f2ScreenXY, fDepth, 1), camera.InvViewProj);
    float3 f3ViewRay = f4RayEndInWorldSpace.xyz / f4RayEndInWorldSpace.w - camera.f3CameraPos;
    float fViewRayLen = length(f3ViewRay);
    f3ViewRay /= fViewRayLen;
    
    float3 f3EarthCenter = float3(0, -atmosphere.bottom_radius, 0);
    float3 f3CameraPos = camera.f3CameraPos - f3EarthCenter;

    float4 f4RMuMuSNu;
    f4RMuMuSNu.xyz = GetRMuMuS(f3CameraPos, f3ViewRay);
    f4RMuMuSNu.w = dot(f3ViewRay, light.f3LightDir);
    
    bool bIsNoScatter, bIsMarchToAtmosphere, bIsMarchToEarth, bIsIntersectEarth;
    float fDistToAtmosphereNear, fDistToAtmosphereFar, fDistToEarthNear;
    float fViewRayLenInWorldSpace = GetRayMarchLen(f4RMuMuSNu, fCamDepth, fViewRayLen,
                                            bIsNoScatter, bIsMarchToAtmosphere, bIsMarchToEarth, bIsIntersectEarth,
                                            fDistToAtmosphereNear, fDistToAtmosphereFar, fDistToEarthNear);
    
    if (fDistToAtmosphereNear > 0)
    {
        float3 f3DistToAtmosphereNear = fDistToAtmosphereNear * f3ViewRay;
        f3CameraPos += f3DistToAtmosphereNear;
        fViewRayLenInWorldSpace -= fDistToAtmosphereNear;
        f4RMuMuSNu.xyz = GetRMuMuS(f3CameraPos, f3ViewRay);
    }
    
    if (bIsMarchToAtmosphere || bIsMarchToEarth)
        return float4(m_f4ExtraLight * GetSkyMultiScatter(f4RMuMuSNu.x, f4RMuMuSNu.y, f4RMuMuSNu.z, f4RMuMuSNu.w), fViewRayLenInWorldSpace);
    else
    {
        float3 f3EndPos = f3CameraPos + fViewRayLenInWorldSpace * f3ViewRay;
        float3 f3EndRMuMuS = GetRMuMuS(f3EndPos, f3ViewRay);
        float3 f3TransmittanceToEnd = bIsIntersectEarth ? 1 : GetTransmittanceIntegralAnalytic(f4RMuMuSNu.x, f4RMuMuSNu.y, fViewRayLenInWorldSpace);
        float3 f3Inscatter = GetSkyMultiScatter(f4RMuMuSNu.x, f4RMuMuSNu.y, f4RMuMuSNu.z, f4RMuMuSNu.w);
        f3Inscatter -= f3TransmittanceToEnd * GetSkyMultiScatter(f3EndRMuMuS.x, f3EndRMuMuS.y, f3EndRMuMuS.z, f4RMuMuSNu.w);
        return float4(m_f4ExtraLight * f3Inscatter, fViewRayLenInWorldSpace);
    }
}

technique11 ComputeEpipolarCoordTex2DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest_IncrStencil, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeEpipolarCoordTex2D()));
    }
}

float4 DoUnshadowedRayMarch(QuadVertexOut In) : SV_Target
{
    uint uiSampleNum = In.m_f4Pos.x;
    uint uiSliceNum = In.m_f4Pos.y;
    float2 f2ScreenXY = g_tex2DEpipolarSample.Load(uint3(uiSampleNum, uiSliceNum, 0));
    float fCamDepth = g_tex2DEpipolarSampleCamDepth.Load(uint3(uiSampleNum, uiSliceNum, 0));
    float fDepth = camera.Proj[2][2] + camera.Proj[3][2] / fCamDepth;
    return ComputeUnshadowedSampleScatter(f2ScreenXY, fCamDepth, fDepth);
}

technique11 ComputeUnshadowedSampleScatterTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, DoUnshadowedRayMarch()));
    }
}

RWTexture2D<uint2> g_rwtex2DInterpolationSource : register(u0);
static const uint SAMPLE_STEP = 16;
static const uint THREAD_GROUP_SIZE = 128;
static const uint g_uiPackNum = THREAD_GROUP_SIZE / 32;
groupshared uint g_uiCamDepthDiffPackFlags[g_uiPackNum];

static const float m_fRefinementThreshold = 0.03;
static const float m_fSampleDense = 2;

[numthreads(THREAD_GROUP_SIZE, 1, 1)]
void RefineSampleCS(uint3 Gid : SV_GroupID,
                    uint3 GTid : SV_GroupThreadID)
{
    uint uiSliceNum = Gid.y;
    uint uiSampleGroupStart = Gid.x * THREAD_GROUP_SIZE;
    uint uiSampleGlobalNum = uiSampleGroupStart + GTid.x;
    
    float2 f2SampleScreenXY = g_tex2DEpipolarSample.Load(uint3(uiSampleGlobalNum, uiSliceNum, 0));

    bool IsValidThread = all(abs(f2SampleScreenXY) < 1 + 1e-4);

    if (GTid.x < g_uiPackNum)
        g_uiCamDepthDiffPackFlags[GTid.x] = 0;
    GroupMemoryBarrierWithGroupSync();

    [branch]
    if (IsValidThread)
    {
#if USE_DEPTH_WEIGHT
        float fSampleCamDepth = g_tex2DEpipolarSampleCamDepth.Load(uint3(uiSampleGlobalNum, uiSliceNum, 0));
        float fSampleCamDepthRight = g_tex2DEpipolarSampleCamDepth.Load(uint3(min(uiSampleGlobalNum + 1, EPIPOLAR_SAMPLE_NUM - 1), uiSliceNum, 0));
        float fMax = max(fSampleCamDepth, fSampleCamDepthRight);
        fMax = max(fMax, 1);
        bool bFlag = abs(fSampleCamDepth - fSampleCamDepthRight) / fMax < 0.2 * m_fRefinementThreshold; // 1 No break
#else
        float3 f3Insctr0 = g_tex2DUnshadowedSampleScatter.Load(uint3(uiSampleGlobalNum, uiSliceNum, 0));
        float3 f3Insctr1 = g_tex2DUnshadowedSampleScatter.Load(uint3(uiSampleGlobalNum + 1, uiSliceNum, 0));
        float3 f3MaxInsctr = max(f3Insctr0, f3Insctr1);

        float fAveLogLum = 0.1;
        const float middleGray = 1.03 - 2 / (2 + log10(fAveLogLum + 1));
        float3 f3MinInsctrThreshold = (0.02f * fAveLogLum.xxx / RGB_TO_LUMINANCE.xyz) / middleGray;

        f3MaxInsctr = max(f3MaxInsctr, f3MinInsctrThreshold);
        bool bFlag = all((abs(f3Insctr0 - f3Insctr1) / f3MaxInsctr) < m_fRefinementThreshold);
#endif
        InterlockedOr(g_uiCamDepthDiffPackFlags[GTid.x / 32], bFlag << (GTid.x % 32));
    }
    GroupMemoryBarrierWithGroupSync();

    if (!IsValidThread)
        return;
    
    uint uiSampleStep = SAMPLE_STEP;
    uint uiSampleLocalNum0 = (GTid.x / uiSampleStep) * uiSampleStep;
    uint uiSampleGlobalNum0 = uiSampleLocalNum0 + uiSampleGroupStart;
    float2 f2SampleScreenXY0 = g_tex2DEpipolarSample.Load(uint3(uiSampleGlobalNum0, uiSliceNum, 0));
    if (length(f2SampleScreenXY0 - light.f4LightScreenPos.xy) < 0.1 &&
        (float) uiSampleGlobalNum0 / (float) EPIPOLAR_SAMPLE_NUM < 0.05)
    {
        uiSampleStep = max(uiSampleStep / m_fSampleDense, 1);
        uiSampleLocalNum0 = (GTid.x / uiSampleStep) * uiSampleStep;
    }
    uint uiSampleLocalNum1 = uiSampleLocalNum0 + uiSampleStep;
    if (Gid.x == EPIPOLAR_SAMPLE_NUM / THREAD_GROUP_SIZE - 1)
        uiSampleLocalNum1 = min(uiSampleLocalNum1, THREAD_GROUP_SIZE - 1);
    uiSampleStep = uiSampleLocalNum1 - uiSampleLocalNum0;

    uint uiSampleLeftBoundary, uiSampleRightBoundary;
    if (GTid.x > uiSampleLocalNum0 && GTid.x < uiSampleLocalNum1)
    {
        uint uiCamDepthDiffPackFlags[g_uiPackNum];
        for (uint i = 0; i < g_uiPackNum; i++)
        {
            uiCamDepthDiffPackFlags[i] = g_uiCamDepthDiffPackFlags[i];
        }
    
        bool bNoBreak = true;

        int iPackNum0 = uiSampleLocalNum0 / 32;
        int iNumInPack0 = uiSampleLocalNum0 % 32;
    
        int iPackNum1 = uiSampleLocalNum1 / 32;
        int iNumInPack1 = uiSampleLocalNum1 % 32;

        for (int j = iPackNum0; j <= iPackNum1; j++)
        {
            if (uiCamDepthDiffPackFlags[j] != 0xFFFFFFFFU)
            {
                bNoBreak = false;
                break;
            }
        }

        if (bNoBreak)
        {
            uiSampleLeftBoundary = uiSampleLocalNum0;
            uiSampleRightBoundary = uiSampleLocalNum1;
        }
        else
        {
            //uint uiSampleLeft = GTid.x;
            int iSampleLeft = GTid.x - 1;
            int iSampleLeftPackNum = uint(iSampleLeft) / 32;

            int iFirstBreakFlag = -1;
            while (iFirstBreakFlag == -1 && iSampleLeftPackNum >= iPackNum0)
            {
                uint uiFlag = uiCamDepthDiffPackFlags[iSampleLeftPackNum];
                int iNumInPackSampleLeft = uint(iSampleLeft) % 32;

                if (iNumInPackSampleLeft < 31)
                {
                    uiFlag |= (uint(0x0FFFFFFFFU) << uint(iNumInPackSampleLeft + 1));
                }
                iFirstBreakFlag = firstbithigh(uint(~uiFlag));
                if (!(iFirstBreakFlag >= 0 && iFirstBreakFlag <= 31))
                    iFirstBreakFlag = -1;
                iSampleLeft -= iNumInPackSampleLeft - iFirstBreakFlag;
                iSampleLeftPackNum--;
            }
            uiSampleLeftBoundary = max(uiSampleLocalNum0, uint(iSampleLeft + 1));

            //uint uiSampleRight = GTid.x - 1; // if last break this right should be this
            int iSampleRight = GTid.x;
            int iSampleRightPackNum = uint(iSampleRight) / 32;
            iFirstBreakFlag = 32;
            while (iFirstBreakFlag == 32 && iSampleRightPackNum <= iPackNum1)
            {
                uint uiFlag = uiCamDepthDiffPackFlags[iSampleRightPackNum];
                int iNumInPackSampleRight = uint(iSampleRight) % 32;

                if (iNumInPackSampleRight > 0)
                {
                    uiFlag |= ((1 << uint(iNumInPackSampleRight)) - 1);
                }
                iFirstBreakFlag = firstbitlow(uint(~uiFlag));
                if (!(iFirstBreakFlag >= 0 && iFirstBreakFlag <= 31))
                    iFirstBreakFlag = 32;
                iSampleRight += iFirstBreakFlag - iNumInPackSampleRight;
                iSampleRightPackNum++;
            
            }
            uiSampleRightBoundary = min(uiSampleLocalNum1, uint(iSampleRight));

            if (uiSampleLeftBoundary == GTid.x || uiSampleRightBoundary == GTid.x)
                uiSampleLeftBoundary = uiSampleRightBoundary = GTid.x;
        }
        //else
            //g_rwtex2DInterpolationSource[uint2(uiSampleGlobalNum, uiSliceNum)] = uint4(uiSampleGroupStart + iSampleLeftBoundary, uiSampleGroupStart + iSampleRightBoundary, 0, 1);
    }
    else
    {
        uiSampleLeftBoundary = uiSampleRightBoundary = GTid.x;
        //g_rwtex2DInterpolationSource[uint2(uiSampleGlobalNum, uiSliceNum)] = uint4(0,0,0,1);
    }
    g_rwtex2DInterpolationSource[uint2(uiSampleGlobalNum, uiSliceNum)] = uint2(uiSampleGroupStart + uiSampleLeftBoundary, uiSampleGroupStart + uiSampleRightBoundary);
}

technique11 RefineSampleTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(NULL);
        SetGeometryShader(NULL);
        SetPixelShader(NULL);
        SetComputeShader(CompileShader(cs_5_0, RefineSampleCS()));
    }
}

float4 ComputeSliceUVOrigDirTex2D(QuadVertexOut In) : SV_Target
{
    float4 f4SliceEnd = g_tex2DSliceEnd.Load(uint3(In.m_f4Pos.x,0,0));
    if (!IsValidScreenLocation(f4SliceEnd.xy))
        return INVALID_EPIPOLAR_LINE;
    float4 f4SliceEndInWorldSpace = mul(float4(f4SliceEnd.zw, 1.0,1.0), camera.InvViewProj);
    f4SliceEndInWorldSpace /= f4SliceEndInWorldSpace.w;
    float4 f4SliceEndInLightSpace = mul(f4SliceEndInWorldSpace, light.ViewProj);
    float2 f2SliceUVEnd = ProjToUV(f4SliceEndInLightSpace.xy / f4SliceEndInLightSpace.w);

    float4 f4SliceStartInLightSpace = mul(float4(camera.f3CameraPos, 1.0), light.ViewProj);
    float2 f2SliceUVOrig = ProjToUV(f4SliceStartInLightSpace.xy / f4SliceStartInLightSpace.w);
    float2 f2SliceDir = f2SliceUVEnd - f2SliceUVOrig;
    float fSliceDirLength = length(f2SliceDir);
    f2SliceDir /= fSliceDirLength;
    
    float2 f2ShadowMapDim = float2(SHADOWMAP_TEXTURE_DIM, SHADOWMAP_TEXTURE_DIM);
    float4 f4Boundary = float4(0, 0, 1, 1) + float4(0.5, 0.5, -0.5, -0.5) / f2ShadowMapDim.xyxy;

    float4 f4DistToBoundary = (f2SliceUVOrig.xyxy - f4Boundary) * float4(1, 1, -1, -1);
    if (any(f4DistToBoundary<0))
    {
        bool4 b4IsCorrectIntersectionFlag = abs(f2SliceDir.xyxy) > 1e-5;
        float4 f4DirDistToBoundary = (f4Boundary - f2SliceUVOrig.xyxy) / (f2SliceDir.xyxy + !b4IsCorrectIntersectionFlag);

        float4 f4InsecBoundary = f2SliceUVOrig.yxyx + f2SliceDir.yxyx * f4DirDistToBoundary;
        b4IsCorrectIntersectionFlag = b4IsCorrectIntersectionFlag && 
                                        (f4InsecBoundary >= f4Boundary.yxyx) && 
                                        (f4InsecBoundary <= f4Boundary.wzwz) &&
                                        (f4DirDistToBoundary >= 0);

        f4DirDistToBoundary = f4DirDistToBoundary * b4IsCorrectIntersectionFlag + 
                                float4(+FLT_MAX, +FLT_MAX, +FLT_MAX, +FLT_MAX) * !b4IsCorrectIntersectionFlag;
        float2 f2FirstDistToBoundary = min(f4DirDistToBoundary.xy, f4DirDistToBoundary.zw);
        float fFirstDistToBoundary = min(f2FirstDistToBoundary.x, f2FirstDistToBoundary.y);
        f2SliceUVOrig += f2SliceDir * fFirstDistToBoundary;
    }
    f2SliceDir /= f2ShadowMapDim.xy;
    return float4(f2SliceUVOrig, f2SliceDir);
}

technique11 ComputeSliceUVOrigDirTex2DTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeSliceUVOrigDirTex2D()));
    }
}

float2 Initial1DMinMaxMipMap(QuadVertexOut In) : SV_Target
{
    float4 f4SliceUVOrigDir = g_tex2DSliceUVOrigDir.Load(uint3(In.m_f4Pos.y, 0, 0));
    float2 f2CurrUV = f4SliceUVOrigDir.xy + f4SliceUVOrigDir.zw * floor(In.m_f4Pos.x) * 2.f;

    float4 f4MinDepth = 1;
    float4 f4MaxDepth = 0;
    for (int i = 0; i <= 1; i++)
    {
        float4 f4Depths = g_tex2DShadowMap.Gather(samLinearBorder0, f2CurrUV + i * f4SliceUVOrigDir.zw);
        if (any(f4Depths == 0))
        {
            uint4 ui4IsCorrect = f4Depths == float4(0, 0, 0, 0);
            float4 f4DepthsNew = f4Depths * !ui4IsCorrect + ui4IsCorrect;
            f4MinDepth = min(f4MinDepth, f4DepthsNew);
        }
        else
            f4MinDepth = min(f4MinDepth, f4Depths);
        f4MaxDepth = max(f4MaxDepth, f4Depths);
    }
    float2 f2MinDepth = min(f4MinDepth.xy, f4MinDepth.zw);
    float2 f2MaxDepth = max(f4MaxDepth.xy, f4MaxDepth.zw);
    float fMinDepth = min(f2MinDepth.x, f2MinDepth.y);
    float fMaxDepth = max(f2MaxDepth.x, f2MaxDepth.y);

    return float2(fMinDepth,fMaxDepth);
}

technique11 Initial1DMinMaxMipMapTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, Initial1DMinMaxMipMap()));
    }
}

float2 Compute1DMinMaxMipMapLevel(QuadVertexOut In) : SV_Target
{
    uint2 ui2SrcMinMaxUV0 = uint2(misc.ui4SrcDstMinMaxOffset.x + (uint(In.m_f4Pos.x) - misc.ui4SrcDstMinMaxOffset.z) * 2, 
                                    In.m_f4Pos.y);
    uint2 ui2SrcMinMaxUV1 = ui2SrcMinMaxUV0 + uint2(1, 0);

    float2 f2MinMaxDepth0 = g_tex2DMinMaxMipMap.Load(uint3(ui2SrcMinMaxUV0, 0));
    float2 f2MinMaxDepth1 = g_tex2DMinMaxMipMap.Load(uint3(ui2SrcMinMaxUV1, 0));

    float fMinDepth = min(f2MinMaxDepth0.x, f2MinMaxDepth1.x);
    float fMaxDepth = max(f2MinMaxDepth0.y, f2MinMaxDepth1.y);

    return float2(fMinDepth, fMaxDepth);
}

technique11 Compute1DMinMaxMipMapLevelTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, Compute1DMinMaxMipMapLevel()));
    }
}

void MarkRayMarchSample(QuadVertexOut In) 
{
    uint2 uiInterpolationSample = g_tex2DInterpolationSample.Load(uint3(In.m_f4Pos.xy,0));
    if(uiInterpolationSample.x != uiInterpolationSample.y)
        discard;
}

technique11 MarkRayMarchSampleTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest_StEqual_IncrStencil, 1);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, MarkRayMarchSample()));
    }
}

float4 ComputeShadowInscatter(float2 f2ScreenXY,float fRayEndCamDepth,float fDepth,bool bIsUseMinMaxMap,uint uiSliceNum)
{
    float4 f4ViewRayEndPos = mul(float4(f2ScreenXY, fDepth, 1.f), camera.InvViewProj);
    float3 f3ViewRay = f4ViewRayEndPos.xyz / f4ViewRayEndPos.w - camera.f3CameraPos;
    float fViewRayLenInWorldSpace = length(f3ViewRay);
    f3ViewRay /= fViewRayLenInWorldSpace;

    float3 f3EarthCenter = float3(0, -atmosphere.bottom_radius, 0);
    float3 f3CameraPos = camera.f3CameraPos - f3EarthCenter;
    
    float4 f4RMuMuSNu = 0;
    f4RMuMuSNu.xyz = GetRMuMuS(f3CameraPos, f3ViewRay);
    f4RMuMuSNu.w = dot(f3ViewRay, light.f3LightDir);

    bool bIsNoScatter, bIsMarchToAtmosphere, bIsMarchToEarth, bIsIntersectEarth;
    float fDistToAtmosphereNear, fDistToAtmosphereFar, fDistToEarthNear;
    fViewRayLenInWorldSpace = GetRayMarchLen(f4RMuMuSNu, fRayEndCamDepth, fViewRayLenInWorldSpace,
                                            bIsNoScatter, bIsMarchToAtmosphere, bIsMarchToEarth, bIsIntersectEarth,
                                            fDistToAtmosphereNear, fDistToAtmosphereFar, fDistToEarthNear);
    
#if USE_SHADOW_OBJECT_TO_EARTH
    if (bIsMarchToAtmosphere && f4RMuMuSNu.w <= 0)
#else 
    if ((bIsMarchToAtmosphere && f4RMuMuSNu.w <= 0) || bIsMarchToEarth)
#endif
    {
        return float4(m_f4ExtraLight * GetSkyMultiScatter(f4RMuMuSNu.x, f4RMuMuSNu.y, f4RMuMuSNu.z, f4RMuMuSNu.w), 1);
    }

    float3 f3StartPos = camera.f3CameraPos;
    float3 f3EndPos = camera.f3CameraPos + fViewRayLenInWorldSpace * f3ViewRay;
    if (fDistToAtmosphereNear > 0.f)
    {
        float3 f3DistToAtmosphereNear = fDistToAtmosphereNear * f3ViewRay;
        f3CameraPos += f3DistToAtmosphereNear;
        f3StartPos += f3DistToAtmosphereNear;
        fViewRayLenInWorldSpace -= fDistToAtmosphereNear;
        fDistToAtmosphereFar -= fDistToAtmosphereNear;
        if (bIsMarchToEarth)
            fDistToEarthNear -= fDistToAtmosphereNear;
        f4RMuMuSNu.xyz = GetRMuMuS(f3CameraPos, f3ViewRay);
    }

    float4 f4StartPosInShadowMap = mul(float4(f3StartPos, 1), light.ViewProj);
    float3 f3StartUVDepthInShadowMap = f4StartPosInShadowMap.w == 0 ? f4StartPosInShadowMap.xyz : 
                                                                        f4StartPosInShadowMap.xyz / f4StartPosInShadowMap.w;
    f3StartUVDepthInShadowMap.xy = ProjToUV(f3StartUVDepthInShadowMap.xy);

    float4 f4EndPosInShadowMap = mul(float4(f3EndPos, 1), light.ViewProj);
    float3 f3EndUVDepthInShadowMap = f4EndPosInShadowMap.w == 0 ? f4EndPosInShadowMap.xyz : 
                                                                    f4EndPosInShadowMap.xyz / f4EndPosInShadowMap.w;
    f3EndUVDepthInShadowMap.xy = ProjToUV(f3EndUVDepthInShadowMap.xy);
    
    float3 f3ViewRayInShadowMap = f3EndUVDepthInShadowMap - f3StartUVDepthInShadowMap;
    float fViewRayLenInShadowMap = length(f3ViewRayInShadowMap.xy);
    fViewRayLenInShadowMap = max(fViewRayLenInShadowMap, 1e-7);
    f3ViewRayInShadowMap /= fViewRayLenInShadowMap;

    float4 f4SliceUVOrigDir = 0;
    float fStepLenInShadowMap = 0;
    float fMaxTraceDirDim = max(abs(f3ViewRayInShadowMap.x), abs(f3ViewRayInShadowMap.y));
    if (bIsUseMinMaxMap)
    {
        f4SliceUVOrigDir = g_tex2DSliceUVOrigDir.Load(uint3(uiSliceNum, 0, 0));
        fStepLenInShadowMap = length(f4SliceUVOrigDir.zw);
    }
    else
    {
        fStepLenInShadowMap = (fMaxTraceDirDim > 0) ? (1 / (float)SHADOWMAP_TEXTURE_DIM / fMaxTraceDirDim):0;
            // Take into account maximum number of steps specified by the g_MiscParams.fMaxStepsAlongRay
        fStepLenInShadowMap = max(fViewRayLenInShadowMap / 256, fStepLenInShadowMap);
    }
    
    float fStepLenInWorldSpace = fViewRayLenInWorldSpace * (fStepLenInShadowMap / fViewRayLenInShadowMap);
    float3 f3StepInShadowMap = f3ViewRayInShadowMap * fStepLenInShadowMap;
    float fMaxStepScale = 1 << misc.uiMinMaxLevelMax;

    float fMaxAllowedWorldStepLen = fViewRayLenInWorldSpace / 10;
    fMaxStepScale = min(fMaxStepScale, fMaxAllowedWorldStepLen / fStepLenInWorldSpace);
    
    if (fStepLenInWorldSpace > fMaxAllowedWorldStepLen)
    {
        fStepLenInWorldSpace = fMaxAllowedWorldStepLen;
        fStepLenInShadowMap = fViewRayLenInShadowMap * fStepLenInWorldSpace / fViewRayLenInWorldSpace;
        f3StepInShadowMap = f3ViewRayInShadowMap * fStepLenInShadowMap;
        fMaxStepScale = -1;
    }

    float fTotalLightLen = 0;
    float fTotalMarchLen = 0;
    float fDistToFirstLight = -1;
    uint uiCurrSamplePos = 0;
       
    uiCurrSamplePos = length(f3StartUVDepthInShadowMap.xy - f4SliceUVOrigDir.xy) / fStepLenInShadowMap + 0.5;
    
    if (any(f3StartUVDepthInShadowMap.xy < 0) || any(f3StartUVDepthInShadowMap.xy > 1))
    {
        fDistToFirstLight = 0;
        fTotalMarchLen = fTotalLightLen = min(uiCurrSamplePos * fStepLenInWorldSpace, fViewRayLenInWorldSpace);
        f3StartUVDepthInShadowMap += uiCurrSamplePos * f3StepInShadowMap;
        uiCurrSamplePos = length(f3StartUVDepthInShadowMap.xy - f4SliceUVOrigDir.xy) / fStepLenInShadowMap + 0.5;
        //return float4(1, 0, 0, 1);
    }

    float3 f3CurrSampleUVDepthInShadowMap = f3StartUVDepthInShadowMap;
    uint uiCurrLevel = 0;
    int iCurrLevelOffsetX = -(int)(MIN_MAX_TEXTURE_DIM);
    float fScale = 1.f;
    
    [loop]
    while (fTotalMarchLen < fViewRayLenInWorldSpace)
    {
        float bIsLight = 0;
        if(bIsUseMinMaxMap)
        {      
            while (uiCurrLevel < misc.uiMinMaxLevelMax && ((uiCurrSamplePos & ((2 << uiCurrLevel) - 1)) == 0)
                    && fScale * 2 < fMaxStepScale)
            {
                iCurrLevelOffsetX += (int) (MIN_MAX_TEXTURE_DIM >> uiCurrLevel);
                fScale *= 2.f;
                uiCurrLevel++;
            }
          
            while(uiCurrLevel > 0)
            {
                float fEndDepthInLevel = f3CurrSampleUVDepthInShadowMap.z + f3StepInShadowMap.z * (fScale - 1);
                float2 f2StartEndDepthInLevel = max(float2(f3CurrSampleUVDepthInShadowMap.z, fEndDepthInLevel),1e-7);
                
                float2 f2MinMaxDepth = g_tex2DMinMaxMipMap.Load(uint3(iCurrLevelOffsetX + (uiCurrSamplePos >> uiCurrLevel), uiSliceNum, 0));
                
                bIsLight = all(f2StartEndDepthInLevel <= f2MinMaxDepth.xx);
                bool bIsShadow = all(f2StartEndDepthInLevel > f2MinMaxDepth.yy);
                if(bIsLight || bIsShadow)
                    break;
                uiCurrLevel--;
                fScale /= 2.f;
                iCurrLevelOffsetX -= (int) (MIN_MAX_TEXTURE_DIM >> uiCurrLevel);
            }

            if (uiCurrLevel == 0)
            {
                bIsLight = g_tex2DShadowMap.SampleLevel(samLinearBorder1, f3CurrSampleUVDepthInShadowMap.xy, 0) >= max(f3CurrSampleUVDepthInShadowMap.z, 1e-7);

            }
        }
        else
        {
            bIsLight = g_tex2DShadowMap.SampleLevel(samLinearBorder1, f3CurrSampleUVDepthInShadowMap.xy, 0) > max(f3CurrSampleUVDepthInShadowMap.z, 1e-7);
        }
        
        if (fDistToFirstLight < 0 && bIsLight > 0)
        {
            fDistToFirstLight = fTotalMarchLen;
        }
        uiCurrSamplePos += (1 << uiCurrLevel);
        f3CurrSampleUVDepthInShadowMap += (fScale * f3StepInShadowMap);

        float fRemainDist = max(fViewRayLenInWorldSpace - fTotalMarchLen, 0);
        float fMarchLenInWorldSpace = min(fStepLenInWorldSpace * fScale, fRemainDist);

        fTotalMarchLen += fMarchLenInWorldSpace;
        fTotalLightLen += (fMarchLenInWorldSpace * bIsLight);
        if (any(f3CurrSampleUVDepthInShadowMap.xy) < 0 ||
            any(f3CurrSampleUVDepthInShadowMap.xy) > 1 ||
            uiCurrSamplePos >= SHADOWMAP_TEXTURE_DIM)
        {
            float fRemainDist = max(fViewRayLenInWorldSpace - fTotalMarchLen, 0);
            if (fDistToFirstLight < 0)
                fDistToFirstLight = fTotalMarchLen;
            fTotalMarchLen += fRemainDist;
            fTotalLightLen += fRemainDist;
            //return float4(1, 0, 0, 1);
        }
    }

    if (fDistToFirstLight < 0)
    {
        fDistToFirstLight = fTotalMarchLen;
    }

    float3 f3MultiScatter = 0;
    float3 f3GroundRadiance = 0;

    //return float4(fDistToFirstLight, fTotalLightLen, fTotalMarchLen, 1);
    [branch]
    if(fTotalLightLen>0)
    {
        //if (bIsMarchToEarth)
        ////if (fViewRayLenInWorldSpace == fDistToEarthNear)
        //{
        //    float3 f3EarthNormal = f3CameraPos + f3ViewRay * fDistToEarthNear;
        //    float3 f3EarthRMuMuS = GetRMuMuS(f3EarthNormal, f3ViewRay);
        //    f3EarthNormal /= f3EarthRMuMuS.x;

        //    float2 f2EarthCoords = float2(f3EarthNormal.x == 0 ? 0 : atan2(f3EarthNormal.z, f3EarthNormal.x), acos(f3EarthNormal.y)) * float2(0.5, 1.0) / PI + float2(0.5, 0.0);
        //    float4 f4Reflentance = g_tex2DEarthGround.SampleLevel(samLinearClamp, f2EarthCoords, 0); //* float4(0.2, 0.2, 0.2, 1.0);
        
        //    float2 f2IrradianceUV = GetIrradianceUVFromRMuS(f3EarthRMuMuS.x, f3EarthRMuMuS.z);

        //     // R[L0] + R[L*]
        //    float3 f3SunColor = g_tex2DDirectIrradianceLUT.SampleLevel(samLinearClamp, f2IrradianceUV,0) * max(f3EarthRMuMuS.z, 0);
        //    float3 f3SkyColor = g_tex2DIndirectIrradianceLUT.SampleLevel(samLinearClamp, f2IrradianceUV,0) * (1.0 + dot(f3EarthNormal, f3CameraPos) / f4RMuMuSNu.x) * 0.5;
          
        //    // S[L]x - T(x,xs)S[L]xs=x0-lv
        //    //float fShadowLength = fDistToEarthNear - fTotalLightLen;
        //    float fShadowLength = 0;
        //    float3 f3TransmittanceToEarth = 1;
        //    f3MultiScatter = GetSkyMultiScatterToGround(f4RMuMuSNu.x, f4RMuMuSNu.y, f4RMuMuSNu.z,
        //                                                f3EarthRMuMuS.x, f3EarthRMuMuS.y, f3EarthRMuMuS.z,
        //                                                f4RMuMuSNu.w, fShadowLength, fDistToEarthNear,f3TransmittanceToEarth);
        //    //f3MultiScatter = GetSkyMultiScatter(f4RMuMuSNu.x, f4RMuMuSNu.y, f4RMuMuSNu.z, f4RMuMuSNu.w);
        //    f3MultiScatter *= m_f4ExtraLight;
        //    //f3GroundRadiance = f4Reflentance.rgb * atmosphere.ground_albedo * (1.0 / PI) * (f3SunColor + f3SkyColor) * f3TransmittanceToEarth;        
        //}
        //else if (bIsMarchToAtmosphere)
        ////else if(fViewRayLenInWorldSpace == fDistToAtmosphereFar)
        //{
        //    // T(x,xs)S[L]xs=x+lv
        //    float fShadowLength = fDistToAtmosphereFar - fTotalLightLen;
        //    f3MultiScatter = GetSkyMultiScatterToAtmosphere(f4RMuMuSNu.x, f4RMuMuSNu.y, f4RMuMuSNu.z, f4RMuMuSNu.w,
        //                                                    fShadowLength);
        //    f3MultiScatter *= m_f4ExtraLight;
        //}
        //else
        //{
            f3StartPos = f3CameraPos + fDistToFirstLight * f3ViewRay;
            f3EndPos = f3StartPos + fTotalLightLen * f3ViewRay;
            
            float3 f3StartRMuMuS = GetRMuMuS(f3StartPos, f3ViewRay);
            float3 f3EndRMuMuS = GetRMuMuS(f3EndPos, f3ViewRay);
            
            //float3 f3TransmittanceToStart = (fDistToFirstLight == 0) ? 1 :
            //                                GetTransmittance(f4RMuMuSNu.x, f4RMuMuSNu.y, f3StartRMuMuS.x, f3StartRMuMuS.y, bIsIntersectEarth);
            float3 f3TransmittanceToStart = (fDistToFirstLight == 0) ? 1 :
                                         GetTransmittanceIntegralAnalytic(f4RMuMuSNu.x, f4RMuMuSNu.y, fDistToFirstLight);
            float fDistToEnd = fDistToFirstLight + fTotalLightLen;
            float3 f3TransmittanceToEnd = (fDistToEnd == 0) ? 1 :
                                        GetTransmittanceIntegralAnalytic(f4RMuMuSNu.x, f4RMuMuSNu.y, fDistToEnd);
            //float3 f3TransmittanceToEnd = (fDistToEnd == 0) ? 1 :
            //                            GetTransmittance(f4RMuMuSNu.x, f4RMuMuSNu.y, f3EndRMuMuS.x, f3EndRMuMuS.y, bIsIntersectEarth);
            
            f3MultiScatter = f3TransmittanceToStart * GetSkyMultiScatter(f3StartRMuMuS.x, f3StartRMuMuS.y, f3StartRMuMuS.z, f4RMuMuSNu.w);

            f3MultiScatter -= f3TransmittanceToEnd * GetSkyMultiScatter(f3EndRMuMuS.x, f3EndRMuMuS.y, f3EndRMuMuS.z, f4RMuMuSNu.w);
            //if (misc.fIsLightInSpaceCorrect)
                f3MultiScatter *= m_f4ExtraLight;
        //}
    }
    return float4(f3MultiScatter, 1);
}

float4 DoRayMarch(QuadVertexOut In) : SV_Target
{
    uint2 ui2XY = uint2(In.m_f4Pos.xy);
    float2 f2XY = g_tex2DEpipolarSample.Load(uint3(ui2XY, 0));
    float fCamDepth = g_tex2DEpipolarSampleCamDepth.Load(uint3(ui2XY, 0));
    float fDepth = camera.Proj[2][2] + camera.Proj[3][2] / fCamDepth;

    return ComputeShadowInscatter(f2XY, fCamDepth, fDepth,true, ui2XY.y);
}

technique11 DoRayMarchTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest_StEqual_KeepStencil, 2);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, DoRayMarch()));
    }
}

float4 InterpolateScatter(QuadVertexOut In) : SV_Target
{
    uint uiSampleNum = In.m_f4Pos.x;
    uint uiSliceNum = In.m_f4Pos.y;
    
    uint2 ui2InterpolationSample = g_tex2DInterpolationSample.Load(uint3(uiSampleNum, uiSliceNum, 0.0));
    float fInterpolated = float(uiSampleNum - ui2InterpolationSample.x) / max(float(ui2InterpolationSample.y - ui2InterpolationSample.x), 1.f);
    
    float3 f3InterpolatedScatter0, f3InterpolatedScatter1;
    if(misc.fEnableLightShaft)
    {
        f3InterpolatedScatter0 = g_tex2DSampleScatter.Load(uint3(ui2InterpolationSample.x, uiSliceNum, 0.0));
        f3InterpolatedScatter1 = g_tex2DSampleScatter.Load(uint3(ui2InterpolationSample.y, uiSliceNum, 0.0));
    }
    else
    {
        f3InterpolatedScatter0 = g_tex2DUnshadowedSampleScatter.Load(uint3(ui2InterpolationSample.x, uiSliceNum, 0.0));
        f3InterpolatedScatter1 = g_tex2DUnshadowedSampleScatter.Load(uint3(ui2InterpolationSample.y, uiSliceNum, 0.0));
    }
    float3 f3InterpolatedScatter = lerp(f3InterpolatedScatter0, f3InterpolatedScatter1, fInterpolated);

    return float4(f3InterpolatedScatter, 1);
}

technique11 InterpolateScatterTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest_StEqual_KeepStencil, 1);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, InterpolateScatter()));
    }
}

bool ComputeIsIntersectEarth(float2 f2XY)
{
    bool bIsIntersectEarth = false;
    float fCamDepth = g_tex2DSpaceLinearDepth.SampleLevel(samLinearClamp, ProjToUV(f2XY), 0);
    float fDepth = camera.Proj[2][2] + camera.Proj[3][2] / fCamDepth;
    float4 f4RayEndInWorldSpace = mul(float4(f2XY, fDepth, 1.f), camera.InvViewProj);
    float3 f3ViewRay = f4RayEndInWorldSpace.xyz / f4RayEndInWorldSpace.w - camera.f3CameraPos;
    float3 f3EarthCenter = float3(0, -atmosphere.bottom_radius, 0);
    float3 f3CameraPos = camera.f3CameraPos - f3EarthCenter;
    float fHeight = length(f3CameraPos);
    float fCosZenithAngle = dot(f3CameraPos, f3ViewRay) / fHeight;
    float fRMu = fHeight * fCosZenithAngle;
    float fIntersectEarthDiscriminant = fRMu * fRMu - fHeight * fHeight + atmosphere.bottom_radius * atmosphere.bottom_radius;
    if (fIntersectEarthDiscriminant >= 0 && fCamDepth > camera.fFarZ && fCosZenithAngle < 0)
        bIsIntersectEarth = true;
    return bIsIntersectEarth;
}

float4 ApplyInterpolateScatter(QuadVertexOut In) : SV_Target
{
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float fCamDepth = g_tex2DSpaceLinearDepth.SampleLevel(samLinearClamp, f2UV, 0);
    float fDepth = camera.Proj[2][2] + camera.Proj[3][2] / fCamDepth;

    bool bIsIntersectEarth = false;
    float4 f4RayEndInWorldSpace = mul(float4(In.m_f2PosPS, fDepth, 1.f), camera.InvViewProj);
    float3 f3ViewRay = f4RayEndInWorldSpace.xyz / f4RayEndInWorldSpace.w - camera.f3CameraPos;
    float fViewRayLen = length(f3ViewRay);
    f3ViewRay /= fViewRayLen;
    float3 f3EarthCenter = float3(0, -atmosphere.bottom_radius, 0);
    float3 f3CameraPos = camera.f3CameraPos - f3EarthCenter;
    float fHeight = length(f3CameraPos);
    float fCosZenithAngle = dot(f3CameraPos, f3ViewRay) / fHeight;
    float fRMu = fHeight * fCosZenithAngle;
    float fIntersectEarthDiscriminant = fRMu * fRMu - fHeight * fHeight + atmosphere.bottom_radius * atmosphere.bottom_radius;
    if (fIntersectEarthDiscriminant >= 0 && fCamDepth > camera.fFarZ && fCosZenithAngle < 0)
        bIsIntersectEarth = true;

    float2 f2RayDir = normalize(In.m_f2PosPS - light.f4LightScreenPos.xy);

    float2 f2ScreenDim = float2(SCREEN_WIDTH, SCREEN_HEIGHT);
    float4 f4Boundary = float4(-1, -1, 1, 1) + float4(1, 1, -1, -1) / f2ScreenDim.xyxy;
    
    bool4 b4IsCorrectIntersectionFlag = abs(f2RayDir.xyxy) > 1e-5;
    float4 f4DistToBoundary = (f4Boundary - In.m_f2PosPS.xyxy) / (f2RayDir.xyxy + !b4IsCorrectIntersectionFlag);
    float4 f4InsecBoundary = In.m_f2PosPS.yxyx + f4DistToBoundary * f2RayDir.yxyx;
    b4IsCorrectIntersectionFlag = b4IsCorrectIntersectionFlag && 
                                (f4DistToBoundary > 0) &&
                                (f4InsecBoundary >= f4Boundary.yxyx) &&
                                (f4InsecBoundary <= f4Boundary.wzwz);
    

    float4 f4EpipolarSlice = float4(0, 0.25, 0.5, 0.75) + (f4InsecBoundary - f4Boundary.wxyz) * float4(-1, 1, 1, -1) / (f4Boundary.wzwz - f4Boundary.yxyx) / 4.f;
    float fEpipolarSlice = frac(dot(f4EpipolarSlice, b4IsCorrectIntersectionFlag));

    //float4 f4HalfSpaceEquationTerms = (In.m_f2PosPS.xxyy - f4Boundary.xzyw/*float4(-1,1,-1,1)*/) * f2RayDir.yyxx;
    //bool4 b4HalfSpaceFlags = f4HalfSpaceEquationTerms.xyyx < f4HalfSpaceEquationTerms.zzww;
    //
    //bool4 b4SectorFlags = b4HalfSpaceFlags.wxyz && !b4HalfSpaceFlags.xyzw;
    //
    //float4 f4DistToBoundaries = (f4Boundary - light.f4LightScreenPos.xyxy) / (f2RayDir.xyxy + float4(abs(f2RayDir.xyxy) < 1e-6));
    //// Select distance to the exit boundary:
    //float fDistToExitBoundary = dot(b4SectorFlags, f4DistToBoundaries);
    //// Compute exit point on the boundary:
    //float2 f2ExitPoint = light.f4LightScreenPos.xy + f2RayDir * fDistToExitBoundary;
    
    //float4 f4EpipolarSlice = float4(0, 0.25, 0.5, 0.75) +
    //    saturate((f2ExitPoint.yxyx - f4Boundary.wxyz) * float4(-1, +1, +1, -1) / (f4Boundary.wzwz - f4Boundary.yxyx)) / 4.0;
    //// Select the right value:
    //float fEpipolarSlice = frac(dot(b4SectorFlags, f4EpipolarSlice));

    float fSliceNum = fEpipolarSlice * ((float) EPIPOLAR_SLICE_NUM - 1.f);
    float fSliceNum0 = floor(fSliceNum);
    
    float fSliceNumInd[2];
    fSliceNumInd[0] = (fSliceNum0 + 0.5) / (float) EPIPOLAR_SLICE_NUM;
    fSliceNumInd[1] = frac(fSliceNumInd[0] + 1 / (float) EPIPOLAR_SLICE_NUM);

    float fSliceWeight[2];
    fSliceWeight[1] = fSliceNum - fSliceNum0;
    fSliceWeight[0] = 1 - fSliceWeight[1];

    float fTotalWeight = 0;
    float3 f3TotalInscatter = 0;
    [unroll]
    for (uint i = 0; i < 2; i++)
    {
        float4 f4SliceStartEnd = g_tex2DSliceEnd.SampleLevel(samLinearClamp, float2(fSliceNumInd[i], 0.5), 0);
        float2 f2SliceRayDir = f4SliceStartEnd.zw - f4SliceStartEnd.xy;
        float fSliceRayLenSqr = dot(f2SliceRayDir, f2SliceRayDir);//instead of len
        
        float fEpipolarInSlice = dot(In.m_f2PosPS - f4SliceStartEnd.xy, f2SliceRayDir) / max(fSliceRayLenSqr, 1e-8);
        float fSampleNum = fEpipolarInSlice * (float)(EPIPOLAR_SAMPLE_NUM - 1);
        float fPreSampleNum = floor(fSampleNum);
        float fSampleWeight = fSampleNum - fPreSampleNum;
        
        float fPreSampleNumInd = (fPreSampleNum + 0.5) / (float)EPIPOLAR_SAMPLE_NUM;
        float2 f2PreEpipolarUV = float2(fPreSampleNumInd, fSliceNumInd[i]);
        
        float2 f2EpipolarDim = float2(EPIPOLAR_SAMPLE_NUM, EPIPOLAR_SLICE_NUM);
        float2 f2LocalCamDepth = g_tex2DEpipolarSampleCamDepth.Gather(samLinearClamp,f2PreEpipolarUV + float2(0.5, 0.5f) / f2EpipolarDim.xy).wz;
        float2 f2MaxDepth = max(f2LocalCamDepth, max(fCamDepth, 1));
        float2 f2DepthWeight = m_fRefinementThreshold * 0.2 / max(abs(fCamDepth - f2LocalCamDepth) / f2MaxDepth, m_fRefinementThreshold * 0.2);
        f2DepthWeight = pow(saturate(f2DepthWeight), 4);
        
        //bool2 b2IsIntersectEarth = bool2(ComputeIsIntersectEarth(g_tex2DEpipolarSample.SampleLevel(samPointClamp, f2PreEpipolarUV, 0, int2(0, 0))),
        //                                  ComputeIsIntersectEarth(g_tex2DEpipolarSample.SampleLevel(samPointClamp, f2PreEpipolarUV, 0, int2(0, 1))));
        //bool2 b2IntersectWeight = bIsIntersectEarth == b2IsIntersectEarth;
        
        //float3 f3Insctr0 = g_tex2DUnshadowedSampleScatter.SampleLevel(samPointClamp, f2PreEpipolarUV, 0, int2(0, 0));
        //float3 f3Insctr1 = g_tex2DUnshadowedSampleScatter.SampleLevel(samPointClamp, f2PreEpipolarUV, 0, int2(0, 1));
        //float3 f3MaxInsctr = max(f3Insctr0, f3Insctr1);
        
        //float fAveLogLum = 0.1;
        //const float middleGray = 1.03 - 2 / (2 + log10(fAveLogLum + 1));
        //float3 f3MinInsctrThreshold = (0.02f * fAveLogLum.xxx / RGB_TO_LUMINANCE.xyz) / middleGray;
        
        //f3MaxInsctr = max(f3MaxInsctr, f3MinInsctrThreshold);
        //bool bFlag = all(abs(f3Insctr0 - f3Insctr1) / f3MaxInsctr < m_fRefinementThreshold);
        //if(!bFlag)
        //    discard;
        
        float2 f2BilateralWeight = f2DepthWeight * float2(1 - fSampleWeight, fSampleWeight) * fSliceWeight[i];
        //float2 f2BilateralWeight = b2IntersectWeight * f2DepthWeight * float2(1 - fSampleWeight, fSampleWeight) * fSliceWeight[i];
        f2BilateralWeight *= (abs(fEpipolarInSlice - 0.5) < 0.5 + 0.5 * f2EpipolarDim.x);

        float fCurrTotalWeight = f2BilateralWeight.x + f2BilateralWeight.y;
        float fOffsetU = f2BilateralWeight.y / max(fCurrTotalWeight, 0.001) / f2EpipolarDim.x;
        float2 f2SampleUV = f2PreEpipolarUV + float2(fOffsetU, 0);
        
        f3TotalInscatter += fCurrTotalWeight * g_tex2DInterpolatedScatter.SampleLevel(samLinearClamp, f2SampleUV, 0);
        //f3TotalInscatter += float3(1,0,0);
        fTotalWeight += fCurrTotalWeight;
    }
    if (fTotalWeight < 1e-2)
    {
        //// Discarded pixels will keep 0 value in stencil and will be later
        //// processed to correct scattering
        //return float4(1, 0, 0, 1);
        discard;
    }
    float3 f3Inscatter = f3TotalInscatter / fTotalWeight;

    float3 f3BackColor = g_tex2DColorBuffer.SampleLevel(samLinearClamp, f2UV, 0);
    //if (misc.fIsLightInSpaceCorrect)
        f3BackColor *= (fCamDepth > camera.fFarZ) ? m_f4ExtraLight : 1;

    float4 f4PosInWorld = mul(float4(In.m_f2PosPS, fDepth,1), camera.InvViewProj);
    f4PosInWorld /= f4PosInWorld.w;
    float3 f3RayInWorld = f4PosInWorld.xyz - camera.f3CameraPos;
    float f3RayLenInWorld = length(f3RayInWorld);
    f3RayInWorld /= max(f3RayLenInWorld, 1e-8);
    
    float3 f3Transmittance = f3RayLenInWorld == 0 ? 1:GetTransmittanceIntegralAnalytic(fHeight, fCosZenithAngle, f3RayLenInWorld);
    f3BackColor *= f3Transmittance;
    
    float3 f3Irradiance = 0;
    if(fCamDepth>camera.fFarZ && !bIsIntersectEarth)
    {
        float fCosSunAngle = dot(f3CameraPos, light.f3LightDir) / fHeight;
        float fSunViewAngle = dot(f3RayInWorld, light.f3LightDir);
        f3Irradiance = GetIrradianceDirectFromSun(fHeight, fCosSunAngle, fSunViewAngle);
    }

    return float4(ToneMap(f3BackColor + f3Inscatter + f3Irradiance), 1);

}

technique11 ApplyInterpolateScatterTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest_IncrStencil, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ApplyInterpolateScatter()));
    }
}


float4 FixInterpolateScatter(QuadVertexOut In) : SV_Target
{
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float fCamDepth = g_tex2DSpaceLinearDepth.SampleLevel(samLinearClamp, f2UV, 0);
    float fDepth = camera.Proj[2][2] + camera.Proj[3][2] / fCamDepth;

    bool bIsIntersectEarth = false;
    float4 f4RayEndInWorldSpace = mul(float4(In.m_f2PosPS, fDepth, 1.f), camera.InvViewProj);
    float3 f3ViewRay = f4RayEndInWorldSpace.xyz / f4RayEndInWorldSpace.w - camera.f3CameraPos;
    float fViewRayLen = length(f3ViewRay);
    f3ViewRay /= fViewRayLen;
    float3 f3EarthCenter = float3(0, -atmosphere.bottom_radius, 0);
    float3 f3CameraPos = camera.f3CameraPos - f3EarthCenter;
    float fHeight = length(f3CameraPos);
    float fCosZenithAngle = dot(f3CameraPos, f3ViewRay) / fHeight;
    float fRMu = fHeight * fCosZenithAngle;
    float fIntersectEarthDiscriminant = fRMu * fRMu - fHeight * fHeight + atmosphere.bottom_radius * atmosphere.bottom_radius;
    if (fIntersectEarthDiscriminant >= 0 && fCamDepth > camera.fFarZ && fCosZenithAngle < 0)
        bIsIntersectEarth = true;

    float3 f3Inscatter = 0;
    if(misc.fEnableLightShaft)
    {
        f3Inscatter = ComputeShadowInscatter(In.m_f2PosPS, fCamDepth, fDepth,false, 0).xyz;
    }
    else
    {
        f3Inscatter = ComputeUnshadowedSampleScatter(In.m_f2PosPS, fCamDepth,fDepth).xyz;
    }

    float3 f3BackColor = g_tex2DColorBuffer.SampleLevel(samLinearClamp, f2UV, 0);
    //if (misc.fIsLightInSpaceCorrect)
        f3BackColor *= (fCamDepth > camera.fFarZ) ? m_f4ExtraLight : 1;

    float4 f4PosInWorld = mul(float4(In.m_f2PosPS, fDepth, 1), camera.InvViewProj);
    f4PosInWorld /= f4PosInWorld.w;
    float3 f3RayInWorld = f4PosInWorld.xyz - camera.f3CameraPos;
    float f3RayLenInWorld = length(f3RayInWorld);
    f3RayInWorld /= max(f3RayLenInWorld, 1e-8);

    float3 f3Transmittance = f3RayLenInWorld == 0 ? 1 : GetTransmittanceIntegralAnalytic(fHeight, fCosZenithAngle, f3RayLenInWorld);
    f3BackColor *= f3Transmittance;

    float3 f3Irradiance = 0;
    if (fCamDepth > camera.fFarZ && !bIsIntersectEarth)
    {
        float fCosSunAngle = dot(f3CameraPos, light.f3LightDir) / fHeight;
        float fSunViewAngle = dot(f3RayInWorld, light.f3LightDir);
        f3Irradiance = GetIrradianceDirectFromSun(fHeight, fCosSunAngle, fSunViewAngle);
    }

    return float4(ToneMap(f3BackColor + f3Inscatter + f3Irradiance), 1);

}

technique11 FixInterpolateScatterTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest_StEqual_KeepStencil, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, FixInterpolateScatter()));
    }
}

static const float fSunAngularRadius = 32.f / 2.f / 60.f * ((2.f * PI) / 180); // Sun angular DIAMETER is 32 arc minutes
static const float fTanSunAngularRadius = tan(fSunAngularRadius);

QuadVertexOut RenderSunVS(in uint VertexId : SV_VertexID,
                          in uint InstID : SV_InstanceID)
{
    float2 fCotanHalfFOV = float2(camera.Proj[0][0], camera.Proj[1][1]);
    float2 f2SunScreenPos = light.f4LightScreenPos.xy;
    float2 f2SunScreenSize = fTanSunAngularRadius * fCotanHalfFOV;
    float4 MinMaxUV = f2SunScreenPos.xyxy + float4(-1, -1, 1, 1) * f2SunScreenSize.xyxy;
 
    QuadVertexOut Verts[4] =
    {
        { float4(MinMaxUV.xy, 1.0, 1.0), MinMaxUV.xy, InstID },
        { float4(MinMaxUV.xw, 1.0, 1.0), MinMaxUV.xw, InstID },
        { float4(MinMaxUV.zy, 1.0, 1.0), MinMaxUV.zy, InstID },
        { float4(MinMaxUV.zw, 1.0, 1.0), MinMaxUV.zw, InstID }
    };

    return Verts[VertexId];
}

float4 RenderSunPS(QuadVertexOut In) : SV_Target
{
    float2 fCotanHalfFOV = float2(camera.Proj[0][0], camera.Proj[1][1]);
    float2 f2SunScreenSize = fTanSunAngularRadius * fCotanHalfFOV;
    float2 f2dXY = (In.m_f2PosPS - light.f4LightScreenPos.xy) / f2SunScreenSize;

    float4 f4RayEndInWorldSpace = mul(float4(In.m_f2PosPS, 1.f, 1.f), camera.InvViewProj);
    float3 f3ViewRay = f4RayEndInWorldSpace.xyz / f4RayEndInWorldSpace.w - camera.f3CameraPos;
    float fViewRayLen = length(f3ViewRay);
    f3ViewRay /= fViewRayLen;
    float3 f3EarthCenter = float3(0, -atmosphere.bottom_radius, 0);
    float3 f3CameraPos = camera.f3CameraPos - f3EarthCenter;
    float fHeight = length(f3CameraPos);
    float fCosZenithAngle = dot(f3CameraPos, f3ViewRay) / fHeight;
    float fRMu = fHeight * fCosZenithAngle;
    float fIntersectEarthDiscriminant = fRMu * fRMu - fHeight * fHeight + atmosphere.bottom_radius * atmosphere.bottom_radius;
    if (fIntersectEarthDiscriminant >= 0 && fCosZenithAngle < 0)
        discard;
    return sqrt(saturate(1 - dot(f2dXY, f2dXY)));
}

technique11 RenderSunTech
{
    pass
    {
        SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_EnableDepthEqTest, 0);

        SetVertexShader(CompileShader(vs_5_0, RenderSunVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, RenderSunPS()));
    }
}