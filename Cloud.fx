#include "Common.fx"

float3 hash(float3 p) // replace this by something better. really. do
{
    p = float3(dot(p, float3(127.1, 311.7, 74.7)),
			  dot(p, float3(269.5, 183.3, 246.1)),
			  dot(p, float3(113.5, 271.9, 124.6)));

    return -1.0 + 2.0 * frac(sin(p) * 43758.5453123);
}

// return value noise (in x) and its derivatives (in yzw)
float4 Perlin(in float3 x)
{
    // grid
    float3 p = floor(x);
    float3 w = frac(x);
    
#if 1
    // quintic interpolant
    float3 u = w * w * w * (w * (w * 6.0 - 15.0) + 10.0);
    float3 du = 30.0 * w * w * (w * (w - 2.0) + 1.0);
#else
    // cubic interpolant
    float3 u = w*w*(3.0-2.0*w);
    float3 du = 6.0*w*(1.0-w);
#endif    
    
    // gradients
    float3 ga = hash(p + float3(0.0, 0.0, 0.0));
    float3 gb = hash(p + float3(1.0, 0.0, 0.0));
    float3 gc = hash(p + float3(0.0, 1.0, 0.0));
    float3 gd = hash(p + float3(1.0, 1.0, 0.0));
    float3 ge = hash(p + float3(0.0, 0.0, 1.0));
    float3 gf = hash(p + float3(1.0, 0.0, 1.0));
    float3 gg = hash(p + float3(0.0, 1.0, 1.0));
    float3 gh = hash(p + float3(1.0, 1.0, 1.0));
    
    // projections
    float va = dot(ga, w - float3(0.0, 0.0, 0.0));
    float vb = dot(gb, w - float3(1.0, 0.0, 0.0));
    float vc = dot(gc, w - float3(0.0, 1.0, 0.0));
    float vd = dot(gd, w - float3(1.0, 1.0, 0.0));
    float ve = dot(ge, w - float3(0.0, 0.0, 1.0));
    float vf = dot(gf, w - float3(1.0, 0.0, 1.0));
    float vg = dot(gg, w - float3(0.0, 1.0, 1.0));
    float vh = dot(gh, w - float3(1.0, 1.0, 1.0));
	
    // interpolations
    return float4(va + u.x * (vb - va) + u.y * (vc - va) + u.z * (ve - va) + u.x * u.y * (va - vb - vc + vd) + u.y * u.z * (va - vc - ve + vg) + u.z * u.x * (va - vb - ve + vf) + (-va + vb + vc - vd + ve - vf - vg + vh) * u.x * u.y * u.z, // value
                 ga + u.x * (gb - ga) + u.y * (gc - ga) + u.z * (ge - ga) + u.x * u.y * (ga - gb - gc + gd) + u.y * u.z * (ga - gc - ge + gg) + u.z * u.x * (ga - gb - ge + gf) + (-ga + gb + gc - gd + ge - gf - gg + gh) * u.x * u.y * u.z + // derivatives
                 du * (float3(vb, vc, ve) - va + u.yzx * float3(va - vb - vc + vd, va - vc - ve + vg, va - vb - ve + vf) + u.zxy * float3(va - vb - ve + vf, va - vb - vc + vd, va - vc - ve + vg) + u.yzx * u.zxy * (-va + vb + vc - vd + ve - vf - vg + vh)));
}

float rand3(float3 co)
{
    return frac(sin(dot(co.xyz, float3(12.9898, 78.233, 42.1897))) * 43758.5453);
}

float3 SphericalRandom(float3 co)
{
    float r = 1.f - pow(rand3(co), 4.0);
    float az = rand3(43.1138 * co) * 2 * PI;
    float sine_el = rand3(17.981 * co) * 2.0 - 1.0;
    float el = asin(sine_el);
    float cos_el = cos(el);
    float3 v;
    v.x = r * cos(az) * cos_el;
    v.y = r * sine_el;
    v.z = r * sin(az) * cos_el;
    return v;
}

float Worley(float3 uvw, int grid, int seed)
{
    float n = float(grid);
    float3 pos = n * uvw;
    float3 fraction, intpart;
    fraction = modf(pos, intpart);
    int3 integer = int3(intpart);
    float3 loc = fraction - 0.5;
    for (int i = -1; i < 2; ++i)
    {
        for (int j = -1; j < 2; ++j)
        {
            for (int k = -1; k < 2; ++k)
            {
                float3 c = float3(i, j, k);
                int3 u = integer + int3(i, j, k);
                u %= grid;
                float3 random = SphericalRandom(u * seed);
                c += random.xyz;
                float k1 = length(loc - c);
                n = min(n, k1);
            }
        }
    }
    float p = saturate(1 - pow(1.5 * n, 4.0)) / 1.0;
    return p;
}

float PerlinfBm(float3 uvw, int octaves, float persistence)
{
    float result = 0;
    float total = 0;
    float amplitude = 1;
   
    for (int i = 0; i < octaves; ++i)
    {
        uvw /= amplitude;
        result += amplitude * Perlin(uvw).x;
        total += amplitude;
        amplitude *= persistence;
    }
    result /= total;
    return result;
}

float WorleyfBm(float3 uvw, int octaves, float persistence, float grid, float seed)
{
    float result = 0;
    float total = 0;
    float amplitude = 1;
    for (int i = 0; i < octaves; ++i)
    {
        uvw /= amplitude;
        result += amplitude * Worley(uvw, grid, seed);
        total += amplitude;
        amplitude *= persistence;
    }
    result /= total;
    return result;
}

float4 ComputePerlinWorleyNoiseTexture(QuadVertexOut In) : SV_Target
{
    float3 f3UVW = float3(ProjToUV(In.m_f2PosPS), misc.f2WQ.x);
    float4 result;
    
    result.x = (saturate(PerlinfBm(f3UVW * 32, 4, 0.5).x * 0.5 + 0.5) + WorleyfBm(f3UVW, 4, 0.5, 8, 1)) / 2;
    result.y = WorleyfBm(f3UVW, 4, 0.5, 8, 2);
    result.z = WorleyfBm(f3UVW, 4, 0.5, 16, 1);
    result.w = WorleyfBm(f3UVW, 4, 0.5, 32, 1);

    return result;
}

technique11 ComputePerlinWorleyNoiseTex3DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0, 0, 0, 0), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputePerlinWorleyNoiseTexture()));
    }
}

float4 ComputeWorleyNoiseTexture(QuadVertexOut In) : SV_Target
{
    float3 f3UVW = float3(ProjToUV(In.m_f2PosPS), misc.f2WQ.x);
    float4 result;
    
    result.x = WorleyfBm(f3UVW, 4, 0.5, 5, 1);
    result.y = WorleyfBm(f3UVW, 4, 0.5, 9, 2);
    result.z = WorleyfBm(f3UVW, 4, 0.5, 16, 3);
    result.w = WorleyfBm(f3UVW, 4, 0.5, 24, 4);
    //result.x = Worley(f3UVW, 5, 1);
    //result.y = Worley(f3UVW, 9, 2);
    //result.z = Worley(f3UVW, 16, 3);
    //result.w = Worley(f3UVW, 24, 4);

    return result;
}

technique11 ComputeWorleyNoiseTex3DTech
{
    pass P0
    {
        SetBlendState(NoBlending, float4(0, 0, 0, 0), 0xFFFFFFFF);
        SetRasterizerState(RS_SolidFill_NoCull);
        SetDepthStencilState(DSS_NoDepthTest, 0);

        SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
        SetGeometryShader(NULL);
        SetPixelShader(CompileShader(ps_5_0, ComputeWorleyNoiseTexture()));
    }
}

struct CloudTypeLayer
{
    float2 f2LayerHeightScale;
    float2 f2LayerDensityPoint;
};

struct CloudParams
{
    CloudTypeLayer mCloudTypeLayer[3];

    float2 f2CloudLayerHeightScale;
    float fTransition;
    float fUpperDensity;
};

cbuffer cbCloudParams
{
    CloudParams cloud;
};

float SampleNoise(Texture3D<float4> tex3DNoiseTex,float3 f3UVW,float fMipLevel)
{
    float4 f4Noise = tex3DNoiseTex.SampleLevel(samLinearClamp, f3UVW, fMipLevel);
    float fLowFreqFBM = f4Noise.y * 0.625 + f4Noise.z * 0.25 + f4Noise.w * 0.125;
    float fNoise = ReMap(f4Noise.x, fLowFreqFBM - 1, 1.0, 0.0, 1.0); //[-1,1] - [0,1]  
    return fNoise;
}

float HGPhaseFunction(float w,float g)
{
    float g2 = g * g;
    return 1 / (4 * PI) * (1 - g2) / pow(max((1 + g2 - 2 * g * cos(w)), 1e-20), 1.5);
}

float GetHumidity(float fHeightScale)
{
    float fInterpolate = pow(saturate((fHeightScale - cloud.f2CloudLayerHeightScale.x) / cloud.fTransition), 0.5);
    float fHumdity = lerp(1, cloud.fUpperDensity, fInterpolate); 
    fHumdity *= pow(saturate((1.0 - fHeightScale) / cloud.fTransition), 0.5); 
    return fHumdity;
}

float GetCloudHeightGradientType(float fHeightScale,int type)
{
    //if(0 == type)
    //{
    //    return ReMap(fHeightScale, 0.0, 0.1, 0.0, 1.0) * ReMap(fHeightScale, 0.2, 0.3, 1.0, 0.0); // Stratus
    //}
    //else if(1 == type)
    //{
    //    return ReMap(fHeightScale, 0.0, 0.3, 0.0, 1.0) * ReMap(fHeightScale, 0.4, 0.6, 1.0, 0.0); // Stratuscumulus
    //}
    //else
    //{
    //    return ReMap(fHeightScale, 0.5, 0.7, 0.0, 1.0) * ReMap(fHeightScale, 0.9, 1.0, 1.0, 0.0); // Cumulus
    //}
    return ReMap(fHeightScale, 
                cloud.mCloudTypeLayer[type].f2LayerHeightScale.x, 
                cloud.mCloudTypeLayer[type].f2LayerDensityPoint.x, 0.0, 1.0) *
            ReMap(fHeightScale,
                cloud.mCloudTypeLayer[type].f2LayerHeightScale.y, 
                cloud.mCloudTypeLayer[type].f2LayerDensityPoint.y, 1.0, 0.0);

}

float3 GetNoiseUVW(float3 f3Pos, float fHeightScale)
{
    //float2 f2XZ = f3Pos.xz - camera.f3CameraPos.xz;
    //float fRadius = length(f3Pos.xz - camera.f3CameraPos.xz);
    //f2XZ /= fRadius;
    float2 f2XZ = f3Pos.xz;
    f2XZ = saturate(f2XZ * 0.5 + 0.5);
	return float3(f2XZ.x, fHeightScale, f2XZ.y);
}

float GetBaseCloudDensity(float3 f3Pos, float fHeight, float fHeightScale,float fMipLevel,out float fType)
{    
    float fCoverage = 0.8; // cloudmap.r
    fType = 0;//cloumap.b

    float3 f3UVW = GetNoiseUVW(f3Pos, fHeightScale);
    //float fNoise = SampleNoise(g_tex3DPerlinWorleyNoise, f3UVW, fMipLevel);
    float4 f4Noise = g_tex3DPerlinWorleyNoise.SampleLevel(samLinearClamp, f3UVW, fMipLevel);
    float fLowFreqFBM = f4Noise.y * 0.625 + f4Noise.z * 0.25 + f4Noise.w * 0.125;
    float fNoise = ReMap(f4Noise.x, fLowFreqFBM - 1, 1.0, 0.0, 1.0); //[-1,1] - [0,1]  
    //float fHumidity = GetHumidity(fHeightScale);
    //float fDiffusivity = 1;
    //float fDensity = saturate((fNoise + fHumidity - 1) / fDiffusivity);
    
    float fBaseCloud = fNoise * GetCloudHeightGradientType(fHeightScale, fType);
    return fBaseCloud;
    float fBaseCloudWithCoverage = ReMap(fBaseCloud, fCoverage, 1.0, 0.0, 1.0);
    fBaseCloudWithCoverage *= fCoverage;

    return fBaseCloudWithCoverage;
}

float GetFullCloudDensity(float fBaseCloudDensity, float3 f3Pos, float fMipLevel)
{
    float fFullCloudDensity = fBaseCloudDensity;
    return fFullCloudDensity;

}

float GetCloudTransmittance(float fSampleDensity)
{
    float fRainAbsorption = 1; // cloudmap.g
    float fBeerLaw = exp(-fRainAbsorption * fSampleDensity);
    float fPowerEffect = 1 - exp(-2 * fSampleDensity);
    float fTransmittance = 2 * fBeerLaw * fPowerEffect;
    return fTransmittance;
}

void GetHeightAndScale(float3 f3Pos, out float fHeight, out float fHeightScale)
{
    fHeight = length(f3Pos) - atmosphere.bottom_radius;
    fHeightScale = fHeight * atmosphere.radius_scale;
}

float GetCloudDensityToLight(float3 f3Pos, float fMipLevel)
{
    float fDensity = 0;
    int SAMPLE_COUNT = 6;
    float fHeight;
    float fHeightScale;
    float fType;

    for (int i = 0; i < 6;i++)
    {
        float f3Step = light.f3LightDir * SphericalRandom(f3Pos);
        f3Pos = f3Pos + f3Step;
        GetHeightAndScale(f3Pos, fHeight, fHeightScale);
        float fBaseDensity = GetBaseCloudDensity(f3Pos, fHeight, fHeightScale, fMipLevel + 1, fType);
        if(fDensity<0.3)
        {
            fDensity += GetFullCloudDensity(fBaseDensity, f3Pos, fMipLevel + 1);
        }
        else
        {
            fDensity += fBaseDensity;
        }
    }
    return saturate(fDensity);

}

float4 DrawCloud(QuadVertexOut In) : SV_Target
{
    uint2 ui2XY = In.m_f4Pos.xy;
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float fCamDepth = g_tex2DSpaceDepth.Load(uint3(ui2XY, 0));
    float fDepth = camera.Proj[2][2] + camera.Proj[3][2] / fCamDepth;

    float3 f3BackColor = g_tex2DColorBuffer.SampleLevel(samLinearClamp, f2UV, 0);

    float4 f4RayEndPos = mul(float4(In.m_f2PosPS, fDepth, 1), camera.InvViewProj);
    //float4 f4RayEndPos = mul(float4(In.m_f2PosPS, 1, 1), camera.InvViewProj);
    float3 f3ViewRay = f4RayEndPos.xyz / f4RayEndPos.w - camera.f3CameraPos;
    float fViewRayLen = length(f3ViewRay);
    f3ViewRay /= fViewRayLen;
    
    float3 f3EarthCenter = float3(0, -atmosphere.bottom_radius, 0);
    float3 f3CameraPos = camera.f3CameraPos - f3EarthCenter;
    float4 f4RMuMuSNu;
    f4RMuMuSNu.xyz = GetRMuMuS(f3CameraPos, f3ViewRay);
    f4RMuMuSNu.w = dot(f3ViewRay, light.f3LightDir);

    //bool bIsNoScatter, bIsMarchToAtmosphere, bIsMarchToEarth, bIsIntersectEarth;
    //float fDistToAtmosphereNear, fDistToAtmosphereFar, fDistToEarthNear;
    //float fViewRayLenInWorldSpace = GetRayMarchLen(f4RMuMuSNu, fCamDepth, fViewRayLen,
    //                                        bIsNoScatter, bIsMarchToAtmosphere, bIsMarchToEarth, bIsIntersectEarth,
    //                                        fDistToAtmosphereNear, fDistToAtmosphereFar, fDistToEarthNear);
    //float3 f3StartPos = f3CameraPos;
    //float3 f3EndPos = f3CameraPos + fViewRayLenInWorldSpace * f3ViewRay;
    //if (fDistToAtmosphereNear > 0.f)
    //{
    //    float3 f3DistToAtmosphereNear = fDistToAtmosphereNear * f3ViewRay;
    //    f3CameraPos += f3DistToAtmosphereNear;
    //    f3StartPos += f3DistToAtmosphereNear;
    //    fViewRayLenInWorldSpace -= fDistToAtmosphereNear;
    //    fDistToAtmosphereFar -= fDistToAtmosphereNear;
    //    if (bIsMarchToEarth)
    //        fDistToEarthNear -= fDistToAtmosphereNear;
    //    f4RMuMuSNu.xyz = GetRMuMuS(f3CameraPos, f3ViewRay);
    //}

    float fRMu = f4RMuMuSNu.x * f4RMuMuSNu.y;
    float fRSqr = f4RMuMuSNu.x * f4RMuMuSNu.x;
    float fRMuSqr = fRMu * fRMu;

    float2 f2LayerHeight = atmosphere.bottom_radius + cloud.f2CloudLayerHeightScale / atmosphere.radius_scale;
    float2 f2IntersectLayerDiscriminant = fRMuSqr - fRSqr + f2LayerHeight * f2LayerHeight;

    float fHorizonMu = -sqrt((fRSqr - atmosphere.bottom_radius * atmosphere.bottom_radius) / fRSqr);

    float fStartLen, fEndLen;
    float2 f2IntersectLayerDiscriminantSqrt = float2(SafeSqrt(f2IntersectLayerDiscriminant.x), SafeSqrt(f2IntersectLayerDiscriminant.y));
    float2 f2IntersectLayerFar = -fRMu + f2IntersectLayerDiscriminantSqrt;
    float2 f2IntersectLayerNear = -fRMu - f2IntersectLayerDiscriminantSqrt;
    if (f4RMuMuSNu.x >= f2LayerHeight.x && f4RMuMuSNu.x <= f2LayerHeight.y)
    {
        fStartLen = 0;
    }
    else if (f4RMuMuSNu.x > f2LayerHeight.y)
    {
        if (f2IntersectLayerDiscriminant.y < 0 ||
		(f2IntersectLayerDiscriminant.y >= 0 && f4RMuMuSNu.y > 0))
        {
            return float4(f3BackColor, 1);
        }
        fStartLen = f2IntersectLayerNear.x;
    }
    else
    {
        if (f4RMuMuSNu.y <= fHorizonMu)
            return float4(f3BackColor, 1);
        fStartLen = f2IntersectLayerFar.x;
    }

    if (f4RMuMuSNu.y <= fHorizonMu)
    {
        fEndLen = f2IntersectLayerNear.x;
    }
    else
    {
        fEndLen = f2IntersectLayerFar.y;
    }

    if (fEndLen <= fStartLen)
        return float4(1, 0, 0, 1);
    float3 f3StartPos = f3CameraPos + f3ViewRay * fStartLen;
    float3 f3EndPos = f3CameraPos + f3ViewRay * fEndLen;
    float fViewRayLenInWorldSpace = fEndLen - fStartLen;

    int SAMPLE_COUNT = ((f4RMuMuSNu.x > fHorizonMu) && (f4RMuMuSNu.x < 0.2 + fHorizonMu)) ? 128 : 64;
    float fStepLen = fViewRayLenInWorldSpace / SAMPLE_COUNT;
    float3 f3Pos = f3StartPos;
    float3 f3Step = f3ViewRay * fStepLen;
    
    float fDensityToCam = 0;
    float fCloudTest = 0;
    float fMipLevel = 0;
    float3 f3TotalInscatter = 0;
    float3 f3PreInscatter = 0;
    float3 f3CurInscatter = 0;
    int fZeroSampleCount = 0;
    float fType;
    float fHeight;
    float fHeightScale;
    
    //return float4(fStartLen, fEndLen, 0, 1);
    for (int i = 0; i < SAMPLE_COUNT; ++i)
    {
        GetHeightAndScale(f3Pos, fHeight, fHeightScale);
        if (fHeightScale < cloud.f2CloudLayerHeightScale.x)
        {
            f3Pos += f3Step;
            continue;
        }
        if (fCloudTest > 0.f)
        {
            float fSampleDensity = GetFullCloudDensity(fCloudTest, f3Pos, fMipLevel);
        
            fDensityToCam += fSampleDensity;
		//float fDensityToLight = GetCloudDensityToLight(f3Pos, fMipLevel);
		//float fTransmittance = GetCloudTransmittance(fDensityToCam) * GetCloudTransmittance(fDensityToLight);
            float fTransmittance = GetCloudTransmittance(fDensityToCam);
		//f3CurInscatter = fTransmittance;
		//f3PreInscatter = f3CurInscatter;
            f3TotalInscatter += fTransmittance;
        
            fCloudTest = 0.f;
            f3Pos += f3Step;
        }
        else
        {
            fCloudTest = GetBaseCloudDensity(f3Pos, fHeight, fHeightScale, fMipLevel, fType);
            if (fCloudTest == 0.f)
            {
                f3Pos += f3Step;
            }
        }
    }
    f3TotalInscatter *= HGPhaseFunction(f4RMuMuSNu.w, atmosphere.mie_g) * fStepLen;

    if (fCamDepth > camera.fFarZ)
    {
        return float4(f3TotalInscatter, 1);
    }
    else
    {
        return float4(f3BackColor, 1); 
    }
}


technique11 DrawCloudTech
{
	pass
	{
		SetBlendState(NoBlending, float4(0.0f, 0.0f, 0.0f, 0.0f), 0xFFFFFFFF);
		SetRasterizerState(RS_SolidFill_NoCull);
		SetDepthStencilState(DSS_NoDepthTest, 0);

		SetVertexShader(CompileShader(vs_5_0, GenerateScreenSizeQuadVS()));
		SetGeometryShader(NULL);
		SetPixelShader(CompileShader(ps_5_0, DrawCloud()));
	}
}