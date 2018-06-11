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

    float2 f2CloudLayerHeight;
    float fScale;
    float padding;
};

cbuffer cbCloudParams
{
    CloudParams cloud;
};

float SampleNoise(float3 f3Pos, Texture2D<float4> tex)
{
    float2 f2Dim;
    tex.GetDimensions(f2Dim.x, f2Dim.y);
    float fSize = 1.f / f2Dim.y;
    
    float3 f3UVW = f3Pos * fSize;
    f3UVW = f3UVW - floor(f3UVW);
    
    float fU = f3UVW.x / fSize - 0.5; // unit
    float fStrench = floor(fU);
    float fInterpolate = fU - fStrench;
    
    fU = (fStrench + f3UVW.z) * fSize; // (fStrench / fSize + f3UVW.z / fSize) 
    
    float2 fNoise = float2(tex.SampleLevel(samLinearWrap, float2(fU, f3UVW.y), 0).x,
					  tex.SampleLevel(samLinearWrap, float2(fU + fSize, f3UVW.y), 0).x);
    
    return lerp(fNoise.x, fNoise.y, fInterpolate);
}

float HGPhaseFunction(float w,float g)
{
    float g2 = g * g;
    return 1 / (4 * PI) * (1 - g2) / pow(max((1 + g2 - 2 * g * cos(w)), 1e-20), 1.5);
}

//float GetHumidity(float fHeightScale)
//{
//    float fInterpolate = pow(saturate(fHeightScale / cloud.fTransition), 0.5);
//    float fHumdity = lerp(1, cloud.fUpperDensity, fInterpolate); 
//    fHumdity *= pow(saturate((1.0 - fHeightScale) / cloud.fTransition), 0.5); 
//    return fHumdity;
//}

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
                cloud.mCloudTypeLayer[type].f2LayerDensityPoint.y,
                cloud.mCloudTypeLayer[type].f2LayerHeightScale.y, 1.0, 0.0);

}

float GetCloudDensity(float3 f3Pos, float fHeightScale, float fCoverage, float fType, bool bCheap)
{
    float3 f3WindDir = float3(1.f, 0.f, 0.f);
    float fCloudSpeed = 1.0;
    float fCloudTopOffset = 0.5;
    f3Pos += fHeightScale * f3WindDir * fCloudTopOffset;
    f3Pos += (f3WindDir + float3(0, 0.1, 0)) * fCloudSpeed * misc.fTime;

    float fNoise = SampleNoise(f3Pos, g_tex2DNoiseBasePacked);
    float fBaseCloud = fNoise * GetCloudHeightGradientType(fHeightScale, fType);
    float fBaseCloudDensity = ReMap(fBaseCloud, fCoverage, 1.0, 0.0, 1.0);
    fBaseCloudDensity *= fCoverage;

    if (fBaseCloudDensity == 0)
        return 0;

    if (!bCheap)
    {
        float fHightFreqNoise = SampleNoise(f3Pos * 0.1, g_tex2DNoiseDetailPacked).x;
	// Transition from wispy shapes to billowy shapes over height
        float fHightFreqNoiseModifier = lerp(fHightFreqNoise, 1 - fHightFreqNoise, saturate(10 * fHeightScale));
    
        float fFullCloudDensity = ReMap(fBaseCloudDensity, fHightFreqNoiseModifier * 0.2, 1.0, 0.0, 1.0);
        return fFullCloudDensity;
    }
    else
    {
        return fBaseCloudDensity;

    }

}

float GetFullCloudDensity(float3 f3Pos, float fBaseCloudDensity, float fHeightScale, float fCoverage, float fType)
{
    float fHightFreqNoise = SampleNoise(f3Pos * 0.1, g_tex2DNoiseDetailPacked).x;
	// Transition from wispy shapes to billowy shapes over height
    float fHightFreqNoiseModifier = lerp(fHightFreqNoise, 1 - fHightFreqNoise, saturate(10 * fHeightScale));
    
    float fFullCloudDensity = ReMap(fBaseCloudDensity, fHightFreqNoiseModifier * 0.2, 1.0, 0.0, 1.0);

    return fFullCloudDensity;
}


float GetCloudTransmittance(float fSampleDensity, float fRainAbsorption)
{
    float fBeerLaw = exp(-fRainAbsorption * fSampleDensity);
    float fPowerEffect = 1 - exp(-2 * fSampleDensity);
    float fTransmittance = 2 * fBeerLaw * fPowerEffect;
    return fTransmittance;
}

void GetHeightAndScale(float3 f3Pos, float3 f3EarthCenter ,out float fHeight, out float fHeightScale)
{
    fHeight = length(f3Pos - f3EarthCenter) - atmosphere.bottom_radius;
    fHeightScale = ReMap(fHeight, cloud.f2CloudLayerHeight.x, cloud.f2CloudLayerHeight.y, 0.0, 1.0);
}

void SampleWeatherTexture(float3 f3Pos, out float fCoverage, out float fRainAbsorption, out float fType)
{
    fCoverage = 0.9; // weatherMap.r
    fRainAbsorption = 1; //  weatherMap.g
    fType = 0; //  weatherMap.b
}

float GetCloudDensityToLight(float3 f3Pos, float3 f3EarthCenter)
{
    float fDensity = 0;
    int SAMPLE_COUNT = 6;
    float fHeight, fHeightScale;
    float fCoverage, fRainAbsorption, fType;

    float fRadius = length(f3Pos - f3EarthCenter);
    float fCosSunZenithAngle = dot(f3Pos, light.f3LightDir) / fRadius;

    float fRMu = fRadius * fCosSunZenithAngle;
    float fRSqr = fRadius * fRadius;
    float fRMuSqr = fRMu * fRMu;

    float fUpLayerHeight = atmosphere.bottom_radius + cloud.f2CloudLayerHeight.y;
    float fIntersectUpLayerDiscriminant = fRMuSqr - fRSqr + fUpLayerHeight * fUpLayerHeight;
    float fIntersectLayerFar = -fRMu + SafeSqrt(fIntersectUpLayerDiscriminant);
    float fLightStepLength = length(fIntersectLayerFar) / 6;
    if (fLightStepLength == 0)
        return 0;


    for (int i = 0; i < 6;i++)
    {
        //float3 f3Dir = float3(rand3(float3(0, i, 0)), rand3(float3(1, i, 0)), rand3(float3(0, i, 1)));
        //f3Dir = 2 * (light.f3LightDir * 2 + normalize(f3Dir));
        //float3 f3Step = f3Dir;
        
        //float3 f3Step = light.f3LightDir * SphericalRandom(f3Pos) * fLightStepLength * float(i);
        float3 f3Step = light.f3LightDir * float(i);

        f3Pos = f3Pos + f3Step;
        GetHeightAndScale(f3Pos, f3EarthCenter, fHeight, fHeightScale);
        SampleWeatherTexture(f3Pos, fCoverage, fRainAbsorption, fType);

        float fBaseDensity = GetCloudDensity(f3Pos, fHeightScale, fCoverage, fType,true);
        if(fDensity<0.3)
        {
            fDensity += GetCloudDensity(f3Pos, fHeightScale, fCoverage, fType,false);
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
    float2 f2UV = ProjToUV(In.m_f2PosPS);
    float fCamDepth = g_tex2DSpaceLinearDepth.SampleLevel(samLinearClamp, f2UV, 0);
    float fDepth = camera.Proj[2][2] + camera.Proj[3][2] / fCamDepth;

    float3 f3BackColor = g_tex2DColorBuffer.SampleLevel(samLinearClamp, f2UV, 0);

    float4 f4RayEndPos = mul(float4(In.m_f2PosPS, fDepth, 1), camera.InvViewProj);
    float3 f3ViewRay = f4RayEndPos.xyz / f4RayEndPos.w - camera.f3CameraPos;
    float fViewRayLen = length(f3ViewRay);
    f3ViewRay /= fViewRayLen;

    float3 f3EarthCenter = float3(0, -atmosphere.bottom_radius, 0);
    float3 f3CameraPos = camera.f3CameraPos - f3EarthCenter;
    float4 f4RMuMuSNu;
    f4RMuMuSNu.xyz = GetRMuMuS(f3CameraPos, f3ViewRay);
    f4RMuMuSNu.w = dot(f3ViewRay, light.f3LightDir);

    float fRMu = f4RMuMuSNu.x * f4RMuMuSNu.y;
    float fRSqr = f4RMuMuSNu.x * f4RMuMuSNu.x;
    float fRMuSqr = fRMu * fRMu;

    float2 f2LayerHeight = atmosphere.bottom_radius + cloud.f2CloudLayerHeight;
    float2 f2IntersectLayerDiscriminant = fRMuSqr - fRSqr + f2LayerHeight * f2LayerHeight;
  
    float fHorizonMu = -SafeSqrt((fRSqr - atmosphere.bottom_radius * atmosphere.bottom_radius) / fRSqr);

    float fStartLen, fEndLen;
    float2 f2IntersectLayerDiscriminantSqrt = float2(SafeSqrt(f2IntersectLayerDiscriminant.x),
                        SafeSqrt(f2IntersectLayerDiscriminant.y));
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
        fStartLen = f2IntersectLayerNear.y;
    }
    else
    {
        if (f4RMuMuSNu.y <= fHorizonMu || fHorizonMu == 0)
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

    float fViewRayLenInWorldSpace = fEndLen - fStartLen;
    fViewRayLenInWorldSpace = clamp(fViewRayLenInWorldSpace,0,10);
    fEndLen = fStartLen + fViewRayLenInWorldSpace;
    float3 f3StartPos = camera.f3CameraPos + f3ViewRay * fStartLen;
    float3 f3EndPos = camera.f3CameraPos + f3ViewRay * fEndLen;
  
    bool bIsNearHorizonMu = (f4RMuMuSNu.y > fHorizonMu) && (f4RMuMuSNu.y < 0.2 + fHorizonMu);
    int SAMPLE_COUNT = (bIsNearHorizonMu ? 128 : 64) * 3;
    //int SAMPLE_COUNT = 64;
    float fStepLen = fViewRayLenInWorldSpace / SAMPLE_COUNT;

    //float fMaxStepLen = 0.05;
    //if (fStepLen > fMaxStepLen)
    //{
    //    SAMPLE_COUNT = floor(fViewRayLenInWorldSpace / fMaxStepLen);
    //    fStepLen = fViewRayLenInWorldSpace / SAMPLE_COUNT;

    //}

    float3 f3Step = f3ViewRay * fStepLen;

    float3 f3Pos = f3StartPos;
    float fDensityToCam = 0, fTransmittanceToCam = 1, fEnergy = 0;
    float4 f4TotalInscatter = 0;
    float fCoverage, fRainAbsorption, fType;
    float fHeight, fHeightScale;
    float fCloudTest = 0.0;
    int iSampleCountZero = 0;

    for (int i = 0,j = 0; i * 3 + j < SAMPLE_COUNT;)
    {
        //if (fEnergy > 0.99)
        //{
        //    break;
        //}
        if (fDensityToCam>1)
        {
            break;
        }
        GetHeightAndScale(f3Pos, f3EarthCenter, fHeight, fHeightScale);
        SampleWeatherTexture(f3Pos, fCoverage, fRainAbsorption, fType);
        if (fHeightScale < cloud.mCloudTypeLayer[fType].f2LayerHeightScale.x ||
            fHeightScale > cloud.mCloudTypeLayer[fType].f2LayerHeightScale.y)
        {
            f3Pos += f3Step * 3;
            i++;
            continue;
        }
      

        if(fCloudTest>0.0)
        {
            float fSampleDensity = GetCloudDensity(f3Pos, fHeightScale, fCoverage, fType,false);
            if (fSampleDensity == 0.0)
            {
                iSampleCountZero++;
            }
            if (iSampleCountZero != 6)
            {                
                fSampleDensity *= fStepLen;
                fDensityToCam += fSampleDensity;
                fTransmittanceToCam = GetCloudTransmittance(fDensityToCam, fRainAbsorption);
                float fDensityToLight = GetCloudDensityToLight(f3Pos, f3EarthCenter);
                fEnergy += (1 - fEnergy) * (GetCloudTransmittance(saturate(fDensityToLight + fSampleDensity), fRainAbsorption));
                //float fTransimttance = fTransmittanceToCam * GetCloudTransmittance(fDensityToLight, fRainAbsorption);
                //f4TotalInscatter += fSampleDensity * fTransmittanceToCam;
            }
            else
            {
                fCloudTest = 0.0;
                iSampleCountZero = 0;
            }
            f3Pos += f3Step;
            j++;
        }
        else
        {
            fCloudTest = GetCloudDensity(f3Pos, fHeightScale, fCoverage, fType,true);
            if (fCloudTest == 0)
            {
                f3Pos += f3Step * 3;
                i++;
            }
            else
            {
                f3Pos -= f3Step * 3;
                i--;
            }
        }
        
    }
    //f4TotalInscatter.xyz *= HGPhaseFunction(f4RMuMuSNu.w, atmosphere.mie_g);

    //return float4(f4TotalInscatter.xyz, 1);
    return GetCloudTransmittance(fDensityToCam, fRainAbsorption);
    //return fEnergy;

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