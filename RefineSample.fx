#include "Common.fx"

RWTexture2D<uint2> g_rwtex2DInterpolationSource : register(u0);

#ifndef SAMPLE_STEP
#define SAMPLE_STEP 128
#endif

#ifndef THREAD_GROUP_SIZE
#define THREAD_GROUP_SIZE max(SAMPLE_STEP, 32)
#endif

static const uint g_uiPackNum = THREAD_GROUP_SIZE / 32;
groupshared uint g_uiCamDepthDiffPackFlags[g_uiPackNum];

static const float fRefinementThreshold = 0.03;
static const float fSampleDense = 2;

[numthreads(THREAD_GROUP_SIZE,1,1)]
void RefineSampleCS(uint3 Gid : SV_GroupID,
                    uint3 GTid : SV_GroupThreadID)
{
    uint uiSliceNum = Gid.y;
    uint uiSampleGroupStart = Gid.x * THREAD_GROUP_SIZE;
    uint uiSampleGlobalNum = uiSampleGroupStart + GTid.x;
    
    float2 f2SampleScreenXY = g_tex2DEpipolarSample.Load(uint3(uiSampleGlobalNum, uiSliceNum, 0));

    bool IsValidThread = all(abs(f2SampleScreenXY) < 1 + 1e-4); 

    if(Gid.x < g_uiPackNum)
        g_uiCamDepthDiffPackFlags[Gid.x] = 0;
    GroupMemoryBarrierWithGroupSync();

    [branch]
    if(IsValidThread)
    {
        float fSampleCamDepth = g_tex2DEpipolarSampleCamDepth.Load(uint3(uiSampleGlobalNum, uiSliceNum, 0));
        float fSampleCamDepthRight = g_tex2DEpipolarSampleCamDepth.Load(uint3(min(uiSampleGlobalNum + 1, EPIPOLAR_SAMPLE_NUM - 1), uiSliceNum, 0));
        float fMax = max(fSampleCamDepth, fSampleCamDepthRight);
        fMax = max(fMax, 1 - 1e-4) + 1e-4;
        bool bFlag = abs(fSampleCamDepth - fSampleCamDepthRight) / fMax < 0.2 * fRefinementThreshold; // 1 No break
                                                                                                      // 0 Depth break
        InterlockedOr(g_uiCamDepthDiffPackFlags[GTid.x / 32], bFlag << (GTid.x % 32));
    }
    GroupMemoryBarrierWithGroupSync();

    if(!IsValidThread)
        return;

    uint uiSampleStep = SAMPLE_STEP;
    uint uiSampleLocalNum0 = (GTid.x / uiSampleStep) * uiSampleStep;
    uint uiSampleGlobalNum0 = uiSampleLocalNum0 + uiSampleGroupStart;
    float2 f2SampleScreenXY0 = g_tex2DEpipolarSample.Load(uint3(uiSampleGlobalNum0, uiSliceNum, 0));
    if (length(f2SampleScreenXY0 - light.f2LightScreenPos) < 0.1 &&
        (float)uiSampleGlobalNum0 / EPIPOLAR_SAMPLE_NUM < 0.05)
    {
        uiSampleStep = max(uiSampleStep / fSampleDense, 1);
        uiSampleLocalNum0 = (GTid.x / uiSampleStep) * uiSampleStep;
    }
    uint uiSampleLocalNum1 = uiSampleLocalNum0 + uiSampleStep;
    if (Gid.x == EPIPOLAR_SAMPLE_NUM / THREAD_GROUP_SIZE - 1)
        uiSampleLocalNum1 = min(uiSampleLocalNum1, THREAD_GROUP_SIZE - 1);
    uiSampleStep = uiSampleLocalNum1 - uiSampleLocalNum0;

    uint uiSampleLeftBoundary = uiSampleLocalNum0;
    uint uiSampleRightBoundary = uiSampleLocalNum1;
    if (GTid.x > uiSampleLocalNum0 && GTid.x < uiSampleLocalNum1)
    {
        uint uiCamDepthDiffPackFlags[g_uiPackNum];
        for (int i = 0; i < g_uiPackNum; i++)
        {
            uiCamDepthDiffPackFlags[i] = g_uiCamDepthDiffPackFlags[i];
        }
        
        bool bNoBreak = true;

        uint uiPackNum0 = uiSampleLocalNum0 / 32;
        uint uiNumInPack0 = uiSampleLocalNum0 % 32;

        uint uiPackNum1 = uiSampleLocalNum1 / 32;
        uint uiNumInPack1 = uiSampleLocalNum1 % 32;

        for (int i = uiPackNum0; i <= uiPackNum1 && bNoBreak ; i++)
        {
            if (uiCamDepthDiffPackFlags[i] != 0xFFFFFFFFU)
                bNoBreak = false;
        }

        if(bNoBreak)
        {
            uint uiSampleLeft = GTid.x - 1; // GTid.x != 0
            uint uiSampleLeftPackNum = uiSampleLeft / 32;

            int iFirstBreakFlag = -1;
            while (iFirstBreakFlag == -1 && uiSampleLeftPackNum >= uiPackNum0)
            {
                uint uiFlag = uiCamDepthDiffPackFlags[uiSampleLeftPackNum];
                uint uiNumInPackSampleLeft = uiSampleLeft % 32;

                if (uiNumInPackSampleLeft < 31)
                {
                    uiFlag |= (uint(0xFFFFFFFFU) << (uiNumInPackSampleLeft + 1));
                }
                uint iFirstBreakFlagNum = firstbithigh(~uiFlag);
                if (!(iFirstBreakFlagNum >= 0 && iFirstBreakFlagNum <= 31))
                    iFirstBreakFlag = -1;
                uiSampleLeft -= uiNumInPackSampleLeft - iFirstBreakFlag;
                uiSampleLeftPackNum--;
            }


        }
        //if (uiPackNum0 == uiPackNum1)
        //{
        //    uint uiFlagMask = (1 << uiSampleStep) - 1;
        //    uint uiFlag = (uiCamDepthDiffPackFlags[uiPackNum0] >> uiPackNum0) & uiFlagMask;
        //    if (uiFlag!=uiFlagMask)
        //        bNoBreak = false;
        //}
        //else
        //{
        //    //uint uiFlagMask = (1 << (32 - uiNumInPack0)) - 1;
        //    //uint uiFlag = (uiCamDepthDiffPackFlags[uiPackNum0] >> uiPackNum0) & uiFlagMask;
        //    //if (uiFlag != uiFlagMask)
        //    //    bNoBreak = false;
        //    //for (int i = uiPackNum0 + 1; i < uiPackNum1 && bNoBreak; i++)
        //    //{
        //    //    if (uiCamDepthDiffPackFlags[i] != 0xFFFFFFFFU)
        //    //        bNoBreak = false;
        //    //}
        //}
        

    }
    else
    {
        uiSampleLeftBoundary = uiSampleRightBoundary = GTid.x;
    }
    g_rwtex2DInterpolationSource[uint2(uiSampleGlobalNum, uiSliceNum)] = uint2(uiSampleGroupStart + uiSampleLeftBoundary, uiSampleGroupStart + uiSampleRightBoundary);
}

technique11 RefineSampleTech
{
    pass
    {
        SetVertexShader(NULL);
        SetGeometryShader(NULL);
        SetPixelShader(NULL);
        SetComputeShader(CompileShader(cs_5_0, RefineSampleCS()));
    }
}