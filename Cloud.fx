#include "Common.fx"

//static const int perm[256] =
//{
//    151, 160, 137, 91, 90, 15,
//    131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
//    190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
//    88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
//    77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
//    102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
//    135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
//    5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
//    223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
//    129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
//    251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
//    49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
//    138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
//};                                                                                        

//float grad(int hash, float x, float y, float z)
//{
//    int h = hash & 15; // Convert low 4 bits of hash code into 12 simple
//    float u = h < 8 ? x : y; // gradient directions, and compute dot product.
//    float v = h < 4 ? y : h == 12 || h == 14 ? x : z; // Fix repeats at h = 12 to 15
//    return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
//}

//float Perlin(float x, float y, float z)
//{
//    float n0, n1, n2, n3; // Noise contributions from the four corners

//    // Skewing/Unskewing factors for 3D
//    static const float F3 = 1.0f / 3.0f;
//    static const float G3 = 1.0f / 6.0f;

//    // Skew the input space to determine which simplex cell we're in
//    float s = (x + y + z) * F3; // Very nice and simple skew factor for 3D
//    int i = floor(x + s);
//    int j = floor(y + s);
//    int k = floor(z + s);
//    float t = (i + j + k) * G3;
//    float X0 = i - t; // Unskew the cell origin back to (x,y,z) space
//    float Y0 = j - t;
//    float Z0 = k - t;
//    float x0 = x - X0; // The x,y,z distances from the cell origin
//    float y0 = y - Y0;
//    float z0 = z - Z0;

//    // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
//    // Determine which simplex we are in.
//    int i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
//    int i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
//    if (x0 >= y0)
//    {
//        if (y0 >= z0)
//        {
//            i1 = 1;
//            j1 = 0;
//            k1 = 0;
//            i2 = 1;
//            j2 = 1;
//            k2 = 0; // X Y Z order
//        }
//        else if (x0 >= z0)
//        {
//            i1 = 1;
//            j1 = 0;
//            k1 = 0;
//            i2 = 1;
//            j2 = 0;
//            k2 = 1; // X Z Y order
//        }
//        else
//        {
//            i1 = 0;
//            j1 = 0;
//            k1 = 1;
//            i2 = 1;
//            j2 = 0;
//            k2 = 1; // Z X Y order
//        }
//    }
//    else
//    { // x0<y0
//        if (y0 < z0)
//        {
//            i1 = 0;
//            j1 = 0;
//            k1 = 1;
//            i2 = 0;
//            j2 = 1;
//            k2 = 1; // Z Y X order
//        }
//        else if (x0 < z0)
//        {
//            i1 = 0;
//            j1 = 1;
//            k1 = 0;
//            i2 = 0;
//            j2 = 1;
//            k2 = 1; // Y Z X order
//        }
//        else
//        {
//            i1 = 0;
//            j1 = 1;
//            k1 = 0;
//            i2 = 1;
//            j2 = 1;
//            k2 = 0; // Y X Z order
//        }
//    }

//    // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
//    // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
//    // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
//    // c = 1/6.
//    float x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
//    float y1 = y0 - j1 + G3;
//    float z1 = z0 - k1 + G3;
//    float x2 = x0 - i2 + 2.0f * G3; // Offsets for third corner in (x,y,z) coords
//    float y2 = y0 - j2 + 2.0f * G3;
//    float z2 = z0 - k2 + 2.0f * G3;
//    float x3 = x0 - 1.0f + 3.0f * G3; // Offsets for last corner in (x,y,z) coords
//    float y3 = y0 - 1.0f + 3.0f * G3;
//    float z3 = z0 - 1.0f + 3.0f * G3;

//    // Work out the hashed gradient indices of the four simplex corners
//    int gi0 = perm[i + perm[j + perm[k]]];
//    int gi1 = perm[i + i1 + perm[j + j1 + perm[k + k1]]];
//    int gi2 = perm[i + i2 + perm[j + j2 + perm[k + k2]]];
//    int gi3 = perm[i + 1 + perm[j + 1 + perm[k + 1]]];

//    // Calculate the contribution from the four corners
//    float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
//    if (t0 < 0)
//    {
//        n0 = 0.0;
//    }
//    else
//    {
//        t0 *= t0;
//        n0 = t0 * t0 * grad(gi0, x0, y0, z0);
//    }
//    float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
//    if (t1 < 0)
//    {
//        n1 = 0.0;
//    }
//    else
//    {
//        t1 *= t1;
//        n1 = t1 * t1 * grad(gi1, x1, y1, z1);
//    }
//    float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
//    if (t2 < 0)
//    {
//        n2 = 0.0;
//    }
//    else
//    {
//        t2 *= t2;
//        n2 = t2 * t2 * grad(gi2, x2, y2, z2);
//    }
//    float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
//    if (t3 < 0)
//    {
//        n3 = 0.0;
//    }
//    else
//    {
//        t3 *= t3;
//        n3 = t3 * t3 * grad(gi3, x3, y3, z3);
//    }
//    // Add contributions from each corner to get the final noise value.
//    // The result is scaled to stay just inside [-1,1]
//    return 32.0f * (n0 + n1 + n2 + n3);
//}

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
        result += amplitude * Perlin(uvw);
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
    
    result.x = (saturate(PerlinfBm(f3UVW * 32, 4, 0.5) * 0.5 + 0.5) + WorleyfBm(f3UVW, 4, 0.5, 8, 1)) / 2;
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