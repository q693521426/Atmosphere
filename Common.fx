
static const float Pi = 3.141592654;

static const int nSamples = 2;
static const float fSamples = (float)nSamples;

float scale(float fCos,float fScaleDepth)
{
	float x = 1.0 - fCos;
	return fScaleDepth * exp(-0.00287 + x*(0.459 + x*(3.83 + x*(-6.80 + x*5.25))));
}

float getMiePhase(float fCos, float fCos2, float g, float g2)
{
	//return 1.5 * ((1.0 - g2) / (2.0 + g2)) * (1.0 + fCos2) / pow(1.0 + g2 - 2.0*g*fCos, 1.5);

	float k = 1.55f *g - 0.55f * g2;
	return 1.f/(4.f * Pi)*(1.f - k*k)/((1.f-k*fCos)*(1.f-k*fCos));
}

float getRayleighPhase(float fCos2)
{
	return 0.75 + 0.75*fCos2;
}

float2 getIntersection(float3 v3Pos,float3 v3Ray,float fDistance2, float fRadius2)
{
	float b = dot(v3Pos,v3Ray);
	float c = fDistance2 - fRadius2;
	float fDet = sqrt(max(0.0,b*b - c));
	return float2(-b-fDet,-b+fDet);
}