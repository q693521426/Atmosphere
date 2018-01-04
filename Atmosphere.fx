//  'exp_term' * exp('exp_scale' * h) + 'linear_term' * h + 'constant_term'
// clamp [0,1]
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

cbuffer cbAtmosphere
{

};

Texture2D transmittance_texture;

float GetLayerDesity(DensityProfileLayer layer, float altitude)
{
	float desity = layer.exp_term * exp(layer.scale_term * altitude) + layer.linear_term * altitude + layer.const_term;
	return clamp(desity, 0.f, 1.f);
}

float GetProfileDesity(DensityProfile profile, float altitude)
{
	return altitude < profile.layer[0].width ?
		GetLayerDesity(profile.layer[0], altitude) : GetLayerDesity(profile.layer[1], altitude);
}

float DistanceToTopAtmosphereBoundary(AtmosphereParameters atmosphere, float r, float mu)
{
	float discriminant = r * r * (mu * mu - 1) + atmosphere.top_radius * atmosphere.top_radius;
	return max(-r * mu + sqrt(max(discriminant, 0.f)), 0.f);
}

float DistanceToBottomAtmosphereBoundary(AtmosphereParameters atmosphere, float r, float mu)
{
	float discriminant = r * r * (mu * mu - 1) + atmosphere.bottom_radius * atmosphere.bottom_radius;
	return max(-r * mu - sqrt(max(discriminant, 0.f)), 0.f);
}

float ComputeOpticalLengthToTopAtmosphereBoundary(AtmosphereParameters atmosphere, DensityProfile profile,
												float r, float mu)
{
	const int SAMPLE_COUNT = 500;
	float dx = DistanceToTopAtmosphereBoundary(atmosphere, r, mu);

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

float ComputeTransmittanceToTopAtmosphereBoundary(AtmosphereParameters atmosphere, float r, float mu)
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
				ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere, atmosphere.rayleigh_density, r, mu) +
		atmosphere.mie_scattering * 
				ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere, atmosphere.mie_density, r, mu)));
}

float GetTransmittance(AtmosphereParameters atmosphere, float r, float mu)
{

}