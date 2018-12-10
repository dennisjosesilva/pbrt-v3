#include "materials/skin.h"
#include "textures/constant.h"
#include "spectrum.h"
#include "texture.h"
#include "interpolation.h"
#include "paramset.h"
#include "interaction.h"

namespace pbrt
{
	// SkinMaterial Method Definitions
	void SkinMaterial::ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
		bool allowMultipleLobes) const 
	{
		/* Perform bump mapping with _bumpMap_, if present */ 
		if (bumpMap) Bump(bumpMap, si);

		/* Let's represent the oil reflections of the skin using a microfacet brdf with fixed RMS slope in 0.35 */
		si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
		Spectrum ks = oilness->Evaluate(*si).Clamp();

		/* Check about the value of eta - I think it is not mentioned on the article - I am usin 1.0 and 1.0. */
		Float r = roughness->Evaluate(*si);
		/* convert roughness into alpha parameter from the pbrt Beckmann distribution model. */
		Float alpha = BeckmannDistribution::RoughnessToAlpha(r);

		/* I am not sure whether I am using the correct paramenters for Fresnel term (actually I am using "plastic" indices). 
		 * TODO: Check it out. */
		Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.5f, 1.f);
		MicrofacetDistribution *distrib = ARENA_ALLOC(arena, BeckmannDistribution)(alpha, alpha);
		BxDF *spec = ARENA_ALLOC(arena, MicrofacetReflection)(ks, distrib, fresnel);
		si->bsdf->Add(spec);


		/* BSSRDF which represent light interaction with the layers derm and epiderm of the skin. */
		Float Ch = hemoglobinFraction->Evaluate(*si);
		Float Cm = melaninFraction->Evaluate(*si);
		Float Bm = melaninBlend->Evaluate(*si);

		Spectrum sig_a, sig_s;
		ComputeSigmas(Ch, Cm, Bm, &sig_a, &sig_s);
		/* The fourth parameter (1.5f) is a eta parameter I must understand it better. I do not know 
		if it is correct. */ 
		si->bssrdf = ARENA_ALLOC(arena, TabulatedBSSRDF)(*si, this, mode, 1.5f, sig_a, sig_s, table);
	}

	void SkinMaterial::ComputeSigmas(const Float &Ch, const Float &Cm, const Float &Bm, 
		Spectrum *sigma_a, Spectrum *sigma_s) const
	{
		Float gamma = 0.75f;
		Spectrum em, pm, baseline, epiderm;
		int c = 0;
		int step = (sampledLambdaEnd - sampledLambdaStart)/Spectrum::nSamples;

		for (int lambda = sampledLambdaStart; lambda < sampledLambdaEnd; lambda += step) {
			
			/* Compute sigma_a and sigma_s from epiderm */
			em[c] = 6.6e10 * std::pow(lambda, -3.33);
			pm[c] = 2.9e14 * std::pow(lambda, -4.75);
			baseline[c] = 0.0244 + 8.53 * std::exp(-(lambda-154) / 66.2);
			
			(*sigma_a)[c] = Cm*(Bm*em[c] + (1 - Bm)*pm[c]) + (1-Cm)*baseline[c];
			(*sigma_s)[c] = 14.74 * std::pow(lambda, -0.22) + 2.2e11 * std::pow(lambda, -4);
		}
	}

	SkinMaterial *CreateSkinMaterial(const TextureParams &mp)
	{
		std::shared_ptr<Texture<Spectrum>> oilness = mp.GetSpectrumTexture("oilness", Spectrum(0.5f));
		std::shared_ptr<Texture<Float>> hemoglobinFraction = mp.GetFloatTexture("hemoglobinFraction", 0.5f);
		std::shared_ptr<Texture<Float>> melaninFraction = mp.GetFloatTexture("melaninFraction", 0.02f);
		std::shared_ptr<Texture<Float>> melaninBlend = mp.GetFloatTexture("melaninBlend", 0.7f);

		std::shared_ptr<Texture<Float>> roughness = mp.GetFloatTexture("roughness", 0.35f);
		std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");

		return new SkinMaterial(oilness, hemoglobinFraction, melaninFraction, melaninBlend, roughness, bumpMap);
	}

	void SimpleSkinMaterial::ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
		bool allowMultipleLobes) const
	{
		if (bumpMap) Bump(bumpMap, si);

		/* Let's represent the oil reflections of the skin using a microfacet brdf with fixed RMS slope in 0.35 */
		si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, eta); 
		Spectrum ks = oilness->Evaluate(*si).Clamp();

		/* Check about the value eta - I think it is not mentioned on the article - I am using 1.0 and 1.0 */
		Float r = roughness->Evaluate(*si);
		/* Convert roughness into alpha parameter from the pbrt Beckmann distribution model. */
		Float alpha = BeckmannDistribution::RoughnessToAlpha(r);

		/* I am not sure whether I am using the correct parameters for Fresnel term (actually I am using "plastic" indices).
		 * TODO: Check it out. */
		Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric(1.f, eta));
		MicrofacetDistribution *distrib = ARENA_ALLOC(arena, BeckmannDistribution)(alpha, alpha);
		BxDF *spec = ARENA_ALLOC(arena, MicrofacetReflection)(ks, distrib, fresnel);
		si->bsdf->Add(spec);

		si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetTransmission)(
                    Spectrum(1.f), distrib, 1.f, eta, mode));

		Spectrum sig_a = scale * sigma_a->Evaluate(*si).Clamp();
		Spectrum sig_s = scale * sigma_s->Evaluate(*si).Clamp();

		/* GetMediumScatteringProperties("skin2", &sig_a,  &sig_s); */
		si->bssrdf = ARENA_ALLOC(arena, TabulatedBSSRDF)(*si, this, mode, eta, sig_a, sig_s, table);
	}

	SimpleSkinMaterial *CreateSimpleSkinMaterial(const TextureParams &mp)
	{
		Spectrum sigma_a, sigma_s;
		GetMediumScatteringProperties("Skin1", &sigma_a, &sigma_s);

		std::shared_ptr<Texture<Spectrum>> sig_a = mp.GetSpectrumTexture("sigma_a", sigma_a);
		std::shared_ptr<Texture<Spectrum>> sig_s = mp.GetSpectrumTexture("sigma_s", sigma_s);

		std::shared_ptr<Texture<Spectrum>> oilness = mp.GetSpectrumTexture("oilness", Spectrum(1.f));
		std::shared_ptr<Texture<Float>> roughness = mp.GetFloatTexture("roughness", 0.35f);
		Float g = mp.FindFloat("g", 0.0f);
		Float eta = mp.FindFloat("eta", 1.33f);
		Float scale = mp.FindFloat("scale", 1.f);
		std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");

		return new SimpleSkinMaterial(oilness, roughness, sig_a, sig_s, g, eta, scale, bumpMap);
	}

} // namespace pbrt