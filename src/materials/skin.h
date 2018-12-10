/* The material will represet the human skin. It will be developed based on 
 * "A Spectral BSSRDF for Shading Human Skin" paper by Craig Donner and 
 * Henrik Wann Jensen.
*/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_SKIN_H
#define PBRT_MATERIALS_SKIN_H

#include "pbrt.h"
#include "material.h"
#include "reflection.h"
#include "bssrdf.h"

namespace pbrt 
{
	/* SkinMaterial Declarations */
	class SkinMaterial : public Material 
	{
	public:
		/* SkinMaterial Public Methods */
		SkinMaterial(const std::shared_ptr<Texture<Spectrum>> &oilness,
					 const std::shared_ptr<Texture<Float>> &hemoglobinFraction,
					 const std::shared_ptr<Texture<Float>> &melaninFraction,
					 const std::shared_ptr<Texture<Float>> &melaninBlend,
					 const std::shared_ptr<Texture<Float>> &roughness,
					 const std::shared_ptr<Texture<Float>> &bumpMap):
		oilness(oilness),
		hemoglobinFraction(hemoglobinFraction),
		melaninFraction(melaninFraction),
		melaninBlend(melaninBlend),
		roughness(roughness),
		bumpMap(bumpMap),
		table(100, 64) {
			/* g = 0.0, eta = 1.33*/
			ComputeBeamDiffusionBSSRDF(0, 1.33, &table);
		}

		void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, 
			bool allowMultipleLobes) const;

	private:
		void ComputeSigmas(const Float &Ch, const Float &Cm, const Float &Bm, 
			Spectrum *sigma_a, Spectrum *sigma_s) const;

	private:
		std::shared_ptr<Texture<Spectrum>> oilness;
		std::shared_ptr<Texture<Float>> hemoglobinFraction;
		std::shared_ptr<Texture<Float>> melaninFraction;
		std::shared_ptr<Texture<Float>> melaninBlend;

		std::shared_ptr<Texture<Float>> roughness;		
		std::shared_ptr<Texture<Float>> bumpMap;
		BSSRDFTable table;
	};

	SkinMaterial *CreateSkinMaterial(const TextureParams &mp);



	/* Simple Skin Material */
	class SimpleSkinMaterial : public Material
	{
	public:
		SimpleSkinMaterial(const std::shared_ptr<Texture<Spectrum>> oilness,
						   const std::shared_ptr<Texture<Float>> roughness,
						   const std::shared_ptr<Texture<Spectrum>> sigma_a,
						   const std::shared_ptr<Texture<Spectrum>> sigma_s,
						   Float g,
						   Float eta,
						   Float scale,
						   const std::shared_ptr<Texture<Float>> bumpMap):
		 oilness(oilness),
		 roughness(roughness),
		 sigma_a(sigma_a),
		 sigma_s(sigma_s),
		 g(g),
		 eta(eta),
		 scale(scale),
		 bumpMap(bumpMap),
		 table(100, 64)
		 {
		 	ComputeBeamDiffusionBSSRDF(g, eta, &table);
		 }

		void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
			bool allowMultipleLobes) const;
	private:
		std::shared_ptr<Texture<Spectrum>> oilness;
		std::shared_ptr<Texture<Float>> roughness;
		std::shared_ptr<Texture<Spectrum>> sigma_a;
		std::shared_ptr<Texture<Spectrum>> sigma_s;
		Float g;
		Float eta;
		Float scale;
		std::shared_ptr<Texture<Float>> bumpMap;
		BSSRDFTable table;
	};

	SimpleSkinMaterial *CreateSimpleSkinMaterial(const TextureParams &mp);
}

#endif