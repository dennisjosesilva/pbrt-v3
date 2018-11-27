#pragma once

#ifndef MULTIPOLE_HPP_INCLUDED
#define MULTIPOLE_HPP_INCLUDED

// --------------------------- Dipole Solver ---------------------------------------------------

class DipoleSolver
{
public:
  DipoleSolver(Float eta0, Float eta1, Float sigma_a, Float sigma_s_prime, Float thickness);

  Float R(Float r);
  Float T(Float r);

private:
	Float A(Float rho_d);
	Float FresnelDiffuseReflectance(float eta);

private:
  Float sigma_a;
  Float sigma_s_prime;
  Float sigma_t_prime;
  Float alpha_prime;
  Float sigma_tr;
  Float Fdr;
  Float zr;
  Float zv;
  Float dr;
  Float dv;
  Float d;
};

// ------------------------ Multipole Layer ------------------------------------------------------
struct MultipoleLayer
{
	Float eta0;
	Float eta1;
	Float thickness;
	Float sigma_a;
	Float sigma_s_prime;
};


// ----------------- Multipole options ----------------------------------------------------------
struct MultipoleOptions
{
	Float desiredStepSize;
	std::size_t desiredLength;
};

// -------------------- Multipole Table ---------------------------------------------------------
class MultipoleTable
{
public:
	MultipoleTable(std::size_t nSamples);

	inline Float &transmitance(int index) { return m_transmitance[index]; }
	inline Float transmitance(int index) const { return m_transmitance[index]; }
	Float &transmitance(Float squaredDistance);
	Float transmitance(Float squaredDistance) const;

	inline Float &reflectance(int index) { return m_reflectance[index]; }
	inline Float reflectance(int index) const { return m_reflectance[index]; }
	Float &reflectance(Float squaredDistance);
	Float reflectance(Float squaredDistance) const;

	inline Float &squaredDistance(int index) { return m_squaredDistance[index]; }
	inline Float squaredDistance(int index) const { return m_squaredDistance[index]; }

private:
	std::vector<Float> m_transmitance;
	std::vector<Float> m_reflectance;
	std::vector<Float> m_squaredDistance;
}


#endif