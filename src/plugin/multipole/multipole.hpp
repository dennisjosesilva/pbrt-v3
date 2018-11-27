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


// -------------------- Multipole Table ---------------------------------------------------------


#endif