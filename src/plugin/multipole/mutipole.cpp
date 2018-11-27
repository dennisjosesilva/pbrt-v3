#include "multipole.hpp"
#include "../../core/pbrt.h"
#include <cmath>


// --------------------------- Dipole Solver -----------------------------------------------------------------------
DipoleSolver::DipoleSolver(Float eta0, Float eta1, Float sigma_a, Float sigma_s_prime, Float thickness):
  sigma_a{a},
  sigma_s_prime{sigma_s_prime},
  d{thickness}
{
  sigma_t_prime = sigma_a + sigma_s_prime;
  alpha_prime = sigma_s_prime / sigma_t_prime;
  sigma_tr = sqrtf(3*sigma_a*sigma_t_prime);
  
  // TODO: check with eta = eta0/eta1 (from the term ratio of index of refraction as described on the paper.  
  Fdr = FresnelDiffuseReflectance(eta0/eta1);

  zr = (1.f / sigma_t_prime);
  zv = (1.f + 4.f*A(Fdr)/3.f)/sigma_t_prime;

  dr = sqrt(r*r + zr*zr);
  dv = sqrt(r*r + zv*zv);
}

Float DipoleSolver::A(Float Fdr)
{
  return (1.f + Fdr) / (1.f - Fdr);
}

Float FresnelDiffuseReflectance(float eta)
{
  if (eta >= 1.f) {
    return -(1.4399f/(eta*eta)) + (0.7099f/eta) + 0.6681f + 0.0636f*eta; 
  }
  else {
    Float eta2 = eta*eta;
    return -0.4399f + (0.7099/eta) - (0.3319/eta2) + (0.0636/(eta2*eta));
  }
}

Float DipoleSolver::R(Float r)
{
  Float D = 1.f/(3.f*sigma_t_prime);

  Float denTerm1 = alpha_prime * zr * (1 + sigma_tr * dr) * exp(-sigma_tr*dr);
  Float numTerm1 = 4*pbrt::Pi * dr*dr*dr;
  Float denTerm2 = alpha_prime * zv * (1 + sigma_tr * dv) * exp(-sigma_tr*dv);
  Float numTerm2 = 4*pbrt::Pi * dv*dv*dv;

  return (denTerm1/numTerm1) - (denTerm2/numTerm2);
}

Float DipoleSolver::T(Float r)
{
  Float denTerm1 = alpha_prime*(d - zr) * (1.f + sigma_tr*dr) * exp(-sigma_tr*dr);
  Float numTerm1 = 4*pbrt::Pi * dr*dr*dr;
  Float denTerm2 = alpha_prime*(d - zv) * (1.f + sigma_tr*dv) * exp(-sigma_tr*dv);
  Float numTerm2 = 4*pbrt::Pi * dv*dv*dv;

  return (denTerm1 / numTerm1) - (denTerm2 / numTerm2);
}


// -------------------------- Multipole Solver -----------------------------------------------------------
