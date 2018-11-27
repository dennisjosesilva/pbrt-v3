#include "multipole.hpp"
#include "../../core/pbrt.h"
#include <cmath>

DipoleSolver::DipoleSolver(Float sigma_a, Float sigma_s_prime, Float d):
  sigma_a{a},
  sigma_s_prime{sigma_s_prime},
  d{d}
{}

Float DipoleSolver::R(Float r)
{
  Float sigma_t_prime = sigma_a + sigma_s_prime;
  Float alpha_prime = sigma_s_prime / sigma_t_prime;
  Float sigma_tr = sqrtf(3*sigma_a*sigma_t_prime);
  Float zr = (1.f / sigma_t_prime);
  Float zv = (1.f + 4.f*A/3.f)/sigma_t_prime;

  Float dr = sqrt(r*r + zr*zr);
  Float dv = sqrt(r*r + zv*zv);
  Float D = 1.f/(3.f*sigma_t_prime);

  Float denTerm1 = alpha_prime * zr * (1 + sigma_tr * dr) * exp(-sigma_tr*dr);
  Float numTerm1 = 4*pbrt::Pi * dr*dr*dr;
  Float denTerm2 = alpha_prime * zv * (1 + sigma_tr * dv) * exp(-sigma_tr*dv);
  Float numTerm2 = 4*pbrt::Pi * dv*dv*dv;

  return (denTerm1/numTerm1) - (denTerm2/numTerm2);
}

Float DipoleSolver::T(Float r)
{
  Float alpha_prime = sigma_s_prime / sigma_t_prime;
  Float sigma_tr = sqrtf(3*sigma_a*sigma_t_prime);
  Float zr = (1.f / sigma_t_prime);
  Float zv = (1.f + 4.f*A/3.f)/sigma_t_prime;
  Float dr = sqrt(r*r + zr*zr);
  Float dv = sqrt(r*r + zv*zv);

  Float denTerm1 = alpha_prime*(d - zr) * (1.f + sigma_tr*dr) * exp(-sigma_tr*dr);
  Float numTerm1 = 4*pbrt::Pi * dr*dr*dr;
  Float denTerm2 = alpha_prime*(d - zv) * (1.f + sigma_tr*dv) * exp(-sigma_tr*dv);
  Float numTerm2 = 4*pbrt::Pi * dv*dv*dv;

  return (denTerm1 / numTerm1) - (denTerm2 / numTerm2);
}
