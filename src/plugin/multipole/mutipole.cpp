#include "multipole.hpp"
#include "external/simple_fft/fft.h"

#include <cmath>
#include <complex>

// static constant
// From the original code - 11 dipole pair is enough to get a good approximation.
static const int NUM_DIPOLE_PAIR = 11;

// --------------------------- Dipole Solver -----------------------------------------------------------------------
DipoleSolver::DipoleSolver(Float eta0, Float eta1, Float sigma_a, Float sigma_s_prime, int zi, Float thickness):
  sigma_a{sigma_a},
  sigma_s_prime{sigma_s_prime},
  d{thickness},
  zi{zi}
{
  sigma_t_prime = sigma_a + sigma_s_prime;
  alpha_prime = sigma_s_prime / sigma_t_prime;
  sigma_tr = sqrtf(3*sigma_a*sigma_t_prime);

  // TODO: check with eta = eta0/eta1 (from the term ratio of index of refraction as described on the paper.  
  Fdr = FresnelDiffuseReflectance(eta0/eta1);
  Float D = 1.f / (3.f * sigma_t_prime);
  Float zb = 2*A(Fdr)*D;
  Float l = 1.f / sigma_t_prime;

  zr = 2.f*zi*(d + 2.f*zb) + l;
  zv = 2.f*zi*(d + 2.f*zb) - l - 2.f*zb;
}

Float DipoleSolver::A(Float Fdr) const
{
  return (1.f + Fdr) / (1.f - Fdr);
}

Float DipoleSolver::FresnelDiffuseReflectance(float eta) const
{
  if (eta >= 1.f) {
    return -(1.4399f/(eta*eta)) + (0.7099f/eta) + 0.6681f + 0.0636f*eta; 
  }
  else {
    Float eta2 = eta*eta;
    return -0.4399f + (0.7099/eta) - (0.3319/eta2) + (0.0636/(eta2*eta));
  }
}

Float DipoleSolver::R(Float r) const
{
  Float D = 1.f/(3.f*sigma_t_prime);

  Float dr = sqrt(r*r + zr*zr);
  Float dv = sqrt(r*r + zv*zv);

  Float denTerm1 = alpha_prime * zr * (1 + sigma_tr * dr) * exp(-sigma_tr*dr);
  Float numTerm1 = 4*pbrt::Pi * dr*dr*dr;
  Float denTerm2 = alpha_prime * zv * (1 + sigma_tr * dv) * exp(-sigma_tr*dv);
  Float numTerm2 = 4*pbrt::Pi * dv*dv*dv;

  return (denTerm1/numTerm1) - (denTerm2/numTerm2);
}

Float DipoleSolver::T(Float r) const
{
  Float dr = sqrt(r*r + zr*zr);
  Float dv = sqrt(r*r + zv*zv);

  Float denTerm1 = alpha_prime*(d - zr) * (1.f + sigma_tr*dr) * exp(-sigma_tr*dr);
  Float numTerm1 = 4*pbrt::Pi * dr*dr*dr;
  Float denTerm2 = alpha_prime*(d - zv) * (1.f + sigma_tr*dv) * exp(-sigma_tr*dv);
  Float numTerm2 = 4*pbrt::Pi * dv*dv*dv;

  return (denTerm1 / numTerm1) - (denTerm2 / numTerm2);
}


// -------------------------- Multipole Table -----------------------------------------------------------
MultipoleTable::MultipoleTable(std::size_t nSample):
  m_transmitance(nSample), 
  m_reflectance(nSample), 
  m_squaredDistance(nSample)
{}

MultipoleTable::MultipoleTable()
{}

void MultipoleTable::PushBack(Float reflectance, Float transmitance, Float squaredDistance)
{
  m_reflectance.push_back(reflectance);
  m_transmitance.push_back(transmitance);
  m_squaredDistance.push_back(squaredDistance);
}


// ---------------------- Compute Diffusion Profile ------------------------------------------------------
// Extracted from the original code 
// https://github.com/patwonder/pbrt-v2-skin/blob/master/src/multipole/MultipoleProfileCalculator/MultipoleProfileCalculator.cpp
unsigned int RoundUpPow2(unsigned int v) 
{
  v--;
  v |= v >> 1; v |= v >> 2;
  v |= v >> 4; v |= v >> 8;
  v |= v >> 16; 
  return v+1;
}

MatrixProfile::MatrixProfile(unsigned int length):
  reflectance(length, length), transmitance(length, length)
{}

MatrixProfile ComputeLayerProfile(const MultipoleLayer &layer, Float stepSize, unsigned int profile_length)
{
  MatrixProfile profile{profile_length};
  unsigned int center = (profile_length - 1) / 2;
  unsigned int extent = center;

  // fill dipole solvers
  std::vector<DipoleSolver> dss;
  dss.reserve(NUM_DIPOLE_PAIR);
  int half_num_dipoles = (NUM_DIPOLE_PAIR-1) / 2;
  for (auto i = -half_num_dipoles; i < half_num_dipoles; ++i) {
    dss.push_back(DipoleSolver{layer.eta0, layer.eta1, layer.sigma_a, layer.sigma_s_prime, i, layer.thickness});
  }

  Float normalizeFactor = stepSize * stepSize;
  for (unsigned int row = 0; row <= extent; row++) {
    for (unsigned int col = row; col < extent; col++) {
      real_type row2 = row * row;
      real_type col2 = col * col;
      real_type r2 = (row2 + col2) * (normalizeFactor);
      for (const DipoleSolver &ds : dss) {
        real_type r = ds.R((Float)r2) * normalizeFactor;
        real_type t = ds.T((Float)r2) * normalizeFactor;
        profile.reflectance(center + row, center + col) += r;
        profile.transmitance(center + row, center + col)  += t;
      }
    }
  }

  // Fill the rest of the matrix  through simmetry - also very similar with the original code.
  for (unsigned int row = 1; row <= extent; row++) {
    for (unsigned int col = 0; col < row; col++) {
      profile.reflectance(center + row, center + col) = profile.reflectance(center + col, center + row);
      profile.transmitance(center + row, center + col) = profile.transmitance(center + col, center + row);
    }
  }

  for (unsigned int row = 0; row <= extent; row++) {
    for (unsigned int col = 0; col <= extent; col++) {
      real_type r = profile.reflectance(center + row, center + col);
      profile.reflectance(center - row, center + col) = r;
      profile.reflectance(center + row, center + col) = r;
      profile.reflectance(center - row, center - col) = r;

      real_type t = profile.transmitance(center + row, center + col);
      profile.transmitance(center - row, center + col) = t;
      profile.transmitance(center + row, center - col) = t;
      profile.transmitance(center - row, center - col) = t;
    }
  }

  return profile;
}

// The below two functions are also very similar with the original source:
FFTMatrix<complex_type> runFFT(FFTMatrix<Float> &matrix) 
{
  unsigned int length = matrix.NRows();
  unsigned int convolution_length = convolution_length / 2;
  unsigned int center = (length - 1) / 2;
  std::cout << convolution_length << " -- " << length << " ---- " << center << "\n";
  FFTMatrix<complex_type> out(convolution_length, convolution_length / 2 + 1);
  const char *error_message = nullptr;
  bool success = simple_fft::FFT(matrix.ScaleAndShift(convolution_length, convolution_length, center, center),
    out, convolution_length, convolution_length, error_message);
  return out; 
}

FFTMatrix<Float> runIFFT(FFTMatrix<complex_type> &matrix) 
{
  unsigned int convolution_length = matrix.NRows();
  unsigned int length = convolution_length / 2;
  unsigned int center = (length - 1) / 2;
  FFTMatrix<complex_type> temp(convolution_length, convolution_length);
  FFTMatrix<Float> out(convolution_length, convolution_length);
  const char *error_message = nullptr;
  bool success = simple_fft::IFFT(matrix, temp, convolution_length, convolution_length, error_message);
  
  for (unsigned int i = 0; i < temp.NElements(); ++i) {
    out[i] = temp[i].real();
  }

  return out.ScaleAndShiftReversed(length, length, center, center);
}

MatrixProfile CombineProfiles(MatrixProfile &layer1, MatrixProfile &layer2)
{
  MatrixProfile combined{layer1.length()};

  FFTMatrix<complex_type> fR1 = runFFT(layer1.reflectance);
  FFTMatrix<complex_type> fR2 = runFFT(layer2.reflectance);
  FFTMatrix<complex_type> fT1 = runFFT(layer1.transmitance);
  FFTMatrix<complex_type> fT2 = runFFT(layer2.transmitance);

  FFTMatrix<complex_type> f1MinusR1TimesR2 = fR2*fR1;
  f1MinusR1TimesR2.OneMinusSelf();
  
  FFTMatrix<complex_type> fR12 = fR1 + (fT1*fR2*fT1/f1MinusR1TimesR2);
  FFTMatrix<complex_type> fT12 = (fT1 * fT2)/f1MinusR1TimesR2;

  /*combined.reflectance = runIFFT(fR12);
  combined.transmitance = runIFFT(fT12);*/

  return combined;
} 


MultipoleTable ComputeMultipoleDiffusionProfile(const std::vector<MultipoleLayer> &layers, 
  const MultipoleOptions &options)
{
  unsigned int length = RoundUpPow2(options.desiredLength);
  MatrixProfile mp0{length * 2};

  mp0 = ComputeLayerProfile(layers[0], options.desiredLength, mp0.length());
  for (unsigned int i = 1; i < layers.size(); ++i) {
    MatrixProfile mp1 = ComputeLayerProfile(layers[i], options.desiredLength, mp1.length());
    mp0 = CombineProfiles(mp0, mp1);
  }

  MultipoleTable table;
  unsigned int center = length - 1;
  unsigned int extent = center;
  float denormalizeFactor = 1.f / (options.desiredLength * options.desiredLength);
  for (unsigned int  i = 0; i <= extent; ++i) {
    for (unsigned int j = i; (i*i + j*j) <= extent; j++) {
      Float r = mp0.reflectance(center + i, center + j)  * denormalizeFactor;
      Float t = mp0.transmitance(center + i, center + j) * denormalizeFactor;
      Float sd = (float)(i*i + j*j) * options.desiredLength * options.desiredLength;
      table.PushBack(r, t, sd);
    }
  }

  return table;
}