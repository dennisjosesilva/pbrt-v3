#pragma once

#ifndef MULTIPOLE_HPP_INCLUDED
#define MULTIPOLE_HPP_INCLUDED

#include "pbrt.h"
#include <vector>
#include "external/simple_fft/fft_settings.h"
#include <algorithm>

// --------------------------- Dipole Solver ---------------------------------------------------

using Float = pbrt::Float;
typedef Float real_type;


class DipoleSolver
{
public:
  DipoleSolver(Float eta0, Float eta1, Float sigma_a, Float sigma_s_prime, int zi, Float thickness);

  Float R(Float r) const;
  Float T(Float r) const;

private:
	Float A(Float rho_d) const;
	Float FresnelDiffuseReflectance(float eta) const;

private:
  Float sigma_a;
  Float sigma_s_prime;
  Float sigma_t_prime;
  Float alpha_prime;
  Float sigma_tr;
  Float Fdr;
  Float zr;
  Float zv;
  int zi;
  Float d;
};

// ------------------------ Multipole Layer ------------------------------------------------------
struct MultipoleLayer
{
	MultipoleLayer(Float eta0, Float eta1, Float thickness, Float sigma_a, Float sigma_s_prime)
		:eta0{eta0}, eta1{eta1}, thickness{thickness}, sigma_a{sigma_a}, sigma_s_prime{sigma_s_prime}
	{}

	Float eta0;
	Float eta1;
	Float thickness;
	Float sigma_a;
	Float sigma_s_prime;
};


// ----------------- Multipole options ----------------------------------------------------------
struct MultipoleOptions
{
	MultipoleOptions(Float desiredStepSize, std::size_t desiredLength)
		:desiredStepSize{desiredStepSize}, desiredLength{desiredLength}
	{}

	Float desiredStepSize;
	std::size_t desiredLength;
};

// -------------------- Multipole Table ---------------------------------------------------------
class MultipoleTable
{
public:
	MultipoleTable(std::size_t nSamples);
	MultipoleTable();

	inline Float &transmitance(int index) { return m_transmitance[index]; }
	inline Float transmitance(int index) const { return m_transmitance[index]; }
	/*Float &transmitance(Float squaredDistance);
	Float transmitance(Float squaredDistance) const;*/

	inline Float &reflectance(int index) { return m_reflectance[index]; }
	inline Float reflectance(int index) const { return m_reflectance[index]; }
	/*Float &reflectance(Float squaredDistance);
	Float reflectance(Float squaredDistance) const;*/

	inline Float &squaredDistance(int index) { return m_squaredDistance[index]; }
	inline Float squaredDistance(int index) const { return m_squaredDistance[index]; }

	void PushBack(Float reflectance, Float transmitance, Float squaredDistance);
	std::size_t NSamples() const { return m_reflectance.size(); }

private:
	std::vector<Float> m_transmitance;
	std::vector<Float> m_reflectance;
	std::vector<Float> m_squaredDistance;
};


// ----------------------- Fast Fourier Matrix "adapter" ---------------------------------------------
template<typename T>
class FFTMatrix
{
public:
  FFTMatrix(unsigned int nrows, unsigned int ncols)
  	:m_nrows(nrows), m_ncols(ncols), m_data(m_nrows*m_ncols)
  {}

  inline T& operator[](unsigned int index) { return m_data[index]; }
  inline T operator[](unsigned int index) const { return m_data[index]; }

  inline T& operator() (unsigned int row, unsigned int col) { return m_data[row*m_ncols + col]; }
  inline T operator() (unsigned int row, unsigned int col) const { return m_data[row*m_ncols + col]; }

  inline unsigned int NRows() const { return m_nrows; }
  inline unsigned int NCols() const { return m_ncols; }
  inline unsigned int NElements() const { return m_ncols*m_nrows; }


  // TODO: I must implement also division and sum.
  inline FFTMatrix<T> operator+ (const FFTMatrix<T> &other) {
  	FFTMatrix<T> ret{m_nrows, m_ncols};
  	for (unsigned int i = 0; i < NElements(); ++i) {
  		ret[i] = (*this)[i] + other[i];
  	}
  	return ret;
  }

  inline void operator+= (const FFTMatrix<T> &other) {
  	for (unsigned int i = 0; i < NElements(); ++i)	{
  		(*this)[i] += other[i];
  	}
  }

  inline FFTMatrix<T> operator/ (const FFTMatrix<T> &other) {
  	FFTMatrix<T> ret{m_nrows, m_ncols};
  	for(unsigned int i = 0; i < NElements(); ++i) {
  		ret[i] = (*this)[i] / other[i];
  	}
  	return ret;
  }

  inline void operator/= (const FFTMatrix<T> &other) {
  	for (unsigned int i = 0; i < NElements(); ++i) {
  		(*this)[i] /= other[i];
  	}
  }

  inline void operator*= (const FFTMatrix<T> &other) { 
  	for (unsigned int i = 0; i < m_data.size(); i++)
  		m_data[i] *= other[i];
  }

  inline FFTMatrix<T> operator* (const FFTMatrix<T> &other) {
  	FFTMatrix<T> ret(m_nrows, m_ncols);
  	for (unsigned int i = 0; i < NElements(); i++) 
  		ret[i] = (*this)[i] * other[i];
  	return ret;
  }

  // The two Functions below have been taken from the original code.
  FFTMatrix<T> ScaleAndShift(unsigned int new_rows, unsigned int new_cols, unsigned int sh_row, unsigned int sh_col) const
  {
  	FFTMatrix ret(new_cols, new_rows);
  	unsigned int min_rows = std::min(m_nrows, new_rows);
  	unsigned int min_cols = std::min(m_ncols, new_cols);

  	for (unsigned int i = 0; i < min_rows; i++) {
  		unsigned int ii = m_nrows - sh_row + i;
  		if (ii >= m_nrows) ii -= m_nrows;
  		for (unsigned int j = 0; j < min_cols; j++) {
  			unsigned int jj = new_cols - sh_col + j;
  			if (jj >= new_cols) jj -= new_cols;
  			ret(ii, jj) = (*this)(i, j);
  		}
  	}

  	return ret;
  }

  FFTMatrix<T> ScaleAndShiftReversed(unsigned int new_rows, unsigned int new_cols, unsigned int sh_row, unsigned int sh_col) const
  {
  	FFTMatrix<T> ret(new_rows, new_cols);
		unsigned int min_rows = std::min(m_nrows, new_rows);
  	unsigned int min_cols = std::min(m_ncols, new_cols);
		
		for (unsigned int i = 0; i < min_rows; i++) {
			unsigned int ii = m_nrows - sh_row + i;
			if (ii >= m_nrows) ii -= m_nrows;
			for (unsigned int j = 0; j < min_cols; j++) {
				unsigned int jj = m_nrows - sh_col + j;
				if (jj >= m_ncols) jj -= m_ncols;
				ret(i, j) = (*this)(ii, jj);
			}
		}
		return ret;
  }

  // Method taken from the original code.
  void OneMinusSelf()
  {
  	for (T &elem : m_data) {
  		elem = T(1) - elem;
  	}
  }

private:
  unsigned int m_nrows;
  unsigned int m_ncols;
  std::vector<T> m_data;
}; 

//  This struct is heavy based on MatrixProfile from the original code.
struct MatrixProfile
{
	MatrixProfile(unsigned int length);

	FFTMatrix<real_type> reflectance;
	FFTMatrix<real_type> transmitance;
	unsigned int length() const { return reflectance.NRows(); }
};

// --------------------------- Multipole Solver -------------------------------------------------------
MultipoleTable ComputeMultipoleDiffusionProfile(const std::vector<MultipoleLayer> &layers, 
	const MultipoleOptions &options);

#endif