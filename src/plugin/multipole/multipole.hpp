#pragma once

#ifndef MULTIPOLE_HPP_INCLUDED
#define MULTIPOLE_HPP_INCLUDED


#include <vector>
#include "external/simple_fft/fft_settings.h"

// --------------------------- Dipole Solver ---------------------------------------------------

class DipoleSolver
{
public:
  DipoleSolver(Float eta0, Float eta1, Float sigma_a, Float sigma_s_prime, int zi, Float thickness);

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
  int zi
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


// ----------------------- Fast Fourier Matrix "adapter" ---------------------------------------------
template<typename T>
class FFTMatrix
{
public:
  FFTMatrix(unsigned int nrows, unsigned int ncols)
  	:m_nrows(nrows), m_ncols(ncols), m_data(m_nrows*m_ncols)
  {}

  inline T& operator[](unsigned int index) { return m_data[index]; }
  inline T operator[] const (unsigned int index) { return m_data[index]; }

  inline T& operator() (unsigned int row, unsigned int col) { return m_data[row*m_ncols + col]; }
  inline T operator() const (unsigned int row, unsigned int col) { return m_data[row*m_ncols + col]; }

  inline unsigned int NRows() const { return m_nrows; }
  inline unsigned int NCols() const { return m_ncols; }
  inline unsigned int NElements() const { return m_ncols*m_nrows; }

  inline Matrix<T> operator*= (const Matrix<T> &other) { 
  	for (unsigned int i = 0; i < m_data.size(); i++)
  		m_data[i] *= other[i];
  }

  inline Matrix<T> operator* (const Matrix<T> &other) {
  	Matrix<T> ret(m_nrows, m_ncols);
  	for (unsigned int i = 0; i < NElements(); i++) 
  		ret[i] = (*this)[i] * other[i];
  	return ret;
  }

  // The two Functions below have been taken from the original code.
  Matrix<T> ScaleAndShift(unsigned int new_rows, unsigned int new_cols, unsigned int sh_row, unsigned int sh_col) const
  {
  	Matrix ret(new_cols, new_rows);
  	unsigned int min_rows = min(m_nrows, new_rows);
  	unsigned int min_cols = min(m_ncols, new_cols);

  	for (unsigned int i = 0; i < min_rows; i++) {
  		unsigned int ii = m_nrows - sh_row + i;
  		if (ii >= m_nrows) ii -= m_nrows;
  		for (unsigned int j = ; j < min_cols j++) {
  			unsigned int jj = new_cols - sh_col + j;
  			if (jj >= new_cols) jj -= new_cols;
  			ret(ii, jj) = (*this)(i, j);
  		}
  	}

  	return ret;
  }

  Matrix<T> ScaleAndShiftReversed(unsigned int new_rows, unsigned int new_cols, unsigned int sh_row, unsigned int sh_col) const
  {
  	Matrix<T> ret(new_rows, new_cols);
		unsigned int min_rows = min(m_nrows, new_rows);
  	unsigned int min_cols = min(m_ncols, new_cols);
		
		for (unsigned int i = 0; i < min_rows; i++) {
			unsigned int ii = m_nrows - sh_row + i;
			if (ii >= new_rows) ii -= m_nrows;
			for (unsigned int j = 0; j < min_cols; j++) {
				unsigned int jj = m_nrows - sh_col + j;
				if (jj >= m_ncols) jj -= m_ncols;
				ret(i, j) = (*this)(ii, jj);
			}
		}
		return ret;
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

	FFTMatrix reflectance;
	FFTMatrix transmitance;
	unsigned int length() { return reflectance.length(); }
};

// --------------------------- Multipole Solver -------------------------------------------------------
MultipoleTable ComputeMultipoleDiffusionProfile(const std::vector<MultipoleLayer> &layers, 
	const MultipoleOptions &options);

#endif