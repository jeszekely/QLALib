/**
 * \file "wvfxn.hpp"
 * \author  J. Szekely
 */

#ifndef QLALIB_WVFXN
#define QLALIB_WVFXN

#include "arrays.hpp"

typedef std::complex<double> cplx;

class wvfxn1D : public Array1D<cplx>
{
public:
  double hbar;
  double mass;

  wvfxn1D(const int, const double, const double);
  wvfxn1D(const int, const cplx, const cplx);
  wvfxn1D(const int, const double, const double, const double, const double);
  wvfxn1D(const int, const cplx, const cplx, const double, const double);

  wvfxn1D(const wvfxn1D&);
  wvfxn1D(wvfxn1D&&);
  wvfxn1D& operator=(const wvfxn1D&);

//Wavefunction specific operations
  double getNorm();
  void normalize();
  double flux(const int xx);
  double hb();
  double m();
  template <typename T> wvfxn1D operator|(T &o)
  {
    assert (o.Nx() == nx);
    wvfxn1D out(o.Nx(),xinit,xstep);
    for (int ii = 0; ii < out.Nx(); ii++)
      out(ii) = std::conj(vals[ii])*o(ii);
    return out;
  }
};

class wvfxn2D : public Array2D<cplx>
{
public:
  double hbar;
  double mass1;
  double mass2;

  wvfxn2D(const int, const int, const double, const double);
  wvfxn2D(const int, const int, const cplx, const cplx);
  wvfxn2D(const int, const int, const double, const double, const double, const double, const double);
  wvfxn2D(const int, const int, const cplx, const cplx, const double, const double, const double);

  wvfxn2D(const wvfxn2D&);
  wvfxn2D(wvfxn2D&&);
  wvfxn2D& operator=(const wvfxn2D&);

//Wavefunction specific operations
  double getNorm();
  void normalize();
  double flux_x(const int xx);
  double flux_y(const int yy);
  double hb();
  double m1();
  double m2();
  template <typename T> wvfxn2D operator|(T &o)
  {
    assert (o.Nx() == nx && o.Ny() == ny);
    wvfxn2D out(o.Nx(),o.Ny(),xstep,ystep);
    for (int ii = 0; ii < out.Nx()*out.Ny(); ii++)
      out.vals[ii] = std::conj(vals[ii])*o.vals[ii];
    return out;
  };
};

#endif