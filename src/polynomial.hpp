#ifndef QLALIB_POLYNOMIAL
#define QLALIB_POLYNOMIAL

#include <gsl/gsl_sf_bessel.h>
#include "arrays.hpp"

typedef std::complex<double> cplx;

template <typename T>
class polynomial
{
public:
  std::shared_ptr<std::vector<T>> vals;

  polynomial(const int n) : vals(std::make_shared<std::vector<T>>(n))
  {
    std::fill_n(vals->data(),size(),T(0.0));
  }
  polynomial(const polynomial& o) : vals(std::make_shared<std::vector<T>>(o.size()))
  {
    std::copy_n(o.vals->data(), size(), vals->data());
  }
  // polynomial(polynomial&& o) : size()(o.size()), vals(std::move(o.vals)){}


  T& element(const int nn)
  {
    return vals->at(nn);
  }

  const T& element(const int nn) const
  {
    return vals->at(nn);
  }

  T& operator()(const int nn)
  {
    return element(nn);
  }

  const T& operator()(const int nn) const
  {
    return element(nn);
  }

  const int size() const { return vals->size(); }

  polynomial& operator=(const polynomial<T>& o)
  {
    vals = std::make_shared<std::vector<T>>(o.size());
    copy_n(&o(0),o.size(),&element(0));
    return *this;
  }

  polynomial derivative()
  {
    int len = size()-1;
    polynomial<T> out(len);
    for (int ii = 0; ii < len; ii++)
      out(ii) = double(ii+1)*vals->at(ii+1);
    return out;
  }

  /// Finds all roots of the polynomial
  std::shared_ptr<std::vector<cplx>> findRoots()
  {
    const int N = size()-1;
    const cplx I = cplx(0.0,1.0);
    int maxIt = 100;

    polynomial<T> polynomialCopy (*this);
    cplx coeff (polynomialCopy(polynomialCopy.size()-1));
    if (real(coeff) != 1.0) polynomialCopy.scale(1.0/coeff);
    std::cout << polynomialCopy;
    std::vector<cplx> guesses(N,cplx(0.0));
    std::vector<cplx> prevIteration(N,cplx(0.0));
    for (int ii = 0; ii < N; ii++)
        guesses[ii] = ( double(ii)/double(N) + (1.0/(3.0*double(N)))) * exp(2.0*I*M_PI*static_cast<double>(ii)/static_cast<double>(N)) ; //Very sensitive to initial guess
    double converged;
    int iteration = 0;
    do
    {
        converged = 0.0;
        std::copy_n(guesses.begin(), N, prevIteration.begin());
        for (int ii = 0; ii < N; ii++)
        {
            cplx num = polynomialCopy.eval(prevIteration[ii]);
            for (int jj = 0; jj < N; jj++)
            {
                if (ii != jj) num /= (prevIteration[ii] - guesses[jj]);
            }
            guesses[ii] -= num;
        }
        for (int ii = 0; ii < N; ii++)
            converged = std::max(converged,abs(guesses[ii]-prevIteration[ii]));
        iteration++;
        if (iteration >= maxIt) std::cout << "Max Iterations Reached" << std::endl;
    } while (converged > 0.0001 && iteration < maxIt);
    std::cout << "Zeros Found: " << std::endl;
    std::sort(guesses.begin(), guesses.end(),[](cplx a, cplx b){return real(a) < real(b);});
    for (int ii = 0; ii < N; ii++)
        std::cout << std::setprecision(10) << real(guesses[ii]) << " " << abs(polynomialCopy.eval(guesses[ii])) <<  std::endl;
    std::shared_ptr<std::vector<cplx>> v = std::make_shared<std::vector<cplx>>(guesses);
    return v;
  }

  std::shared_ptr<std::vector<cplx>> findRootsAlt()
  {
    const int N = size()-1;
    // const cplx I = cplx(0.0,1.0);
    // int maxIt = 500;

    polynomial<T> polynomialCopy (*this);
    cplx coeff (polynomialCopy(polynomialCopy.size()-1));
    if (real(coeff) != 1.0) polynomialCopy.scale(1.0/coeff);
    std::cout << polynomialCopy;
    // std::vector<cplx> guesses(N,cplx(0.0));
    std::vector<cplx> guesses( *findRoots() );
    std::vector<cplx> prevIteration(N,cplx(0.0));
//  Initial guess currently optimized for Legendre polynomials
    double c = 1.0-pow(2.0/M_PI,2);
    // for (int ii = 0; ii < N; ii++)
    // {
      // double jk = gsl_sf_bessel_zero_J0(ii+1);
      // double denom = sqrt( pow(double(1+ii)+0.5,2)+0.25*c );
      // double lowerbound = jk/denom;
      // double upperbound = jk/(double(ii+1) + 0.5);
      // guesses[ii] = cos(lowerbound);
      // std::cout << jk << " " << guesses[ii] << " " << lowerbound << " " << upperbound << std::endl;
    // }

    polynomial deriv = derivative();
    int maxIts = 1000;
    for (int ii = 0; ii < N; ii++)
    {
      int its = 0;
      cplx currentIteration;
      bool notConverged = true;
      while (notConverged)
      {
        currentIteration = guesses[ii] - eval(guesses[ii])/deriv.eval(guesses[ii]);
        if (abs(eval(currentIteration)) < 1.0e-8 || its > maxIts)
        {
          notConverged = false;
          std::cout << its << std::endl;
        }
        guesses[ii] = currentIteration;
        its++;
      }
    }
    std::cout << "Zeros Found: " << std::endl;
    std::sort(guesses.begin(), guesses.end(),[](cplx a, cplx b){return real(a) < real(b);});
    for (int ii = 0; ii < N; ii++)
        std::cout << std::setprecision(10) << real(guesses[ii]) << " " << abs(polynomialCopy.eval(guesses[ii])) <<  std::endl;
    std::shared_ptr<std::vector<cplx>> v = std::make_shared<std::vector<cplx>>(guesses);
    return v;
  }

  polynomial operator+(const polynomial<T>& o) const
  {
    int len;
    const polynomial<T> *p;
    size() >= o.size() ? (p=this, len=size()) : (p=&o, len=o.size());
    polynomial<T> out(len);
    std::transform(vals->data(),vals->data()+std::min(size(),o.size()),o.vals->data(),out.vals->data(),[](T a, T b){return a+b;});
    if (size() != o.size()) std::copy_n(&p->element(std::min(size(),o.size())),abs(size()-o.size()),&out(std::min(size(),o.size())));
    return out;
  }

  polynomial& operator+=(const polynomial<T>& o)
  {
    *this = *this + o;
    return *this;
  }

  polynomial operator-(const polynomial<T>& o) const
  {
    int len, sign;
    int min = std::min(size(),o.size());
    const polynomial<T> *p;
    size() >= o.size() ? (p=this, len=size(), sign=1) : (p=&o, len=o.size(), sign=-1);
    polynomial<T> out(len);
    std::transform(vals->data(),vals->data()+min,o.vals->data(),out.vals->data(),[](T a, T b){return a-b;});
    if (size() != o.size())
    {
      std::copy_n(&p->element(min),abs(size()-o.size()),&out(min));
      if (sign == -1) std::transform(&out(min),&out(min)+abs(size()-o.size()),&out(min),[](T a){return -1.0*a;});
    }
    return out;
  }

  polynomial& operator-=(const polynomial<T>& o)
  {
    *this = *this - o;
    return *this;
  }

  polynomial operator*(const polynomial<T>& o) const
  {
    int len = size()+o.size()-1;
    polynomial<T> out(len);
    for (int ii = 0; ii<size(); ii++)
    {
      for (int jj = 0; jj<o.size(); jj++)
        out(ii+jj) += (element(ii)*o(jj));
    }
    return out;
  }

  polynomial& operator*=(const polynomial<T>& o)
  {
    *this = *this * o;
    return *this;
  }

  template <typename U>
  void scale(const U a)
  {
    std::for_each(vals->data(), vals->data()+size(), [&a](T& p){p*=a;});
    return;
  }

  Array1D<T> eval(Array1D<T> &o)
  {
    Array1D<T> out(o);
    std::fill_n(out.data(),out.size(),element(0));
    Array1D<T> power(o);
    for (int ii = 1; ii<size(); ii++)
    {
      out += (power*element(ii));
      power *= o;
    }
    return out;
  }

  template <typename U>
  U eval(U a)
  {
    U out = U(element(0));
    for (int ii = 1; ii<size(); ii++)
      out += std::pow(a,ii)*element(ii);
    return out;
  }

};

template <typename T>
std::ostream &operator<<(std::ostream &out, const polynomial<T> &o)
{
  for (int ii = 0; ii < o.size(); ii++)
    out << o(ii) << "\t";
  out << std::endl;
  return out;
};

//"probabilist's" Hermite polynomials
template <typename T>
std::shared_ptr<polynomial<T>> Hermite(int nn)
{
  std::shared_ptr<polynomial<T>> p0,p1,pn;
  polynomial<T> a1(2);
  a1(1)          = 1;
  p0             = std::make_shared<polynomial<T>>(1);
  p0->element(0) = T(1.0);
  p1             = std::make_shared<polynomial<T>>(2);
  p1->element(1) = T(1.0);
  if (nn == 0) return p0;
  if (nn == 1) return p1;
  int kk = 1;
  while (kk < nn)
  {
    kk++;
    p0->scale(kk-1);
    pn = std::make_shared<polynomial<T>>(*p1*a1-(*p0));
    p0 = p1;
    p1 = pn;
  }
  return pn;
};


//"physicists's" Hermite polynomials
template <typename T>
std::shared_ptr<polynomial<T>> physHermite(int nn)
{
  std::shared_ptr<polynomial<T>> p0,p1,pn;
  polynomial<T> a1(2);
  a1(1)          = 2;
  p0             = std::make_shared<polynomial<T>>(1);
  p0->element(0) = T(1.0);
  p1             = std::make_shared<polynomial<T>>(2);
  p1->element(1) = T(2.0);
  if (nn == 0) return p0;
  if (nn == 1) return p1;
  int kk = 1;
  while (kk < nn)
  {
    kk++;
    p0->scale(2.0*(kk-1));
    pn = std::make_shared<polynomial<T>>(*p1*a1-(*p0));
    p0 = p1;
    p1 = pn;
  }
  return pn;
};

template <typename T>
std::shared_ptr<polynomial<T>> Legendre(int nn)
{
  std::shared_ptr<polynomial<T>> p0,p1,pn;
  polynomial<T> a1(2);
  a1(1)          = T(1.0);
  p0             = std::make_shared<polynomial<T>>(1);
  p0->element(0) = T(1.0);
  p1             = std::make_shared<polynomial<T>>(2);
  p1->element(1) = T(1.0);
  if (nn == 0) return p0;
  if (nn == 1) return p1;
  int kk = 1;
  while (kk < nn)
  {
    p0->scale(kk);
    a1(1) = (2*kk+1);
    pn = std::make_shared<polynomial<T>>(*p1*a1-(*p0));
    p0 = p1;
    p1 = pn;
    pn->scale(1.0/(static_cast<T>(kk+1)));

    kk++;
  }
  return pn;
};

//Chebyshev polynomials of the first kind
template <typename T>
std::shared_ptr<polynomial<T>> Chebyshev(int nn)
{
  std::shared_ptr<polynomial<T>> p0,p1,pn;
  polynomial<T> a1(2);
  a1(1)          = 2;
  p0             = std::make_shared<polynomial<T>>(1);
  p0->element(0) = T(1.0);
  p1             = std::make_shared<polynomial<T>>(2);
  p1->element(1) = T(1.0);
  if (nn == 0) return p0;
  if (nn == 1) return p1;
  int kk = 1;
  while (kk < nn)
  {
    kk++;
    pn = std::make_shared<polynomial<T>>(*p1*a1-(*p0));
    p0 = p1;
    p1 = pn;
  }
  return pn;
};

#endif