#include <config.h>
#include "DVonMises.h"

#include <rng/RNG.h>
#include <util/nainf.h>
#include <cmath>
#include <JRmath.h>

using std::vector;
using std::log;
using std::string;
using std::cos;
using std::acos;
using std::sqrt;
using std::floor;
using std::fmod;

namespace jags {
namespace vonmises {

DVonMises::DVonMises() : ScalarDist("dvonmises", 2, DIST_SPECIAL)
{}

DVonMises::DVonMises(string const &name, unsigned int npar) : ScalarDist(name, npar, DIST_SPECIAL)
{}

double DVonMises::logDensity(double x, PDFType type,
  vector<double const *> const &par,
  double const *lower, double const *upper) const
{

  if(x < 0 || x > _2pi()){
    return JAGS_NEGINF;
  }
  else if(type == PDF_PRIOR){
  //omit terms that depend only on parameters
    return *par[1]*cos(x-*par[0]);
  }
  else{
  //full log density
    return *par[1]*cos(x-*par[0]) - _log2pi() - log(modBessel0(*par[1]));
  }
  
}


double DVonMises::l(std::vector<double const *> const &parameters) const
{
  return 0;
}

double DVonMises::u(std::vector<double const *> const &parameters) const
{
  return _2pi();
}


double DVonMises::randomSample(vector<double const *> const &par,
        double const *lower, double const *upper, RNG *rng) const
{

  double theta = 0;
  if(*par[1]==0){
  //if concentration is zero, von Mises becomes a uniform distribution on [0, 2pi)
    return _2pi()*(rng->uniform());
  }
  else if(*par[1]<=0.5){
  //if concentration is small, do rejection sampling with uniform envelope
  //a.k.a. "Siegerstetter method"
    while(true){
      double u1 = _2pi()*(rng->uniform());
      double u2 = rng->uniform();
      if(log(u2) <= (*par[1])*(cos(u1)-1)){
        theta = u1;
        break;
      }
    }
    return fmod(theta + *par[0], _2pi());
    //add the mean, and take result mod 2pi
  }
  else{
  //if concentration is large, do rejection sampling with wrapped Cauchy envelope
  //a.k.a. "Best-Fisher method"
    double k   = *par[1];
    double tau = 1 + sqrt(1 + 4*k*k);
    double rho = (tau - sqrt(2*tau))/(2*k);
    double r   = (1+rho*rho)/(2*rho);
    double theta = 0;
    while(true){
      double z  = cos(M_PI*rng->uniform());
      double f  = (1 + r*z)/(r + z);
      double c  = k*(r - f);
      double u2 = rng->uniform();
      if((c*(2-c) - u2) > 0 || (log(c/u2) + 1 - c) > 0){
      //acceptance pre-check to avoid evaluating log
      //(uses || operator short-circuiting)
        double u3 = rng->uniform();
        theta = (floor(u3+0.5)*2 - 1)*acos(f);
        break;
      }
    }
    //the above algorithm generates a random variate on [-pi, pi]
    //but the modulus takes care of that
    //just need to add the mean 
    return fmod(theta + *par[0], _2pi());
  }
}

double DVonMises::typicalValue(vector<double const *> const &par,
  double const *lower, double const *upper) const
{
  return *par[0];
  //just return the mean value
}

bool DVonMises::checkParameterValue (vector<double const *> const &par) const
{
  return (*par[0]>=0) && (*par[0]<=_2pi()) && (*par[1]>=0);
}

bool DVonMises::isSupportFixed(vector<bool> const &fixmask) const
{
  return fixmask[0] && fixmask[1];
}

//compute the modified Bessel function of first kind, of order zero
//this is equal to \sum_{i=0}^\infty ((x/2)^{2i})/((i!)^2)
double DVonMises::modBessel0(double x) const
{
  double ccrit = 10E-6;
  //convergence criterion
  //stop when the increments in the series are smaller than this
  //the precise value does not affect speed much, in general
  //because the increments eventually start to decrease exponentially

  double s = 1;
  int i = 1;
  double inc = 1;
  double x_2i = 0;
  while (true){
    x_2i = x/(2*i);
    inc  = inc * x_2i * x_2i;
    s    = s + inc;
    i    = i + 1;
    if(inc < ccrit){
      break;
    }
  }
  return s;
}

} //namespace vonmises
} //namespace jags
