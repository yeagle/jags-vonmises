#include <config.h>
#include "DLogVonMisesFun.h"

#include <util/nainf.h>
#include <cmath>
#include <JRmath.h>

using std::vector;
using std::log;
using std::cos;

namespace jags {
namespace vonmises {

DLogVonMisesFun::DLogVonMisesFun() : ScalarFunction("dlogvonmises", 3)
{}

bool DLogVonMisesFun::checkParameterValue (vector<double const *> const &args) const
{
  return (*args[0]>=0) && (*args[0]<=_2pi()) && (*args[1]>=0);
}

double DLogVonMisesFun::evaluate(vector<double const *> const &args) const
{
  if(*args[0] < 0 || *args[0] > _2pi()){
    return JAGS_NEGINF;
  }
  else{
    //full log density
    return *args[2]*cos(*args[0]-*args[1]) - _log2pi() - log(modBessel0(*args[2]));
  }
}

//compute the modified Bessel function of first kind, of order zero
//this is equal to \sum_{i=0}^\infty ((x/2)^{2i})/((i!)^2)
double DLogVonMisesFun::modBessel0(double x) const
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
