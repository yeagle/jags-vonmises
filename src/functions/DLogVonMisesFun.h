/*  The von Mises distribution ("circular equivalent of the normal distribution")
  For a reference for this distribution, see e.g. 

  Forbes, C., Evans, M.,  Hastings, N. and Peacock, B. (2011). Statistical Distributions, 4th       edition (Hoboken: John Wiley and Sons)

  Reference for the sampling method:

  Best, D.J. and Fisher, N.I. (1979). Efficient Simulation of the von Mises distribution. 
    J. Roy. Statist. Soc. Ser. C 28, 152-157
**/

#ifndef DLOGVONMISES_FUNC_H_
#define DLOGVONMISES_FUNC_H_

#include <function/ScalarFunction.h>

using std::vector;

namespace jags {
namespace vonmises {

class DLogVonMisesFun : public ScalarFunction
{
  public:
    DLogVonMisesFun();

    bool checkParameterValue(vector<double const *> const &arguments) const;
    double evaluate(vector<double const *> const &arguments) const;

  private:
    double modBessel0(double x) const;
    inline double _2pi() const
      {return 6.283185307179586;}
    inline double _log2pi() const
      {return 1.837877066409345;}
};

} //namespace vonmises
} //namespace jags

#endif /* DLOGVONMISES_FUNC_H_ */
