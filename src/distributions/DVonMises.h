/*  The von Mises distribution ("circular equivalent of the normal distribution")
  For a reference for this distribution, see e.g. 

  Forbes, C., Evans, M.,  Hastings, N. and Peacock, B. (2011). Statistical Distributions, 4th       edition (Hoboken: John Wiley and Sons)

  Reference for the sampling method:

  Best, D.J. and Fisher, N.I. (1979). Efficient Simulation of the von Mises distribution. 
    J. Roy. Statist. Soc. Ser. C 28, 152-157
**/

#ifndef DVONMISES_H_
#define DVONMISES_H_

#include <distribution/ScalarDist.h>

using std::vector;

namespace jags {
namespace vonmises {

class DVonMises : public ScalarDist 
{
  public:
  DVonMises();
  DVonMises(std::string const &name, unsigned int npar);

  double   logDensity(double x, PDFType type,
    vector<double const *> const &parameters,
    double const *lower, double const *upper) const;
  double   randomSample(vector<double const *> const &parameters,
    double const *lower, double const *upper, RNG *rng) const;
  double   typicalValue(vector<double const *> const &parameters,
    double const *lower, double const *upper) const;
   bool   checkParameterValue(vector<double const *> const &parameters) const;

  double   l(vector<double const *> const &parameters) const;
  double   u(vector<double const *> const &parameters) const;

  bool  isSupportFixed(vector<bool> const &fixmask) const;

  private:
  double modBessel0(double x) const;
  inline double _2pi() const
    {return 6.283185307179586;}
  inline double _log2pi() const
    {return 1.837877066409345;}
};

} //namespace vonmises
} //namespace jags

#endif /* DVONMISES_H_ */
