#include <module/Module.h>
#include <distributions/DVonMises.h>
#include <functions/DLogVonMisesFun.h>

using std::vector;

namespace jags {
namespace vonmises {

class VONMISESModule : public Module {
  public:
    VONMISESModule();
    ~VONMISESModule();
};

VONMISESModule::VONMISESModule() : Module("vonmises")
{
  //load distributions
  insert(new DVonMises);
  insert(new DLogVonMisesFun);
}

VONMISESModule::~VONMISESModule() 
{
  vector<Function*> const &fvec = functions();
  for (unsigned int i = 0; i < fvec.size(); ++i) {
    delete fvec[i];
  }
  vector<Distribution*> const &dvec = distributions();
  for (unsigned int i = 0; i < dvec.size(); ++i) {
    delete dvec[i];
  }
}

} // namespace vonmises
} // namespace jags

jags::vonmises::VONMISESModule _vonmises_module;
