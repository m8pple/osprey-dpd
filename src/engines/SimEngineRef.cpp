#include "SimDefs.h"


#include <cstdint>
#include <vector>
#include <atomic>
#include <vector>
#include <mutex>
#include <cassert>
#include <cstdlib>
#include <set>
#include <cstring>

#include "ISimEngine.h"


#include "AbstractBead.h"
#include "ISimBox.h"
#include "xxBase.h"

struct SimEngineRef
    : public ISimEngine
{
public:
    std::string Name() const override
    {
        return "SimEngineRef";
    }

    std::string CanSupport(const ISimBox *box) const override
    {
        return {};
    }

    bool IsParallel() const override
    { return false; }

    void Run(ISimBox *box, bool modified, unsigned num_steps) override 
    {
        CSimBox *cbox=const_cast<CSimBox*>(box->GetSimBox());
        for(unsigned i=0; i<num_steps; i++){
           cbox->Evolve();
        }
    }
};

static bool reg_SimEngineRef = SimEngineBase<SimEngineRef>::Register();