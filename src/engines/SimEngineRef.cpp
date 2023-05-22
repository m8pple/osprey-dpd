#include "SimDefs.h"


#include <cstdint>
#include <vector>
#include <atomic>
#include <vector>
#include <mutex>
#include "DebugAssert.hpp"
#include <cstdlib>
#include <set>
#include <cstring>

#include "ISimEngine.h"


#include "AbstractBead.h"
#include "ISimBox.h"
#include "xxBase.h"

class SimEngineRef
    : public ISimEngine
{
public:
    std::string Name() const override
    {
        return "SimEngineRef";
    }

    bool IsParallel() const override
    { return false; }

    run_result Run(ISimBox *box, bool modified, unsigned start_sim_time, unsigned num_steps) override 
    {
        CSimBox *cbox=const_cast<CSimBox*>(box->GetSimBox());
        for(unsigned i=0; i<num_steps; i++){
           cbox->Evolve(start_sim_time + i);
           cbox->CNTCellCheck();
        }
        return {Supported, {}, num_steps};
    }
};

static bool reg_SimEngineRef = SimEngineBase<SimEngineRef>::Register();