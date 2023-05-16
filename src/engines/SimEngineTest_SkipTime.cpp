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

#include <random>

/*
A broken simulation engine that skips 5% of time steps
*/
class SimEngineTest_SkipTime
    : public ISimEngine
{
    std::mt19937 rng;

public:
    std::string Name() const override
    {
        return "SimEngineTest_SkipTime";
    }

    bool IsParallel() const override
    { return false; }

    run_result Run(ISimBox *box, bool modified, unsigned num_steps) override 
    {
        CSimBox *cbox=const_cast<CSimBox*>(box->GetSimBox());
        
        for(unsigned i=0; i<num_steps; i++){
            if((rng()%20)!=0){
               cbox->Evolve();
            }
        }
        return {Supported, {}, num_steps};
    }
};

static bool reg_SimEngineTest_SkipTime = SimEngineBase<SimEngineTest_SkipTime>::Register();
