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
#include "AbstractBead.h"

// TODO: Needs adapting for changed ISimEngine interface
#if 0

/*
A broken simulation engine that doesn't update 1% of the beads.
*/
class SimEngineTest_FrozenBeads
    : public ISimEngine
{
public:
    std::string Name() const override
    {
        return "SimEngineTest_FrozenBeads";
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

        std::vector<CAbstractBead*> frozen;
        auto beads=cbox->GetAllBeadsInCNTCells();
        for(unsigned i=0; i<beads.size(); i++){
            auto b=beads[i];
            if(!b->GetFrozen() && ((i%100)==0)){
                b->SetFrozen();
                frozen.push_back(b);
            }
        }
        
        for(unsigned i=0; i<num_steps; i++){
           cbox->Evolve();
        }

        for(auto b : frozen){
            b->SetNotFrozen();
        }
    }
};

static bool reg_SimEngineTest_SkipTime = SimEngineBase<SimEngineTest_FrozenBeads>::Register();

#endif