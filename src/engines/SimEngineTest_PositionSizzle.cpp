#include "SimDefs.h"


#include <cstdint>
#include <vector>
#include <atomic>
#include <vector>
#include <mutex>
#include "DebugAssert.h"
#include <cstdlib>
#include <set>
#include <cstring>
#include <cfloat>

#include "IIntegrationEngine.h"


#include "AbstractBead.h"
#include "ISimBox.h"
#include "xxBase.h"
#include "AbstractBead.h"
#include "CNTCell.h"

#include <random>

/*
A correct simulation engine that adds position sizzle with same std-deviation
as single-precision rounding error.
*/
class SimEngineTest_PositionSizzle
    : public IIntegrationEngine
{
public:
    std::string Name() const override
    {
        return "SimEngineTest_PositionSizzle";
    }

    bool IsParallel() const override
    { return false; }

    run_result Run(ISimBox *box, bool modified, unsigned start_sim_time, unsigned num_steps) override 
    {
        CSimBox *cbox=const_cast<CSimBox*>(box->GetSimBox()); 

        auto beads=cbox->GetAllBeadsInCNTCells();

        std::mt19937_64 rng;
        std::uniform_real_distribution<> dist(-1e-8, 1e-8);

        double lower[3]={box->GetSimBoxXOrigin(), box->GetSimBoxYOrigin(), box->GetSimBoxZOrigin() };
        double upper[3]={ std::nexttoward(box->GetSimBoxXLength(), -DBL_MAX), std::nexttoward(box->GetSimBoxYLength(),-DBL_MAX), std::nexttoward(box->GetSimBoxZLength(), -DBL_MAX) };

        for(unsigned i=0; i<num_steps; i++){
            for(CAbstractBead *b : beads){
                if(!b->GetFrozen()){
                    double dx=dist(rng), dy=dist(rng), dz=dist(rng);
                    b->SetXPos( std::max(lower[0], std::min(upper[0], b->GetXPos()+dx)));
                    b->SetYPos( std::max(lower[1], std::min(upper[1], b->GetYPos()+dy)));
                    b->SetZPos( std::max(lower[2], std::min(upper[2], b->GetZPos()+dz)));

                    b->SetunPBCXPos(b->GetunPBCXPos() + dx);
                    b->SetunPBCYPos(b->GetunPBCYPos() + dy);
                    b->SetunPBCZPos(b->GetunPBCZPos() + dz);
                    
                }
            }
            cbox->Evolve(start_sim_time+i);
        }

        return {Supported, {}, num_steps};
    }
};

static bool reg_SimEngineTest_WrongDiss = SimEngineBase<SimEngineTest_PositionSizzle>::Register();
