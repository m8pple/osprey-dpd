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
A correct simulation engine that adds position sizzle with similar std-deviation
to a few single-precision rounding errors.
*/
class SimEngineTest_Sizzle
    : public IIntegrationEngine
{
public:
    std::string Name() const override
    {
        return "SimEngineTest_Sizzle";
    }

    bool IsParallel() const override
    { return false; }

    run_result Run(ISimBox *box, bool modified, unsigned start_sim_time, unsigned num_steps) override 
    {
        CSimBox *cbox=const_cast<CSimBox*>(box->GetSimBox()); 

        auto beads=cbox->GetAllBeadsInCNTCells();

        std::mt19937_64 rng;
        std::uniform_real_distribution<> dist(-3e-8, 3e-8);

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

                    b->SetXMom( b->GetXMom() * (1+dist(rng)) );
                    b->SetYMom( b->GetYMom() * (1+dist(rng)) );
                    b->SetZMom( b->GetZMom() * (1+dist(rng)) );

                    b->SetXForce( b->GetXForce() * (1+dist(rng)) );
                    b->SetYForce( b->GetYForce() * (1+dist(rng)) );
                    b->SetZForce( b->GetZForce() * (1+dist(rng)) );


                }
            }
            cbox->Evolve(start_sim_time+i);
        }

        return {Supported, {}, num_steps};
    }
};

static bool reg_SimEngineTest_WrongDiss = SimEngineBase<SimEngineTest_Sizzle>::Register();
