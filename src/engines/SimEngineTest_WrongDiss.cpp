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
#include "CNTCell.h"

/*
A broken simulation engine that uses the wrong dissipative parameter.
*/
class SimEngineTest_WrongDiss
    : public ISimEngine
{
public:
    std::string Name() const override
    {
        return "SimEngineTest_WrongDiss";
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

        unsigned bead_types=cbox->GetBeadTypeTotal();
        std::vector<std::vector<double>> original(bead_types);
        for(unsigned i=0; i<bead_types; i++){
            original[i].resize(bead_types);
            for(unsigned j=0; j<=i; j++){
                double v=CCNTCell::GetDPDBeadConsInt(i,j);
                original[i][j]=v;
                if(i==j){
                    v *= 0.95;
                }else{
                    v *= 1.05;
                }
                CCNTCell::SetDPDBeadConsInt(i, j, v);
            }
        }
        

        for(unsigned i=0; i<num_steps; i++){
           cbox->Evolve();
        }

        for(unsigned i=0; i<bead_types; i++){
            for(unsigned j=0; j<=i; j++){
                CCNTCell::SetDPDBeadConsInt(i,j, original[i][j]);
            }
        }
    }
};

static bool reg_SimEngineTest_WrongDiss = SimEngineBase<SimEngineTest_WrongDiss>::Register();
