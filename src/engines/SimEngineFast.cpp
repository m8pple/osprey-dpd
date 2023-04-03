#include <vector>
#include <atomic>
#include <vector>
#include <mutex>
#include <cassert>
#include <cstdlib>
#include <set>
#include <cstring>
#include <unordered_map>

#include "StdAfx.h"
#include "SimDefs.h"
#include "xxBase.h"
#include "AbstractBead.h"

#include "ISimEngine.h"

#include "SimBox.h"

#include "ISimBox.h"


struct SimEngineFast
    : public ISimEngine
{
public:
    std::string Name() const override
    {
        return "SimEngineFast";
    }

    std::string CanSupport(const ISimBox *box) const override
    {
        return {};
    }

    bool IsParallel() const override
    { return false; }

    void Run(ISimBox *box, bool modified, unsigned num_steps) override
    {
        CSimBox *mbox=const_cast<CSimBox*>(box->GetSimBox());
        mbox->EvolveFast(num_steps);
    }
};

static bool reg_SimEngineFast = SimEngineBase<SimEngineFast>::Register();