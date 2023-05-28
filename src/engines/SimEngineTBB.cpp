#include "SimEngineTBB.h"


namespace{
    struct EnginePolicyTBBMorton
        : EnginePolicyTBB
    {
        static const bool USE_MORTON = true;
    };

    static bool reg_SimEngineTBB = SimEngineBase<SimEngineTBB<>>::Register();
    static bool reg_SimEngineTBBMorton = SimEngineBase<SimEngineTBB<EnginePolicyTBBMorton>>::Register();
}
