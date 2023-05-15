#include "SimEngineTBB.hpp"

// TODO: Needs adapting for changed ISimEngine interface
#if 0


namespace{
    struct EnginePolicyTBBMorton
        : EnginePolicyTBB
    {
        static const bool USE_MORTON = true;
    };

    static bool reg_SimEngineTBB = SimEngineBase<SimEngineTBB<>>::Register();
    static bool reg_SimEngineTBBMorton = SimEngineBase<SimEngineTBB<EnginePolicyTBBMorton>>::Register();
}

#endif
