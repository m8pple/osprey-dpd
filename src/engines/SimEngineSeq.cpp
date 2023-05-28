#include "SimEngineSeq.h"


namespace {
    struct EnginePolicyRNG
        : EnginePolicyConcept
    {
        static const RNGPolicy RNG_POLICY = RNGPolicy_AsymmetricRng;
    };

    struct EnginePolicyHashBeadTag
        : EnginePolicyConcept
    {
        static const RNGPolicy RNG_POLICY = RNGPolicy_HashBeadTag;
    };

    struct EnginePolicyRNGMorton
        : EnginePolicyRNG
    {
        static const bool USE_MORTON = true;
    };

    struct EnginePolicyHashBeadTagMorton
        : EnginePolicyHashBeadTag
    {
        static const bool USE_MORTON = true;
    };

    static bool reg_SimEngineSeqRNG = SimEngineBase<SimEngineSeq<EnginePolicyRNG>>::Register();
    static bool reg_SimEngineSeqRNGMorton = SimEngineBase<SimEngineSeq<EnginePolicyRNGMorton>>::Register();
    static bool reg_SimEngineSeqHashBeadTag = SimEngineBase<SimEngineSeq<EnginePolicyHashBeadTag>>::Register();
    static bool reg_SimEngineSeqHashBeadTagMorton = SimEngineBase<SimEngineSeq<EnginePolicyHashBeadTagMorton>>::Register();
};
