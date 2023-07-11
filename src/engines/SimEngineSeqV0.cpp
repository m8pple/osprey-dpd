#include "SimEngineSeqV0.h"

namespace {

    static bool reg_SimEngineSeqV0_C32_A32_P32 = SimEngineBase<SimEngineSeqV0<float,float,float>>::Register();
    static bool reg_SimEngineSeqV0_C32_A32_P64 = SimEngineBase<SimEngineSeqV0<float,float,double>>::Register();
    static bool reg_SimEngineSeqV0_C32_A64_P32 = SimEngineBase<SimEngineSeqV0<float,double,float>>::Register();
    static bool reg_SimEngineSeqV0_C32_A64_P64 = SimEngineBase<SimEngineSeqV0<float,double,double>>::Register();
    static bool reg_SimEngineSeqV0_C64_A64_P64 = SimEngineBase<SimEngineSeqV0<double,double,double>>::Register();
};