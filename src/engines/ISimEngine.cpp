#include "ISimEngine.h"

#include "CNTCell.h"
#include "ISimBox.h"
#include "aeActiveSimBox.h"

ISimEngineCapabilities::support_result ISimEngineCapabilities::StandardSupportAssumptions(const ISimBox *box)
{
    #if (SimDimension!=3)
    return {PermanentProblem,"Compile-time standard SimEngine assumption not met - dimensions!=3."};
    #endif
    #if (EnableStressTensorSphere == SimMiscEnabled)
    return {PermanentProblem,"Compile-time standard SimEngine assumption not met - EnableStressTensorSphere == SimMiscEnabled."};
    #endif
    #if (EnableParallelSimBox != SimMPSDisabled)
    return {PermanentProblem,"Compile-time standard SimEngine assumption not met - EnableParallelSimBox != SimMPSDisabled."};
    #endif
    #if defined(UseDPDBeadRadii)
    // This could be relaxed if needed. Not too difficult to support, but adds a bit of programming/optimisaiton overhead
    return {PermanentProblem,"Compile-time standard SimEngine assumption not met - UseDPDBeadRadii is defined."};        
    #endif

    if(CCNTCell::GetLambda()!=0.5){
        return {PermanentProblem,"Run-time SimBox standard SimEngine assumption not met - Lambda!=0.5"};
    }

    if(box->GetSimBoxXOrigin()!=0 || box->GetSimBoxXOrigin()!=0 || box->GetSimBoxZOrigin()!=0){
        return {PermanentProblem, "Run-time SimBox standard SimEngine assumption not met - Sim box origin is not at zero."}; // Actually fairly easy to support
    }

    if( fmod( box->GetSimBoxXLength(), 1) != 0  || fmod(box->GetSimBoxYLength(), 1) != 0 || fmod(box->GetSimBoxZLength(), 1) != 0){
        return {PermanentProblem,"Run-time SimBox standard SimEngine assumption not met - Box dimensions are not integer"};
    }
    if( box->GetCNTXCellWidth() != 1 || box->GetCNTYCellWidth()!=1 || box->GetCNTZCellWidth()!=1){
        return {PermanentProblem,"Run-time SimBox standard SimEngine assumption not met - Cell dimensions are not integer"};
    }

    CSimBox *mbox=const_cast<CSimBox*>(box->GetSimBox());

    if(mbox->IsEnergyOutputOn()){
        return {TransientProblemFlags,"Run-time SimBox standard SimEngine assumption not met - IsEnergyOutputOn"};
    }
    if(mbox->IsActiveBondsOn() && mbox->m_pShadow && mbox->m_pShadow->IsAnyACNPresent()){
        return {TransientProblemFlags,"Run-time SimBox standard SimEngine assumption not met - IsActiveBondsOn and active ACNs"};
    }
    if(mbox->IsBeadChargeOn()){
        return {TransientProblemFlags,"Run-time SimBox standard SimEngine assumption not met - IsBeadChargeOn"};
    }
    if(mbox->IsGravityOn()){
        return {TransientProblemFlags,"Run-time SimBox standard SimEngine assumption not met - IsGravityOn"};
    }

    if(!mbox->m_ActiveCommandTargets.empty()){
        return {TransientProblemStep,"Run-time SimBox standard SimEngine assumption not met - ActiveCommandTargets"};
    }
    if(!mbox->m_ActiveForceTargets.empty()){
        return {TransientProblemStep,"Run-time SimBox standard SimEngine assumption not met - ActiveForceTargets"};
    }


    return {Supported,{}};
}

ISimEngineCapabilities::support_result ISimEngineCapabilities::CanSupport(const ISimBox *box) const
{
    auto tmp=StandardSupportAssumptions(box);
    if(tmp.status!=Supported){
        return tmp;
    }
    return CanSupport_ExtraConstraints(box);
}