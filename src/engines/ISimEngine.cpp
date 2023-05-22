#include "ISimEngine.h"

#include "CNTCell.h"
#include "ISimBox.h"
#include "aeActiveSimBox.h"

#include "bead_id_hash_rng.hpp"
#include "StateLogger.hpp"

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



class SimEngineRefDiff
    : public ISimEngine
{
public:
    SimEngineRefDiff(std::shared_ptr<ISimEngine> testee)
        : m_testee(testee)
    {}

    std::string Name() const override
    {
        return "SimEngineRefDiff_"+m_testee->Name();
    }

    bool IsParallel() const override
    { return m_testee->IsParallel(); }

    run_result Run(ISimBox *box, bool modified, unsigned start_sim_time, unsigned num_steps) override 
    {
        auto check_pbc_drift = [&](const CAbstractBead *bg){
            double eg[3]={
                std::abs(fmod(bg->GetunPBCXPos(), box->GetSimBoxXLength()) - bg->GetXPos()),
                std::abs(fmod(bg->GetunPBCYPos(), box->GetSimBoxYLength()) - bg->GetYPos()),
                std::abs(fmod(bg->GetunPBCZPos(), box->GetSimBoxZLength()) - bg->GetZPos())
            };
            eg[0] = std::min( eg[0], box->GetSimBoxXLength() - eg[0]);
            eg[1] = std::min( eg[1], box->GetSimBoxYLength() - eg[1]);
            eg[2] = std::min( eg[2], box->GetSimBoxZLength() - eg[2]);

            if(eg[0] > 1e-6 || eg[1] > 1e-6 || eg[2] > 1e-6){
                fprintf(stderr, "PBC drift.\n");
                exit(1);
            }
        };

        //std::ofstream dump("dump.txt");
        std::ostream &dst=std::cout;

        std::vector<std::string> log_lines;

        StateLogger::Enable( [&](const std::string &x){ log_lines.push_back(x); } );

        auto bv=box->GetBeads();
        std::vector<std::unique_ptr<CAbstractBead>> original;
        std::vector<std::unique_ptr<CAbstractBead>> testee_state;

        for(unsigned i=0; i<bv.size(); i++){
            original.emplace_back( bv[i]->Clone() );
            testee_state.emplace_back( bv[i]->Clone() );

            check_pbc_drift(bv[i]);
        }

        auto cells=box->GetCNTCells();
        std::vector<std::vector<CAbstractBead*>> original_lbeads;

        original_lbeads.resize( cells.size() );

        CCNTCell::SetCustomRNGProc(
            bead_id_hash_rng__random_symmetric_uniform,
            bead_id_hash_rng__round_hash,
            0
        );

        CSimBox *cbox=const_cast<CSimBox*>(box->GetSimBox());
        
        for(unsigned i=0; i<num_steps; i++){
            for(unsigned j=0; j<original_lbeads.size(); j++){
                auto &cell=cells[j];
                auto beads=cell->GetBeads();
                original_lbeads[j].assign(beads.begin(), beads.end());
            }

            StateLogger::SetPrefix("dut,");
            StateLogger::BeginStep(start_sim_time+i);

            auto res=m_testee->Run(box, true, start_sim_time+i, 1);
            if(res.completed_steps==0){
                return res;
            }

            std::sort(log_lines.begin(), log_lines.end());
            for(const auto &l : log_lines){
                dst<<l<<"\n";
            }
            log_lines.clear();
            
            // Rewind state back to the original
            for(unsigned j=0; j<bv.size(); j++){
                check_pbc_drift(bv[j]);

                testee_state[j]->Assign( bv[j] );
                bv[j]->Assign( original[j].get() );
            }
            for(unsigned j=0; j<cells.size(); j++){
                const auto &beads=original_lbeads[j];
                cells[j]->RemoveAllBeadsFromCell();
                for(unsigned k=0; k<beads.size(); k++){
                    cells[j]->AddBeadtoCell(beads[k]);
                }
            }

            StateLogger::SetPrefix("ref,");
            StateLogger::BeginStep(start_sim_time+i);
            cbox->Evolve(start_sim_time+i);
            cbox->CNTCellCheck();

            std::sort(log_lines.begin(), log_lines.end());
            for(const auto &l : log_lines){
                dst<<l<<"\n";
            }
            log_lines.clear();

            for(unsigned j=0; j<bv.size(); j++){
                auto br=bv[i], bg=testee_state[i].get();

                check_pbc_drift(br);

                double dx[3] = {
                    br->GetXPos() - bg->GetXPos(),
                    br->GetYPos() - bg->GetYPos(),
                    br->GetZPos() - bg->GetZPos()
                };
                double r=sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
                if(r > 1e-6){
                    fprintf(stderr, "dx = %g, %g, %g\n", dx[0], dx[1], dx[2]);
                    exit(1);
                }
                double dv[3] = {
                    br->GetXMom() - bg->GetXMom(),
                    br->GetYMom() - bg->GetYMom(),
                    br->GetZMom() - bg->GetZMom()
                };
                double v=sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
                if(v > 1e-6){
                    fprintf(stderr, "dv = %g, %g, %g\n", dv[0], dv[1], dv[2]);
                    exit(1);
                }
                double df[3] = {
                    br->GetXForce() - bg->GetXForce(),
                    br->GetYForce() - bg->GetYForce(),
                    br->GetZForce() - bg->GetZForce()
                };
                double f=sqrt(df[0]*df[0] + df[1]*df[1] + df[2]*df[2]);
                if(f > 1e-6){
                    fprintf(stderr, "df = %g, %g, %g\n", df[0], df[1], df[2]);
                    exit(1);
                }

                original[j]->Assign(bv[j]);
            }

            StateLogger::Disable();
            std::cout.flush();
        }

        CCNTCell::SetCustomRNGProc(0,0,0);

        return {Supported, {}, num_steps};
    }
protected:
    support_result CanSupport_ExtraConstraints(const ISimBox *box) const override
    {
        return m_testee->CanSupport(box);
    }

private:
    std::shared_ptr<ISimEngine> m_testee;
};

void ISimEngine::WrapGlobalEngineWithRefDiff()
{
    auto e=ISimEngine::GetGlobalEngine();
    if(e){
        ISimEngine::SetGlobalEngine( std::make_shared<SimEngineRefDiff>(e) );
    }

}
