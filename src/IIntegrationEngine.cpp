#include "IIntegrationEngine.h"

#include "CNTCell.h"
#include "ISimBox.h"
#include "aeActiveSimBox.h"

#include "BeadIdHashRNG.h"
#include "StateLogger.h"

#include <cfloat>
#include <climits>
#include <algorithm>
#include <utility>

IIntegrationEngineCapabilities::support_result IIntegrationEngineCapabilities::StandardSupportAssumptions(const ISimBox *box)
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

IIntegrationEngineCapabilities::support_result IIntegrationEngineCapabilities::CanSupport(const ISimBox *box) const
{
    auto tmp=StandardSupportAssumptions(box);
    if(tmp.status!=Supported){
        return tmp;
    }
    return CanSupport_ExtraConstraints(box);
}



class SimEngineRefDiff
    : public IIntegrationEngine
{
public:
    SimEngineRefDiff(std::shared_ptr<IIntegrationEngine> testee)
        : m_testee(testee)
    {}

    std::string Name() const override
    {
        return "SimEngineRefDiff_"+m_testee->Name();
    }

    bool IsParallel() const override
    { return m_testee->IsParallel(); }

    struct stats_tracker
    {
        stats_tracker(const std::string &_name)
            : name(_name)
        {}

        std::string name;
        double n=0;
        long double sum = 0, sum_sqr = 0, sum_cube = 0, sum_quart=0;
        double max = -DBL_MAX, min=DBL_MAX;
        double max_abs = -DBL_MAX;
        unsigned max_bead_id = 0;

        void add(unsigned bead_id, double x)
        {
            n += 1;
            sum += x;
            sum_sqr += x*x;
            sum_cube += x*x*x;
            sum_quart += x*x*x*x;
            max = std::max(max, x);
            min = std::min(min, x);
            if(std::abs(x) > max_abs){
                max_bead_id=bead_id;
                max_abs=std::abs(x);
            }
        }

        void add_err_mag(unsigned bead_id, const double *ref, const double *dut)
        {
            double dr2=0;
            for(unsigned i=0; i<3; i++){
                double e=ref[i] - dut[i];
                dr2 += e*e;
            }
            add(bead_id, sqrt(dr2));
        }

        void add_rel_err_mag(unsigned bead_id, const double *ref, const double *dut)
        {
            double dr2=0;
            double dref=0;
            for(unsigned i=0; i<3; i++){
                double e=ref[i] - dut[i];
                dr2 += e*e;
                dref += ref[i]*ref[i];
            }
            add(bead_id, sqrt(dr2 / dref));
        }

        void add_err_mag_wrapped(unsigned bead_id, const double *ref, const double *dut, const double *box)
        {
            double dr2=0;
            for(unsigned i=0; i<3; i++){
                double e=ref[i] - dut[i];
                if(e < -box[i]/2){
                    e += box[i];
                }
                if(e > box[i]/2){
                    e -= box[i];
                }
                dr2 += e*e;
            }
            add(bead_id, sqrt(dr2));
        }

        void print_stats(std::ostream &dst) const
        {
            auto mean=sum/n;
            dst<<name<<", "<<min<<","<<mean<<","<<max<<", ";
            auto var=sum_sqr/n - mean*mean;
            if(var <= 0){
                dst<<"nan,nan,nan, -1";
            }else{
                auto std=sqrt(var);
                auto skew=(sum_cube/n-3*mean*std*std-mean*mean*mean)/(std*std*std);
                auto kurt=(sum_quart/n-4*mean*sum_cube/n+6*mean*mean*std*std+3*mean*mean*mean*mean)/(std*std*std*std);
                dst<<std<<","<<skew<<","<<kurt<<", "<<max_bead_id;
            }
            dst<<"\n";
        }
    };

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

        std::vector<std::pair<std::string,bool> > log_lines;

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
        
        unsigned i=0;
        while(i<num_steps){
            unsigned todo=std::min<unsigned>(num_steps-i, 1);

            for(unsigned j=0; j<original_lbeads.size(); j++){
                auto &cell=cells[j];
                auto beads=cell->GetBeads();
                original_lbeads[j].assign(beads.begin(), beads.end());
            }

            StateLogger::Enable( [&](const std::string &x){ log_lines.push_back({x,false}); } );

            auto res=m_testee->Run(box, true, start_sim_time+i, todo);
            if(res.completed_steps<todo){
                return res;
            }

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

            StateLogger::Enable( [&](const std::string &x){ log_lines.push_back({x,true}); } );

            for(unsigned ioff=0; ioff<todo; ioff++){
                StateLogger::BeginStep(start_sim_time+i+ioff);
                cbox->Evolve(start_sim_time+i+ioff);
                cbox->CNTCellCheck();
            }

            // Sigh. We need sort order that only depends on first 5 fields
            auto sort_key=[](const std::pair<std::string,bool> &a, const std::pair<std::string,bool> &b)
            {
                  int pa=a.first.find(',',0), pb=b.first.find(',',0);
                  for(int i=0; i<4; i++){
                    pa=a.first.find(',',pa+1);
                    pb=a.first.find(',',pb+1);
                  }
                  if(pa<pb){
                    return true;
                  }
                  if(pa>pb){
                    return false;
                  }
                  auto sa=a.first.substr(0,pa);
                  auto sb=b.first.substr(0,pb);
                  if(sa < sb){
                    return true;
                  }
                  if(sa > sb){
                    return false;
                  }
                  return a.second < b.second;
            };

            std::sort(log_lines.begin(), log_lines.end(), sort_key);
            for(const auto &l : log_lines){
                if(l.second){
                    dst<<"ref,"<<l.first<<"\n";
                }else{
                    dst<<"dut,"<<l.first<<"\n";
                }
            }
            log_lines.clear();
            dst.flush();

            double dims[3] = {box->GetSimBoxXLength(), box->GetSimBoxYLength(), box->GetSimBoxZLength()};
            stats_tracker stats_pos_err_mag{"pos_err_mag"};
            stats_tracker stats_mom_rel_err_mag{"mom_rel_err_mag"};
            stats_tracker stats_force_rel_err_mag{"force_rel_err_mag"};
            for(unsigned j=0; j<bv.size(); j++){
                auto bead_id=bv[j]->GetId()-1;
                auto br=bv[j], bg=testee_state[j].get();

                check_pbc_drift(br);

                stats_pos_err_mag.add_err_mag_wrapped( bead_id, br->GetPos(), bg->GetPos(), dims );
                stats_mom_rel_err_mag.add_rel_err_mag( bead_id, br->GetMom(), bg->GetMom() );
                stats_force_rel_err_mag.add_rel_err_mag( bead_id, br->GetForce(), bg->GetForce() );

                original[j]->Assign(bv[j]);
            }

            auto &sdest=std::cerr;
            
            sdest<<m_testee->Name()<<","<<start_sim_time<<","<<todo<<",";
            stats_pos_err_mag.print_stats(std::cerr);
            
            sdest<<m_testee->Name()<<","<<start_sim_time<<","<<todo<<",";
            stats_mom_rel_err_mag.print_stats(std::cerr);
            
            sdest<<m_testee->Name()<<","<<start_sim_time<<","<<todo<<",";
            stats_force_rel_err_mag.print_stats(std::cerr);

            i+=todo;

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
    std::shared_ptr<IIntegrationEngine> m_testee;
};

void IIntegrationEngine::WrapGlobalEngineWithRefDiff()
{
    auto e=IIntegrationEngine::GetGlobalEngine();
    if(e){
        IIntegrationEngine::SetGlobalEngine( std::make_shared<SimEngineRefDiff>(e) );
    }

}
