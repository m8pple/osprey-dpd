#ifndef bond_info_hpp
#define bond_info_hpp

#include <cstdint>
#include <vector>
#include <memory>
#include <cmath>
#include <unordered_map>

#include "Polymer.h"
#include "Bond.h"
#include "BondPair.h"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

template<class TCalc=float>
struct BondInfo
{
    struct Bond
    {
        uint32_t head_bead_id;
        uint32_t tail_bead_id;
        TCalc springConst;
        TCalc length;
    };

    struct BondPair
    {
        uint16_t first_bond_off;
        uint16_t second_bond_off;
        TCalc Modulus;
        TCalc m_CosPhi0;
        TCalc m_SinPhi0;
    };
    union BondOrBondPair
    {
        Bond bond;
        BondPair bond_pair;
    };

    struct Polymer
    {
        uint16_t num_bonds;
        uint16_t num_bond_pairs;
        uint32_t colour;
        BondOrBondPair data[1]; // Actually a flexible length array
    };

    struct polymer_free
    {
        void operator()(Polymer *p)
        { free(p); }
    };
    using polymer_ptr = std::unique_ptr<Polymer,polymer_free>;

    std::vector<polymer_ptr> polymers;
    unsigned max_bonds=0;
    std::vector<TCalc> working_space;

    /* For safe parallel updating we need to ensure that each bead is contained in exactly
        one polymer. Each polymer has a colour, and the first polymer to touch a bead gives
        it it's own colour. If a polymer contains a bead with a different colour it can't
        be updated in parallel. */
    std::vector<uint32_t> bead_colours;
    // If true, then all polymers reference distinct beads
    bool polymers_use_disjoint_beads; 

    polymer_ptr import_polymer(uint32_t colour, const CPolymer &src_polymer, std::vector<uint32_t> &working)
    {
        auto &bonds=src_polymer.GetBonds();
        auto &bond_pairs=src_polymer.GetBondPairs();
        unsigned data_count=bonds.size() + bond_pairs.size();
        if(data_count==0){
            return {};
        }

        auto backing=(Polymer*)malloc(sizeof(Polymer) + data_count*sizeof(BondOrBondPair));
        polymer_ptr res(backing);
        res->num_bonds=bonds.size();
        res->num_bond_pairs=bond_pairs.size();
        res->colour=colour;

        working.clear();
        std::unordered_map<const CBond*,unsigned> bond_to_index;
        for(unsigned i=0; i<bonds.size(); i++){
            auto &dst_bond=res->data[i].bond;
            CBond *src_bond=bonds[i];

            dst_bond.head_bead_id = src_bond->GetHead()->GetId()-1;
            working.push_back(dst_bond.head_bead_id);
            
            dst_bond.tail_bead_id = src_bond->GetTail()->GetId()-1;
            working.push_back(dst_bond.tail_bead_id);
            
            dst_bond.length = src_bond->GetUnStrLength();
            dst_bond.springConst = src_bond->GetSprConst();
            bond_to_index[src_bond] = i;
        }


        for(unsigned i=0; i<bond_pairs.size(); i++){
            auto &dst_bond_pair=res->data[bonds.size()+i].bond_pair;
            CBondPair *src_bond_pair=bond_pairs[i];

            dst_bond_pair.first_bond_off=bond_to_index.at(src_bond_pair->GetFirst());
            dst_bond_pair.second_bond_off=bond_to_index.at(src_bond_pair->GetSecond());
            dst_bond_pair.m_CosPhi0=src_bond_pair->GetCosPhi0();
            dst_bond_pair.m_SinPhi0=src_bond_pair->GetSinPhi0();
            dst_bond_pair.Modulus=src_bond_pair->GetModulus();
        }

        // Each bead id appears at least twice.
        // TODO: Is cost of removing redundancy worth it?
        std::sort(working.begin(), working.end());
        working.erase( std::unique(working.begin(), working.end()), working.end() );

        unsigned existing_bead_colour_min = UINT32_MAX;
        for(auto bead_id : working){
            auto existing_bead_colour = bead_colours[bead_id];
            existing_bead_colour_min = std::min(existing_bead_colour_min, existing_bead_colour);
            bead_colours[bead_id] = colour;
        }

        if(existing_bead_colour_min != UINT32_MAX){
            polymers_use_disjoint_beads = false;
            assert(colour < existing_bead_colour_min);
            colour = existing_bead_colour_min;
            res->colour = colour;
        
            for(auto bead_id : working){
                bead_colours[bead_id] = colour;
            }
        }

        return res;
    }

    ISimEngineCapabilities::support_result import_all(ISimBox *box)
    {
        const auto &src_polymers=box->GetPolymers();

        unsigned nbeads=box->GetBeadTotal();
        bead_colours.assign(nbeads, UINT32_MAX);
        polymers_use_disjoint_beads=true;

        polymers.clear();
        max_bonds=0;

        std::vector<uint32_t> working;
        for(unsigned i=0; i<src_polymers.size(); i++){
            const CPolymer *src_polymer=src_polymers[i];
            if(src_polymer->GetSize() <= 1){
                continue;
            }

            auto poly=import_polymer(polymers.size(), *src_polymer, working);
            if(poly){
                max_bonds=std::max<unsigned>(max_bonds, poly->num_bonds);
                polymers.push_back(std::move(poly));
            }
        }

        working_space.resize(max_bonds*4);

        return {ISimEngineCapabilities::Supported};
    }

    template<class TBeadSource>
    void update_polymers_seq(TBeadSource &bead_source, TCalc dims_float[4])
    {
        TCalc local_dims_float[4] = {dims_float[0], dims_float[1], dims_float[2], 0};
        TCalc half_dims_float[4] = {dims_float[0]/2, dims_float[1]/2, dims_float[2]/2, 0};
        for(const auto &p : polymers){
            update_polymer(dims_float, half_dims_float, *p, &working_space[0], bead_source);
        }
    }

    template<class TBeadSource>
    void update_polymers_tbb(TBeadSource &bead_source, TCalc dims_float[4])
    {
        TCalc local_dims_float[4] = {dims_float[0], dims_float[1], dims_float[2], 0};
        TCalc half_dims_float[4] = {dims_float[0]/2, dims_float[1]/2, dims_float[2]/2, 0};
        
        using range_t = tbb::blocked_range<unsigned>;
        tbb::parallel_for( range_t(0, polymers.size()), [&](const range_t &rr){
            TCalc *local_working=0;
            std::vector<TCalc> local_working_space;
            if(working_space.size() <= 1024){
                // Allocate up to 4KB on the current stack. We should be safe for alloca here as:
                // - We don't expect this function to be inlined into a loop, as it is launched as a task by TBB
                // - update_polymer doesn't do anything recursive or use a lot of stack (kind of the point of passing in working space)
                local_working = (TCalc*)alloca(sizeof(TCalc) * working_space.size());
            }else{
                local_working_space.resize(working_space.size());
                local_working = &local_working_space[0];
            }
            
            for(unsigned i=rr.begin(); i<rr.end(); i++){
                update_polymer(dims_float, half_dims_float, *polymers[i], local_working, bead_source);
            }
        });
    }

    template<class TBeadSource>
    void update_polymers(TBeadSource &bead_source, TCalc dims_float[4])
    {
        if(polymers.empty()){
            return;
        }

        if(polymers_use_disjoint_beads){
            update_polymers_tbb(bead_source, dims_float);
        }else{
            update_polymers_seq(bead_source, dims_float);
        }
    }

    /*
    bead_source must map a bead id to a bead reference
    bead reference must have .pos and .force 
    */
    template<class TBeadSource>
    void update_polymer(
        TCalc dims_float[4], TCalc half_dims_float[4],
        const Polymer &polymer, TCalc *working, TBeadSource &bead_source
    )
    {
        for(unsigned i=0; i<polymer.num_bonds; i++){
            const Bond &bond=polymer.data[i].bond;
            auto &head = bead_source(bond.head_bead_id);
            auto &tail = bead_source(bond.tail_bead_id);

            TCalc dx[4];
            for(int d=0; d<3; d++){
                dx[d] = head.pos[d] - tail.pos[d];
               if(dx[d] < -half_dims_float[d]){
                    dx[d] += dims_float[d];
                }else if(dx[d] > half_dims_float[d]){
                    dx[d] -= dims_float[d];
                }
                working[ 4*i + d ] = dx[d];
            }

            TCalc r=sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
            assert(std::abs(r)<4);
            working[4*i+3] = r;

            TCalc fScale = bond.springConst * (bond.length - r) / r;

            for(int d=0; d<3; d++){
                head.force[d] += dx[d] * fScale;
                tail.force[d] -= dx[d] * fScale;
            
                assert(!std::isnan(head.force[d]));
                assert(!std::isnan(tail.force[d]));
            }
            
        }

        for(unsigned i=0; i<polymer.num_bond_pairs; i++){
            const BondPair &bond_pair=polymer.data[polymer.num_bonds + i].bond_pair;
            const TCalc *first = working + bond_pair.first_bond_off*4;
            const TCalc *second = working + bond_pair.second_bond_off*4;
            auto &first_head_bead = bead_source( polymer.data[bond_pair.first_bond_off].bond.head_bead_id );
            auto &first_tail_bead = bead_source( polymer.data[bond_pair.first_bond_off].bond.tail_bead_id );
            auto &second_head_bead = bead_source( polymer.data[bond_pair.second_bond_off].bond.head_bead_id );

            TCalc FirstLength=first[3];
            assert(std::abs(FirstLength) < 4);
            TCalc SecondLength=second[3];
            assert(std::abs(SecondLength) < 4);

            TCalc magProduct = FirstLength * SecondLength;

            if(magProduct > TCalc(0.0001)) {
                TCalc b1MagSq		= FirstLength*FirstLength;
                TCalc b2MagSq		= SecondLength*SecondLength;
                TCalc b1Dotb2		= first[0]*second[0] + first[1]*second[1] + first[2]*second[2];
                TCalc b1b2Overb1Sq	= b1Dotb2/b1MagSq;
                TCalc b1b2Overb2Sq	= b1Dotb2/b2MagSq;
                TCalc cosPhiSq		= b1b2Overb1Sq*b1b2Overb2Sq;

        		TCalc forceMag = bond_pair.Modulus/magProduct;

                // Check that the bond angle is not exactly 90 deg but allow the cosine to be < 0

                if(bond_pair.m_SinPhi0!=0 && std::fabs(b1Dotb2) > TCalc(0.000001))
                {
                    TCalc Prefactor = sqrt(TCalc(1)/cosPhiSq - TCalc(1));

                    // phi0 > 0 -> sin(phi) > 0,   for phi in [0,pi/2)
                    // phi0 ==0 -> sin(phi)==0 && cos(phi)==1
                    /*if(bond_pair.phi0 > 0.0)
                    { forceMag = m_Modulus*(m_CosPhi0 - m_SinPhi0/Prefactor)/magProduct; }
                    else
                    { forceMag = m_Modulus/magProduct; }*/

                    forceMag *= (bond_pair.m_CosPhi0 - bond_pair.m_SinPhi0 / Prefactor);
                    assert(!std::isnan(forceMag));
                }

        		TCalc BeadXForce[3], BeadYForce[3], BeadZForce[3];

                BeadXForce[0] = forceMag*(b1b2Overb1Sq*first[0] - second[0]);
                BeadYForce[0] = forceMag*(b1b2Overb1Sq*first[1] - second[1]);
                BeadZForce[0] = forceMag*(b1b2Overb1Sq*first[2] - second[2]);
                assert(!std::isnan(BeadXForce[0]));

                first_tail_bead.force[0] += BeadXForce[0];
                first_tail_bead.force[1] += BeadYForce[0];
                first_tail_bead.force[2] += BeadZForce[0];

                // Store the force on the third bead for use in calculating stresses

                BeadXForce[2] = forceMag*(first[0] - b1b2Overb2Sq*second[0]);
                BeadYForce[2] = forceMag*(first[1] - b1b2Overb2Sq*second[1]);
                BeadZForce[2] = forceMag*(first[2] - b1b2Overb2Sq*second[2]);

                second_head_bead.force[0] += BeadXForce[2];
                second_head_bead.force[1] += BeadYForce[2];
                second_head_bead.force[2] += BeadZForce[2];

                // Restoring force on the middle bead is chosen so as to make the 
                // total force generated by the bond stiffness vanish, hence the
                // minus sign. Note that the first bond's head bead is the same 
                // as the second's tail.

                BeadXForce[1] = -BeadXForce[0] - BeadXForce[2];
                BeadYForce[1] = -BeadYForce[0] - BeadYForce[2];
                BeadZForce[1] = -BeadZForce[0] - BeadZForce[2];

                first_head_bead.force[0] += BeadXForce[1];
                first_head_bead.force[1] += BeadYForce[1];
                first_head_bead.force[2] += BeadZForce[1];

                for(int d=0; d<3; d++){
                    assert(!std::isnan(BeadXForce[d]));
                    assert(!std::isnan(BeadYForce[d]));
                    assert(!std::isnan(BeadZForce[d]));
                }

            }
        }
    }
};

#endif
