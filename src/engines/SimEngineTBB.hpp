#ifndef SimEngineTBB_hpp
#define SimEngineTBB_hpp

#include "SimEngineSeq.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

struct SimEngineTBB
    : public SimEngineSeq
{
public:
    std::string Name() const override
    {
        return "SimEngineTBB";
    }

    bool IsParallel() const override
    { return true; }



protected:

    template<class TF, class TP>
    void for_each_cell(TF &&f, TP &&prefetch)
    {
        using range_t = tbb::blocked_range<unsigned>;

        unsigned n=cells.size();
        tbb::parallel_for( range_t(n), [&]( const range_t &rr ){
            for(unsigned i : rr){
                if(i+1 < n){
                    prefetch( cells[i] );
                }
                f( cells[i] );
            }
        }
    }


    float RandUnifScaled()
    {
        uint32_t u32 = rng_state>>32;
        rng_state=6364136223846793005ull * rng_state + 1;
        
        int32_t i31;
    	memcpy(&i31, &u32, 4); // Avoid undefined behaviour. Gets number in range [-2^31,2^31)
	    return i31 * rng_scale; // rng_scale = ( CCNTCell::m_invrootdt * 2^-32  )
    }


    template<bool AnyExternal>
    void update_forces(Cell &home_cell)
    {
        home_cell.to_move = home_cell.count;

        if(home_cell.count==0){
            return;
        }

        for(unsigned i0=0; i0<home_cell.count-1; i0++){
            for(unsigned i1=i0+1; i1<home_cell.count; i1++){
                calc_force( home_cell.local[i0], home_cell.local[i1], home_cell.local[i1].pos );
            }
        }

        for(unsigned nhood_index=0; nhood_index<13; nhood_index++){
            Cell &prefetch_cell = cells[home_cell.nhood[nhood_index+1].index];
            // Get the header of the cell
            __builtin_prefetch(&prefetch_cell.count);
            // First bead - note that it is offset from the cell header by quite a distance due to nhood links
            __builtin_prefetch(&prefetch_cell.local[0]);
            
            auto nlink=home_cell.nhood[nhood_index];
            Cell &other_cell=cells[nlink.index];

            float other_delta[4];
            if(AnyExternal){
                for(int d=0; d<3; d++){
                    other_delta[d] = edge_adjustments[d][nlink.edge_tags[d]];
                }
                other_delta[3]=0;
            }

            for(unsigned other_i=0; other_i<other_cell.count; other_i++){
                Bead &other_bead = other_cell.local[other_i];
                float other_x[4];
                if(AnyExternal){
                    for(int d=0; d<CALC_DIM; d++){
                        other_x[d] = other_bead.pos[d] + other_delta[d];
                    }
                }

                for(unsigned home_i=0; home_i<home_cell.count; home_i++){
                    Bead &home_bead = home_cell.local[home_i];

                    calc_forces(home_bead, other_bead, AnyExternal ? other_x : other_bead.pos);
                }
            }
        }
    }

    virtual __attribute__((noinline)) void update_forces()
    {
        // TODO : branching based on any external may introduce data-dependent
        // control and bloat instruction size. Is minor saving worth it?
        for_each_cell([&](Cell &cell){
            if(cell.any_external){
                update_forces<true>(cell);
            }else{
                update_forces<false>(cell);
            }
        },[](Cell &cell){

        });
    }
s
};

#endif
