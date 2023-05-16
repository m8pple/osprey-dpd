#ifndef SimEngineTBB_hpp
#define SimEngineTBB_hpp

#include "SimEngineSeq.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <cstring>


struct EnginePolicyTBB
    : EnginePolicyConcept
{
    static constexpr NhoodType NHOOD_TYPE = (NhoodType)( NhoodType_Size_Full | NHoodType_ForceShare_None );
    static constexpr RNGPolicy RNG_POLICY = RNGPolicy_Hash_BeadTagXorShiftAdd;
};

template<class TPolicy=EnginePolicyTBB>
struct SimEngineTBB
    : public SimEngineSeq<TPolicy>
{
public:
    using base_t = SimEngineSeq<TPolicy>;
    using Cell = typename base_t::Cell;
    using rng_t = typename base_t::rng_t;
    using Bead = typename base_t::Bead;

    using base_t::CALC_DIM;
    using base_t::cells;
    using base_t::edge_adjustments;        
    using base_t::MAX_BEADS_PER_CELL;
    using base_t::dims_float;
    using base_t::pos_to_cell_index;
    using base_t::flush_unpbc_wrap;
    using base_t::bead_locations;
    using base_t::validate_cells;
    using base_t::find;
    using base_t::bonds;
    using base_t::round_id;


    std::string Name() const override
    {
        std::string res="SimEngineTBB_";
        res += rng_t::Name();
        if(base_t::USE_MORTON){
            res += "Morton";
        }
        return res;
    }

    bool IsParallel() const override
    { return true; }

    ISimEngine::run_result Run(ISimBox *box, bool modified, unsigned num_steps) override
    {
        auto err=base_t::import_all(box);
        if(err.status!=ISimEngine::Supported){
            return err;
        }

        auto bead_source=[&](uint32_t bead_id) -> Bead &{
            return find(bead_id);
        };

        #ifndef NDEBUG
        validate_cells();
        #endif

        for(unsigned i=0; i<num_steps; i++){
            update_mom_and_move();

            #ifndef NDEBUG
            validate_cells();
            #endif

            update_forces();

            #ifndef NDEBUG
            validate_cells();
            #endif

            bonds.update_polymers_tbb(
                bead_source,
                dims_float
            );

            #ifndef NDEBUG
            validate_cells();
            #endif

            round_id += 1;
        }

        #ifndef NDEBUG
        validate_cells();
        #endif

        base_t::export_all(box);
    
        return {ISimEngine::Supported, {}, num_steps};
    }
protected:

    std::mutex mutex_barrier;

    template<class TF, class TP>
    void for_each_cell(TF &&f, uint64_t uniquifier_base, TP &&prefetch)
    {
        using range_t = tbb::blocked_range<unsigned>;

        {
            std::unique_lock<std::mutex> lk(mutex_barrier);
        }

        if(1){
            unsigned n=base_t::cells.size();
            tbb::parallel_for( range_t(0,n), [&]( const range_t &rr ){
                rng_t rng(base_t::rng_stddev, base_t::global_seed, base_t::round_id, uniquifier_base + rr.begin());

                for(unsigned i=rr.begin(); i<rr.end(); i++){
                    if(i+1 < rr.end()){
                        prefetch( base_t::cells[i+1] );
                    }
                    f(rng,  base_t::cells[i] );
                }
            });
        }else{

            unsigned n=base_t::cells.size();
            rng_t rng(base_t::rng_stddev, base_t::global_seed, base_t::round_id, 0);

            for(unsigned i=0; i<n; i++){
                if(i+1 < n){
                    prefetch( base_t::cells[i+1] );
                }
                f(rng,  base_t::cells[i] );
            }
        }

        {
            std::unique_lock<std::mutex> lk(mutex_barrier);
        }
    }


    template<bool AnyExternal>
    void update_forces(rng_t &rng, Cell &home_cell)
    {
        home_cell.to_move = home_cell.count;

        if(home_cell.count==0){
            return;
        }

        for(unsigned i0=0; i0<home_cell.count-1; i0++){
            for(unsigned i1=i0+1; i1<home_cell.count; i1++){
                base_t::calc_force(rng, true, home_cell.local[i0], home_cell.local[i1], home_cell.local[i1].pos );
            }
        }

        int nhood_size=base_t::GetNhoodSize(home_cell);
        assert(nhood_size==26);
        for(unsigned nhood_index=0; nhood_index<base_t::GetNhoodSize(home_cell); nhood_index++){
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
                    other_delta[d] = edge_adjustments[d][nlink.offsets[d]];
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

                    base_t::calc_force(rng, false, home_bead, other_bead, AnyExternal ? other_x : other_bead.pos);
                }
            }
        }
    }

    virtual __attribute__((noinline)) void update_forces() override
    {
        // TODO : branching based on any external may introduce data-dependent
        // control and bloat instruction size. Is minor saving worth it?
        for_each_cell(
            [&](rng_t &rng, Cell &cell){
                if(cell.any_external){
                    update_forces<true>(rng, cell);
                }else{
                    update_forces<false>(rng, cell);
                }
            },
            0,
            [](Cell &cell){}
        );
    }

    /*
    Locking strategy is a little complicated. We have
        [0, to_move) - Beads that will be moved (are owned) by this thread
        [to_move,count) - New arrivals from other threads.

        The beads in [0,to_move) are thread private, and can be modified directly.
        So we iterate over the cells as normal.
        
        When a bead leaves, we can insert it into the new location by:
        - Taking a lock on dest cell
            - That gives us read-write access to count and beads in [to_move,MAX_BEADS_PER_CELL)
        - Copy bead into dest cell count
        - Increase dest cell count
        - Unlock cell
        At no point do we interfere with local cells being processed
        Other threads neveer modify to_move, even under lock.

        We end up with a gap in [0,to_move), which needs to be filled. We cannot
        directly change to_move without locking, so we deferr this to the end.
        Then under lock we move beads from [to_move,count) into the gaps.

        Overally this means that:
        - If no beads leave a cell, no locks are needed
        - Each bead movement requires 1 lock per bead (dest cell)
        - Each cell with a bead that leaves requires 1 lock per cell (source cell)

        Assuming beads leaves cells relatively infrequently, most cells don't need to take a lock.

        TODO: this can be optimised if we cluster cells and give each cluster a colour. When
        moving to a cluster with the same colour no locking is needed.
    */
    template<bool AnyFrozen>
    void update_mom_and_move_locked(rng_t &rng, Cell &cell)
    {
        assert(AnyFrozen ? base_t::any_frozen_beads : true);

        #ifndef NDEBUG
        {
            std::unique_lock<std::mutex> lk(cell.mutex);
            for(unsigned i=0; i<cell.count; i++){
                assert(bead_locations[cell.local[i].bead_id] == cell.cell_index*MAX_BEADS_PER_CELL+i);
            }
        }
        #endif

        unsigned i=0;
        unsigned new_to_move=cell.to_move;
        unsigned gap=0;
        
        while(i<new_to_move){
            Bead &bead=cell.local[i];

            if(AnyFrozen && bead.frozen){
                for(int d=0; d<CALC_DIM; d++){
                    bead.force[d] = 0;
                }
                continue;
            }

            bool moved=false;

            float orig_force[4];

            for(int d=0; d<CALC_DIM; d++){
                orig_force[d] = bead.force[d];

                // Mom from previous step (loop skewed)
                bead.mom[d] += base_t::halfdt * bead.force[d];
                bead.pos[d] += base_t::dt * bead.mom[d] + base_t::halfdt2 * bead.force[d];
		        // We know that lambda=0.5
                bead.mom[d] += base_t::halfdt * bead.force[d];
        		bead.force[d] = 0;

                float new_origin_d = std::floor( bead.pos[d] );
                moved |= new_origin_d != cell.origin[d];

                if(rng_t::USES_BEAD_ROUND_TAG){
                    assert(0); //Needs to be thread safe
                    bead.round_tag = rng.MakeBeadTag(bead.bead_id);
                }
            }

            if(!moved){
                ++i;
                continue;
            }


            float pre_wrap_pos[4];
            int new_origin[4];
            for(int d=0; d<CALC_DIM; d++){
                pre_wrap_pos[d] = bead.pos[d];
                if(bead.pos[d] < 0){
                    bead.pos[d] += dims_float[d];
                    bead.unpbcWrap[d] -= 1;
                    if( bead.unpbcWrap[d] < -100 ){
                        // This does not need to be locked, as each bead uniquely maps to a thread
                        flush_unpbc_wrap(bead);
                    }
                }
                // The above could results in bead.pos[d] == dims_float[d] due to
                // rounding, so this is not an else if
                if(bead.pos[d] >= dims_float[d]){
                    bead.pos[d] -= dims_float[d];
                    bead.unpbcWrap[d] += 1;
                    if( bead.unpbcWrap[d] < -100 ){
                        // Doesn't need to be locked
                        flush_unpbc_wrap(bead);
                    }
                }
                new_origin[d] = std::floor( bead.pos[d] );
            }

            unsigned new_index = pos_to_cell_index(new_origin);
            Cell &new_cell = cells[new_index];

            // Move the bead into new cell under lock
            {
                ///////////////////////////////////////////////////////////
                // Enter lock
                std::unique_lock<std::mutex> lk(new_cell.mutex);
                
                // updating location actually needs to be under the lock, as a given bead might get
                // modified twice in a round:
                // - Moved into a new cell (this critical section)
                // - Moved into gap in existing cell (later critical section)
                // We need to worry about WAW hazard, so do it under lock.
                bead_locations[bead.bead_id] = new_index*MAX_BEADS_PER_CELL + new_cell.count;
                
                memcpy( new_cell.local+new_cell.count, &bead, sizeof(Bead) );
                new_cell.count++;

                #ifndef NDEBUG
                for(unsigned j=new_cell.to_move; j<new_cell.count; j++){
                    assert( bead_locations[new_cell.local[j].bead_id] == new_cell.cell_index * MAX_BEADS_PER_CELL + j  );
                }
                #endif

                #ifndef NDEBUG
                // Poison the bead in the gap
                memset( &bead, -1, sizeof(Bead) );
                #endif

                // Exit lock
                ///////////////////////////////////////////////////////////
            }

            // We have created a gap at index i
            // Processed: [0,i)
            // Gap:       [i,i+1)
            // Todo:      [i+1,new_to_move)
            gap += 1;

            assert(bead.bead_id==UINT32_MAX); // Should be poisoned
            assert(&bead == cell.local+i);

            if(i+1 == new_to_move){
                // Processed: [0,i)                             This is all that is left
                // Gap:       [i,i+1==new_to_move)
                // Todo:      [i+1==new_to_move,new_to_move).   So this is an empty set.
                // The gap is the last thing in the local array
                new_to_move -= 1;
                // We now exit loop and patch gaps under local lock
                break;
            }else{
                assert(i < new_to_move);
                // Processed: [0,i)       
                // Gap:       [i,i+1)
                // Todo:      [i+1,new_to_move).   i+1<new_to_move, so the set is not empty
                memcpy( cell.local+i, cell.local+new_to_move-1, sizeof(Bead) );
                bead_locations[cell.local[i].bead_id] = cell.cell_index*MAX_BEADS_PER_CELL+i;

                #ifndef NDEBUG
                // Poison the bead at the end
                memset( cell.local+new_to_move-1, -1, sizeof(Bead) );
                #endif

                // i remains the same. We process the bead just copied on the next iteration.
                new_to_move -= 1;
            }  
        }

        #ifndef NDEBUG
        {
            std::unique_lock<std::mutex> lk(cell.mutex);
            for(unsigned i=0; i<new_to_move; i++){
                assert(bead_locations[cell.local[i].bead_id] == cell.cell_index*MAX_BEADS_PER_CELL+i);
            }
        }
        #endif

        if(gap){
            /////////////////////////////////////////////////////////
            // Enter lock
            std::unique_lock<std::mutex> lk(cell.mutex);

            assert(new_to_move + gap == cell.to_move);
            assert(cell.to_move <= cell.count);

            // Valid cells: [0,new_to_move)
            // Gap:         [new_to_move,to_move)  ==  [new_to_move,new_to_move+gap)
            // Incoming:    [to_move,count)        ==  [new_to_move+gap,count)

            #ifndef NDEBUG
            for(unsigned j=new_to_move; j<cell.to_move; j++){
                assert(cell.local[j].bead_id == UINT32_MAX); // Should all be poisoned
            }
            #endif

            unsigned new_arrivals=cell.count - cell.to_move;

            unsigned moved=0;

            if(new_arrivals==0){
                // Nobody has tried to insert a new bead, so we are already packed
            }else if(gap >= new_arrivals){
                // The gap is the same size or bigger than the new arrivals
                memcpy( cell.local+new_to_move, cell.local+cell.to_move, new_arrivals*sizeof(Bead) );
                moved=new_arrivals;
            }else{
                // The number of new arrivals is larger than the gap
                memcpy( cell.local+new_to_move, cell.local+cell.count-gap, gap*sizeof(Bead) );
                moved=gap;
            }

            // This must be done under lock due to possible WAW hazard
            for(unsigned j=new_to_move; j<new_to_move+moved; j++){
                bead_locations[cell.local[j].bead_id] = cell.cell_index*MAX_BEADS_PER_CELL+j;
            }

            cell.to_move = new_to_move+new_arrivals;
            cell.count = new_to_move+new_arrivals;

            lk.unlock();
            // Exit lock
            /////////////////////////////////////////////////////
        }

        #ifndef NDEBUG
        {
            std::unique_lock<std::mutex> lk(cell.mutex);
            for(unsigned i=0; i<cell.count; i++){
                assert(bead_locations[cell.local[i].bead_id] == cell.cell_index*MAX_BEADS_PER_CELL+i);
            }
        }
        #endif
    }

     __attribute__((noinline)) void update_mom_and_move() override
    {
        if(base_t::any_frozen_beads){
            for_each_cell([&](rng_t &rng, Cell &cell){
                if(cell.any_frozen){
                    update_mom_and_move_locked<true>(rng, cell);
                }else{
                    update_mom_and_move_locked<false>(rng, cell);
                }
            },
            1,
            [&](Cell &next){
                __builtin_prefetch(&next.count);
                __builtin_prefetch(&next.local[0]);
            });
        }else{
            for_each_cell([&](rng_t &rng, Cell &cell){
                update_mom_and_move_locked<false>(rng, cell);
            },
            1,
            [&](Cell &next){
                __builtin_prefetch(&next.count);
                __builtin_prefetch(&next.local[0]);
            });
        }
    }
};

#endif
