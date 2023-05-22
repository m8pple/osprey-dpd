#ifndef SimEngineSeq_hpp
#define SimEngineSeq_hpp

#include <cstdint>
#include <vector>
#include <atomic>
#include <vector>
#include <mutex>
#include "DebugAssert.hpp"
#include <cstdlib>
#include <set>
#include <cstring>
#include <unordered_set>

#include "AbstractBead.h"
#include "ISimBox.h"
#include "xxBase.h"
#include "CNTCell.h"
#include "Polymer.h"
#include "Bond.h"
#include "BondPair.h"

#include "ISimEngine.h"

#include "morton_codec.hpp"

#include "bond_info.hpp"

#include "RNGPolicy.hpp"

static inline bool require_fail_impl(const char *file, int line, const char *cond)
{
    fprintf(stderr, "%s:%d: requirement failed : %s\n", file, line, cond);
    DEBUG_ASSERT(0);
    exit(1);
}

#define require(cond) \
    if(!(cond)){ require_fail_impl(__FILE__, __LINE__, #cond); }

static inline void declare(bool cond)
{
    if(!cond){
        DEBUG_ASSERT(0);
        __builtin_unreachable();
    }
}

enum NhoodType
{
    NhoodType_Size_Half = 0,
    NhoodType_Size_Full = 1,

    NhoodType_Size_Mask = 1,

    NHoodType_ForceShare_None =   0x10,  // No forces shared across boundary
    NHoodType_ForceShare_Direct = 0x20, // Forces added directly with no locking

    NhoodType_ForceShare_Mask =   0xf0
};

struct EnginePolicyConcept
{
    static const RNGPolicy RNG_POLICY = RNGPolicy_Rng_LCG64;
    static constexpr bool USE_MORTON = false;
    static constexpr NhoodType NHOOD_TYPE = (NhoodType)( NhoodType_Size_Half | NHoodType_ForceShare_Direct );
};


template<class TPolicy = EnginePolicyConcept>
struct SimEngineSeq
    : public ISimEngine
{
public:
    static constexpr bool USE_MORTON = TPolicy::USE_MORTON;
    static constexpr RNGPolicy RNG_POLICY = TPolicy::RNG_POLICY;
    static constexpr NhoodType NHOOD_TYPE = TPolicy::NHOOD_TYPE;

    using rng_t = RNGImpl<RNG_POLICY>;

    std::string Name() const override
    {
        std::string res="SimEngineSeq_";
        res += rng_t::Name();
        if(USE_MORTON){
            res += "Morton";
        }
        return res;
    }

    bool IsParallel() const override
    { return false; }

    support_result CanSupport_ExtraConstraints(const ISimBox *box) const override
    {
        if(box->GetBeadTypeTotal() >= 255){
            return {PermanentProblem, "More than 254 bead types"};
        }
        return {Supported};
    }

    run_result Run(ISimBox *box, bool modified, unsigned /*start_sim_time*/, unsigned num_steps) override
    {
        auto res=import_all(box);
        if(res.status!=Supported){
            return res;
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

            bonds.update_polymers(
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

        export_all(box);

        return {Supported, {}, num_steps};
    }

protected:

    static const int MAX_BEADS_PER_CELL = 32;

    static const int CALC_DIM = 3; // Dimension used for maths. Could be 3 for scalar, or 4 to allow for more optimisation in vector modes

    // We expect PBC and unPBC positions to be within this distance of each other
    // at all times (particularly import)
    static constexpr double UNPBC_DRIFT_TOLERANCE = 1e-6;

    struct Bead
    {
        float pos[4];
        float mom[4];
        float force[4];
        
        uint8_t type : 7;
        uint8_t frozen : 1;
        int8_t unpbcWrap[3]; // Offset of bead position in PBC terms 

        uint32_t bead_id;
        uint32_t round_tag; // Only used if hashing RNG policy is used
        uint32_t _pad_[1];
    };
    static_assert(sizeof(Bead)==4*4*4, "Bead is expected to be 64 bytes.");

    enum EdgeTag : uint8_t
    {
        EdgeTag_ForceShare_None = 0,
        EdgeTag_ForceShare_Direct = 1
    };

    struct Neighbour
    {
        uint32_t index;
        EdgeTag tag;
        std::array<uint8_t,3> offsets; // Tells us what boundary to apply to positions in this cell
    };

    Neighbour make_edge_tag(const int32_t *home_pos, const int32_t *neighbour_pos, const int32_t *dims)
    {
        Neighbour local;
        local.index=pos_to_cell_index(neighbour_pos);
        for(int d=0; d<3; d++){
            int delta=neighbour_pos[d] - home_pos[d];
            if(delta==0 || delta==1 || delta==-1){
                local.offsets[d]=0;
            }else if(delta==1-dims[d]){
                local.offsets[d]=1;
            }else if(delta==dims[d]-1){
                local.offsets[d]=2;
            }else{
                DEBUG_ASSERT(0);
                fprintf(stderr, "Logic violation in make_edge_tag");
                exit(1);
            }
            if( (NHOOD_TYPE & NhoodType_ForceShare_Mask) == NHoodType_ForceShare_None){
                local.tag=EdgeTag_ForceShare_None;
            }else if( (NHOOD_TYPE & NhoodType_ForceShare_Mask) == NHoodType_ForceShare_Direct){
                local.tag=EdgeTag_ForceShare_Direct;
            }else{
                require(0); // Don't know how to choose
            }

        }
        return local;
    }

    // It is deliberately 14 = 13+1, so that we can use the final
    // extra entry to act as prefetch for the next cell
    static constexpr unsigned NHOOD_ARRAY_SIZE = (NHOOD_TYPE & NhoodType_Size_Mask) == NhoodType_Size_Half ? 14 : 27;

    struct Cell
    {
        // Hottest variables for force calculation
        uint32_t count;
        uint32_t to_move; // to_move <= count
        bool any_frozen;
        bool any_external;

        float origin[CALC_DIM];
        uint32_t cell_index;

        Neighbour nhood[NHOOD_ARRAY_SIZE];

        Bead local[MAX_BEADS_PER_CELL];

        std::mutex mutex; // Place at the end. Ideally not hot.
        uint8_t __pad__[64];
    };

    unsigned GetNhoodSize(const Cell &c)
    {
        switch( (NHOOD_TYPE&NhoodType_Size_Mask) ){
        case NhoodType_Size_Half: return 13;
        case NhoodType_Size_Full: return 26;
        default:
            require(0);
            return 0;
        }
    }

    float dims_float[4];
    int32_t dims_int[4];
    float half_dims_float[4];
    float dt;
    float halfdt;
    float halfdt2;

    uint64_t global_seed;
    uint64_t round_id;
    double rng_stddev;
    
    float edge_adjustments[3][3];

    uint32_t num_bead_types;
    std::vector<std::pair<float,float>> coefficients;

    std::vector<Cell> cells;
    std::vector<int32_t> bead_locations; // The current location of each slot

    // If we need to flush the small bead->unpbcWrap counters then they go here.
    // We don't want to flush into the AbstractBead in case some kindof error
    // happens and we abort.
    // We also don't want them in the beads, as it increases cache pressur for no
    // reason.
    // Flushing the pbc requires a bead to wrap round the entire box 127 times,
    // so accessing this array is very rare.
    std::vector<int32_t> unpbc_flushed_fixups; // An array of 3*NumBeads integers

    // Indices of cells in order of visiting. Used for morton orderig.
    std::vector<uint32_t> cell_enum_order;

    BondInfo<float> bonds;

    AbstractBeadVector original_beads;
    std::vector<CAbstractBead*> id_to_original_bead;

    bool any_frozen_beads;

    uint64_t bead_round_tag_lcg; // Used to populate bead.round_tag for hashing rngs

    unsigned pos_to_cell_index(const float *pos)
    {
        int parts[3];
        for(int d=0; d<3; d++){
            parts[d] = std::floor(pos[d]);
            DEBUG_ASSERT( 0 <= parts[d] && parts[d] <dims_int[d]);
        }
        return pos_to_cell_index(parts);
    }

    unsigned pos_to_cell_index(const int *parts)
    {
        for(int d=0; d<3; d++){
            DEBUG_ASSERT( 0 <= parts[d] && parts[d] <dims_int[d]);
        }
        unsigned index = parts[0] + dims_int[0]*parts[1] + dims_int[0]*dims_int[1]*parts[2];
        DEBUG_ASSERT(index < cells.size());
        return index;
    }

    template<class TF, class TP>
    void for_each_cell(TF &&f, TP &&prefetch)
    {
        unsigned n=cells.size();
        if(USE_MORTON){
            for(unsigned i=0; i<n; i++){
                unsigned index1 = cell_enum_order[i];
                if(i+1 < n){
                    unsigned index2 = cell_enum_order[i+1];
                    prefetch( cells[index2] );
                }
                f( cells[index1] );
            }
        }else{
            for(unsigned i=0; i<n; i++){
                if(i+1<n){
                    prefetch( cells[i+1] );
                }
                f( cells[i] );
            }
        }
    }

    std::vector<std::array<int,3>> make_half_nhood()
    {
        std::set<std::array<int,3>> res;
        std::array<int,3> l;
        for(l[2]=-1; l[2] < 2; l[2]++){
            for(l[1]=-1; l[1] < 2; l[1]++){
                for(l[0]=-1; l[0] < 2; l[0]++){
                    if( l == std::array<int,3>{0,0,0} ){
                        continue; // Don't add centre
                    }
                    std::array<int,3> neg{ -l[0], -l[1], -l[2]};
                    if(res.find(neg)==res.end()){
                        res.insert(l);
                    }
                }
            }
        }
        require(res.size()==13);
        return {res.begin(), res.end()};
    }

    std::vector<std::array<int,3>> make_full_nhood()
    {
        std::set<std::array<int,3>> res;
        std::array<int,3> l;
        for(l[2]=-1; l[2] < 2; l[2]++){
            for(l[1]=-1; l[1] < 2; l[1]++){
                for(l[0]=-1; l[0] < 2; l[0]++){
                    if( l == std::array<int,3>{0,0,0} ){
                        continue; // Don't add centre
                    }
                    res.insert(l);
                }
            }
        }
        require(res.size()==26);
        return {res.begin(), res.end()};
    }

    virtual support_result on_import(ISimBox *box)
    {
        return {Supported};   
    }

    support_result import_all(ISimBox *box)
    {
        auto err=CanSupport(box);
        if(err.status!=Supported){
            return err;
        }

        global_seed=box->GetRNGSeed();
        round_id=box->GetCurrentTime();
        rng_stddev=sqrt(2 * CCNTCell::GetKt() / box->GetStepSize()); 

        dims_float[0] = box->GetSimBoxXLength();
        dims_float[1] = box->GetSimBoxYLength();
        dims_float[2] = box->GetSimBoxZLength();
        dims_float[3] = 0;
        for(int d=0; d<4; d++){
            dims_int[d] = (int)dims_float[d];
            half_dims_float[d] = dims_float[d] / 2;
        }
        unsigned ncells=dims_int[0] * dims_int[1] * dims_int[2];

        if(cells.size() != ncells){
            std::vector<Cell> buffer(ncells);
            std::swap(cells, buffer);
        }
        for(int d=0; d<3; d++){
            edge_adjustments[d][0] = 0;
            edge_adjustments[d][1] = +dims_float[d];
            edge_adjustments[d][2] = -dims_float[d];
        }

        dt=box->GetStepSize();
        halfdt=dt * 0.5f;
        halfdt2=halfdt * dt;

        const auto &bead_types=box->GetBeadTypes();
        num_bead_types = bead_types.size();
        coefficients.resize( num_bead_types * num_bead_types );
        for(unsigned i=0; i<num_bead_types; i++){
            for(unsigned j=0; j<num_bead_types; j++){
                coefficients[i*num_bead_types+j].first = box->GetDPDBeadConsInt(i, j);
                coefficients[i*num_bead_types+j].second = box->GetDPDBeadDissInt(i, j);
            }
        }

        std::vector<std::array<int,3>> rel_nhood;
        if( (NHOOD_TYPE & NhoodType_Size_Mask) == NhoodType_Size_Half){
            rel_nhood=make_half_nhood();
        }else if( (NHOOD_TYPE & NhoodType_Size_Mask) == NhoodType_Size_Full){
            rel_nhood=make_full_nhood();
        }else{
            require(0);
        }

        unsigned cell_index=0;
        std::array<int,3> cell_pos;
        for(cell_pos[2]=0; cell_pos[2]< dims_int[2]; cell_pos[2]++){
            for(cell_pos[1]=0; cell_pos[1]< dims_int[1]; cell_pos[1]++){
                for(cell_pos[0]=0; cell_pos[0] < dims_int[0]; cell_pos[0]++){
                    auto &cell = cells[cell_index];
                    for(int d=0; d<3; d++){
                        cell.origin[d] = cell_pos[d];
                    }
                    for(int d=3; d<CALC_DIM; d++){
                        cell.origin[d] = 0;
                    }
                    cell.cell_index = cell_index;
                    DEBUG_ASSERT( pos_to_cell_index(cell.origin) == cell_index );

                    cell.any_external=false;

                    unsigned neighbour_offset=0;
                    for(const auto &pos_delta : rel_nhood){
                        int neighbour_pos[3];
                        for(int d=0; d<3; d++){
                            neighbour_pos[d] = ( cell_pos[d] + pos_delta[d] + dims_int[d] ) % dims_int[d];
                        }

                        auto link = make_edge_tag( cell_pos.data(), neighbour_pos, dims_int );
                        cell.nhood[neighbour_offset] = link;
                        cell.any_external |= link.offsets[0] | link.offsets[1] | link.offsets[2];

                        neighbour_offset++;
                    }

                    cell.count = 0;
                    cell.to_move = 0;
                    cell.any_frozen = false;

                    ++cell_index;
                }
            }
        }

        // Build prefetch links. Final nhood entry is the first entry of
        // the following cell
        for(unsigned i=0; i<cells.size()-1; i++){
            cells[i].nhood[rel_nhood.size()] = cells[i].nhood[0];
        }
        cells.back().nhood[rel_nhood.size()] = cells.back().nhood[rel_nhood.size()];

        if(USE_MORTON){
            auto order=morton_codec::make_morton_order(dims_int[0], dims_int[1], dims_int[2]);
            cell_enum_order.resize(order.size());
            for(unsigned i=0; i<order.size(); i++){
                const auto &cp=order[i];
                int pos[3]={(int)std::get<0>(cp), (int)std::get<1>(cp), (int)std::get<2>(cp)};
                cell_enum_order[i] = pos_to_cell_index(pos);
            }
        }

        const auto &beads=box->GetBeads();
        unsigned num_beads=beads.size();
        bead_locations.assign(num_beads, -1); // We are assuming bead ids are contiguous.
        id_to_original_bead.assign(num_beads, nullptr);
        unpbc_flushed_fixups.assign(num_beads*3, 0);

        original_beads=box->GetBeads(); // This method will create a new vector.
        any_frozen_beads = false;

        std::set<long> seen_beed_ids;

        bead_round_tag_lcg = SplitMix64( box->GetRNGSeed() + SplitMix64(box->GetCurrentTime()) );
        for(CAbstractBead * b : original_beads){
            #ifndef NDEBUG
            auto it=seen_beed_ids.find(b->GetId());
            DEBUG_ASSERT(it==seen_beed_ids.end());
            seen_beed_ids.insert(b->GetId());
            #endif
            import_bead(*b, err);
            if(err.status!=Supported){
                return err;
            }
        }

        err = bonds.import_all(box);
        if(err.status!=Supported){
            return err;
        }

        return on_import(box);
    }

    double wrapped_distance(double x0, double x1, double w )
    {
        if(x0 < x1){
            return std::min(x1-x0, (x0+w)-x1);
        }else{
            return std::min(x0-x1, (x1+w)-x0);
        }
    }

    void import_bead(CAbstractBead &b, support_result &err)
    {
        // Osprey bead ids start at 1
        if( !(1 <= b.GetId()) ){
            err.status=PermanentProblem;
            err.reason="Bead ids do not start at 1.";
            return;
        }
        if(b.GetId() >= UINT32_MAX){
            err.status=PermanentProblem;
            err.reason="Bead id is more than 32-bit.";
            return;
        }
        DEBUG_ASSERT(0 <= b.GetType() && b.GetType() <= UINT8_MAX);

        Bead c;
        c.pos[0]=b.GetXPos();
        c.pos[1]=b.GetYPos();
        c.pos[2]=b.GetZPos();
        c.pos[3]=0;
        bool in_bounds = true;
        for(int d=0; d<3; d++){
            if(c.pos[d] == dims_float[d] ){ // Possible due to double->float rounding
                c.pos[d] = 0;  
            }
            require( 0 <= c.pos[d] && c.pos[d] < dims_float[d] );

            c.unpbcWrap[d] = 0;
            
            in_bounds |= (0 <= c.pos[0]) && (c.pos[0] < dims_float[d]);
        }
        if(!in_bounds){
            err.status=TransientProblemStep;
            err.reason="At least one bead is not currently in side the sim box";
            return;
        }

        #ifndef NDEBUG
            double unPBC=b.GetunPBCXPos(), PBC=b.GetXPos();
            DEBUG_ASSERT( wrapped_distance( PBC, fmod(unPBC , dims_float[0] ), dims_float[0]) < UNPBC_DRIFT_TOLERANCE );
        #endif

        if( wrapped_distance( b.GetXPos(), fmod(b.GetunPBCXPos() , dims_float[0] ), dims_float[0]) > UNPBC_DRIFT_TOLERANCE ){
            err.status=PermanentProblem; // Is anyone else going to fix this?
            err.reason="unPBC and PBC positions for one bead differ by more than tolerance";
            return;
        }

        c.mom[0]=b.GetXMom();
        c.mom[1]=b.GetYMom();
        c.mom[2]=b.GetZMom();
        c.mom[3]=0;
        c.force[0]=b.GetXForce();
        c.force[1]=b.GetYForce();
        c.force[2]=b.GetZForce();
        c.force[3]=0;
        c.bead_id=b.GetId() - 1; // Convert from 1-based to 0-based
        c.type=b.GetType();
        c.frozen=b.GetFrozen();
        any_frozen_beads |= c.frozen;
        if(rng_t::USES_BEAD_ROUND_TAG){
            c.round_tag = bead_round_tag_lcg >> 32;
            bead_round_tag_lcg = 6364136223846793005ull * bead_round_tag_lcg + 1;
        }

        for(int d=0; d<3; d++){
            c.mom[d] -= halfdt * c.force[d];
        }

        unsigned cell_index = pos_to_cell_index( c.pos );
        DEBUG_ASSERT(cell_index < cells.size());
        auto &cell = cells[cell_index];
        unsigned cell_offset = cell.count;
        if( cell_offset >= MAX_BEADS_PER_CELL){
            err.status = TransientProblemStep;
            err.reason = "Currently there are more than MAX_BEADS_PER_CELL in at least one cell";
            return;
        }
        cell.count += 1;
        cell.to_move += 1;
        cell.local[cell_offset] = c;
        cell.any_frozen |= c.frozen;

        unsigned bead_location=cell_index * MAX_BEADS_PER_CELL + cell_offset;
        if( c.bead_id >= bead_locations.size() ){
            err.status = PermanentProblem;
            err.reason = "Bead ids do not appear to be contiguous.";
            return;
        }
        if( bead_locations[c.bead_id] != -1){
            err.status = PermanentProblem;
            err.reason = "Two beads with the same id found.";
            return;
        }
        bead_locations[c.bead_id] = bead_location;

        require(id_to_original_bead[c.bead_id]==nullptr);
        id_to_original_bead[c.bead_id] = &b;
    }

    void export_all(ISimBox *box)
    { 
        #ifndef NDEBUG
        validate_cells();
        #endif

        auto simstate=box->GetSimBox();

        for(CAbstractBead * b : original_beads){
            unsigned id = b->GetId() - 1;
            auto &ib = find(id); // Not const, as we flush unPBC
            
            // Apply half mom correction to get out of loop skew
            b->SetXMom( ib.mom[0] + halfdt * ib.force[0] );
            b->SetYMom( ib.mom[1] + halfdt * ib.force[1] );
            b->SetZMom( ib.mom[2] + halfdt * ib.force[2] );

            b->SetXForce( ib.force[0] );
            b->SetYForce( ib.force[1] );
            b->SetZForce( ib.force[2] );

            for(int d=0; d<3; d++){
                DEBUG_ASSERT( 0<=ib.pos[d] && ib.pos[d]<dims_float[d] );
            }

            #ifndef NDEBUG
            {
                double unPBC=b->GetunPBCXPos(), PBC=b->GetXPos();
                DEBUG_ASSERT( wrapped_distance( PBC, fmod(unPBC , dims_float[0] ), dims_float[0]) < UNPBC_DRIFT_TOLERANCE );
            }
            #endif

            simstate->MoveBeadBetweenCNTCells(b, ib.pos[0], ib.pos[1], ib.pos[2]);
            // Recover the complete number of wraps
            int pbc_offset[3];
            for(int d=0; d<3; d++ ){
                pbc_offset[d] = ib.unpbcWrap[d] + unpbc_flushed_fixups[ib.bead_id*3+d];
            }
            
            // We take the PBC position as truth and snap the unPBC back to it.
            // For extremely long simulations this should be more accurate, as
            // if we do huge numbers of steps in here (e.g. 100000+), then the PBC positions will
            // be close to equilibrium for both DPD and bond forces. Trying to
            // snap back to unPBC would potentially jerk all the beads around
            // and disrupt that.
            b->SetunPBCXPos( b->GetunPBCXPos() - fmod(b->GetunPBCXPos(), dims_float[0]) + pbc_offset[0]*dims_float[0] + ib.pos[0] );
            b->SetunPBCYPos( b->GetunPBCYPos() - fmod(b->GetunPBCYPos(), dims_float[1]) + pbc_offset[1]*dims_float[1]+ ib.pos[1] );
            b->SetunPBCZPos( b->GetunPBCZPos() - fmod(b->GetunPBCZPos(), dims_float[2]) + pbc_offset[2]*dims_float[2]+ ib.pos[2] );

            // Check we preserve the invariant
            #ifndef NDEBUG
            {
                double unPBC=b->GetunPBCXPos(), PBC=b->GetXPos();
                DEBUG_ASSERT( wrapped_distance( PBC, fmod(unPBC , dims_float[0] ), dims_float[0]) < UNPBC_DRIFT_TOLERANCE );
            }
            #endif
        }
    }

    // Do sweep of all cells and validate beads. Pretty expensive.
    void validate_cells() const
    {
        std::unordered_set<long> ids;

        for(unsigned i=0; i<bead_locations.size(); i++){
            int32_t loc=bead_locations[i];
            require(0 <= loc);

            unsigned cell_index = loc / MAX_BEADS_PER_CELL;
            unsigned cell_off = loc % MAX_BEADS_PER_CELL;

            require(cell_index < cells.size());
            const auto &cell = cells[cell_index];

            const auto &bead = cell.local[cell_off];
            require(bead.bead_id == i);

            const auto &obead = find(i);
            require( &bead == &obead);

            require( ids.insert(i).second );
        }
    }

    const Bead &find(uint32_t bead_id) const
    {
        int32_t loc=bead_locations[bead_id];
        DEBUG_ASSERT(0 <= loc);

        unsigned cell_index=loc/MAX_BEADS_PER_CELL;
        unsigned cell_offset=loc%MAX_BEADS_PER_CELL;
        auto &cell=cells[cell_index];
        auto &res=cell.local[cell_offset];
        DEBUG_ASSERT(res.bead_id==bead_id);
        return res;
    }

    Bead &find(uint32_t bead_id)
    {
        int32_t loc=bead_locations[bead_id];
        DEBUG_ASSERT(0 <= loc);

        unsigned cell_index=loc/MAX_BEADS_PER_CELL;
        unsigned cell_offset=loc%MAX_BEADS_PER_CELL;
        auto &cell=cells[cell_index];
        auto &res=cell.local[cell_offset];
        DEBUG_ASSERT(res.bead_id==bead_id);
        return res;
    }

    bool calc_force(
        rng_t &rng,
        Bead &home,
        Bead &other,
        float other_x[4], // This includes any adjustment
        float newForce[4]  // The force acting on home bead
    ){
        float dx[4], dx2[4];

        for(int d=0; d<CALC_DIM; d++){
            dx[d] = home.pos[d] - other_x[d];
            DEBUG_ASSERT( std::abs(dx[d]) <= 2);  // Something has gone wrong if they are more than two apart

            dx2[d] = dx[d] * dx[d];
        }
	
		float dr2 = dx2[0] + dx2[1] + dx2[2];

        static const float min_r = 0.000000001f;
        static const float min_r2 = min_r * min_r;

        if( dr2 >= 1.0f || dr2 < min_r2)
        {		
            return false;
        }
        float dv[4], dx_dv[4];

        for(int d=0; d<3; d++){
            dv[d] = home.mom[d] - other.mom[d];
            dx_dv[d] = dx[d] * dv[d];
        }

        declare(dr2 > 0);
        float dr = std::sqrt(dr2);
        float wr = (1.0f - dr);
        float wr2 = wr*wr;
        float inv_dr = 1.0f/dr;

        const auto &coeffs = coefficients[ home.type * num_bead_types + other.type ];

        float conForce = coeffs.first * wr;

        float rdotv		= (dx_dv[0] + dx_dv[1] + dx_dv[2]) * inv_dr;
        float gammap	= coeffs.second * wr2;

        float dissForce	= -gammap*rdotv;
        declare(gammap > 0);
        float randForce	= std::sqrt(gammap)*rng.NextUnifScaled(home,other);
        float normTotalForce = (conForce + dissForce + randForce) * inv_dr;

        for(int d=0; d<CALC_DIM; d++){
            newForce[d] = normTotalForce * dx[d];
        }
        return true;
    }

    void calc_force(
        rng_t &rng,
        bool reflect_force,
        Bead &home,
        Bead &other,
        float other_x[4] // This includes any adjustment
    ){
        float newForce[4];

        if(calc_force(rng, home, other, other_x, newForce)){
            for(int d=0; d<CALC_DIM; d++){
                home.force[d] += newForce[d];
                if( reflect_force ){
                    other.force[d] -= newForce[d];
                }
            }
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
                calc_force(rng, /* reflect_force */ true, home_cell.local[i0], home_cell.local[i1], home_cell.local[i1].pos );
            }
        }

        for(unsigned nhood_index=0; nhood_index<GetNhoodSize(home_cell); nhood_index++){
            Cell &prefetch_cell = cells[home_cell.nhood[nhood_index+1].index];
            // Get the header of the cell
            __builtin_prefetch(&prefetch_cell.count);
            // First bead - note that it is offset from the cell header by quite a distance due to nhood links
            __builtin_prefetch(&prefetch_cell.local[0]);
            
            auto nlink=home_cell.nhood[nhood_index];
            Cell &other_cell=cells[nlink.index];
            bool reflect_forces = (NHOOD_TYPE & NhoodType_ForceShare_Mask) == NHoodType_ForceShare_Direct;

            DEBUG_ASSERT( nlink.tag==EdgeTag_ForceShare_Direct || nlink.tag==EdgeTag_ForceShare_None );
            reflect_forces |= nlink.tag==EdgeTag_ForceShare_Direct;
            
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
                    calc_force(rng, reflect_forces, home_bead, other_bead, AnyExternal ? other_x : other_bead.pos);
                }
            }
        }
    }

    virtual __attribute__((noinline)) void update_forces()
    {
        rng_t rng(rng_stddev, global_seed, round_id, 0);

        // TODO : branching based on any external may introduce data-dependent
        // control and bloat instruction size. Is minor saving worth it?
        for_each_cell([&](Cell &cell){
            if(cell.any_external){
                update_forces<true>(rng, cell);
            }else{
                update_forces<false>(rng, cell);
            }
        },[](Cell &cell){

        });
    }

    virtual __attribute__((noinline)) void update_mom_and_move()
    {
        rng_t rng(rng_stddev, global_seed, round_id, 1);

        if(any_frozen_beads){
            for_each_cell([&](Cell &cell){
                if(cell.any_frozen){
                    update_mom_and_move<true>(rng, cell);
                }else{
                    update_mom_and_move<false>(rng, cell);
                }
            },[&](Cell &next){
                __builtin_prefetch(&next.count);
                __builtin_prefetch(&next.local[0]);
            });
        }else{
            for_each_cell([&](Cell &cell){
                update_mom_and_move<false>(rng, cell);
            },[&](Cell &next){
                __builtin_prefetch(&next.count);
                __builtin_prefetch(&next.local[0]);
            });
        }
    }

    void flush_unpbc_wrap(Bead &ib)
    {
        for(int d=0; d<3; d++){
            unpbc_flushed_fixups[ib.bead_id*3+d] += ib.unpbcWrap[d];
            ib.unpbcWrap[d] = 0;
        }
    }

    template<bool AnyFrozen>
    void update_mom_and_move(rng_t &rng, Cell &cell)
    {
        DEBUG_ASSERT(AnyFrozen ? any_frozen_beads : true);
        
        unsigned i=0;
        while(i<cell.to_move){
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
                bead.mom[d] += halfdt * bead.force[d];
                bead.pos[d] += dt * bead.mom[d] + halfdt2 * bead.force[d];
		        // We know that lambda=0.5
                bead.mom[d] += halfdt * bead.force[d];
        		bead.force[d] = 0;

                float new_origin_d = std::floor( bead.pos[d] );
                moved |= new_origin_d != cell.origin[d];

                if(rng_t::USES_BEAD_ROUND_TAG){
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
                        flush_unpbc_wrap(bead);
                    }
                }
                // The above could results in bead.pos[d] == dims_float[d] due to
                // rounding, so this is not an else if
                if(bead.pos[d] >= dims_float[d]){
                    bead.pos[d] -= dims_float[d];
                    bead.unpbcWrap[d] += 1;
                    if( bead.unpbcWrap[d] > +100 ){
                        flush_unpbc_wrap(bead);
                    }
                }
                new_origin[d] = std::floor( bead.pos[d] );
            }

            unsigned new_index = pos_to_cell_index(new_origin);
            Cell &new_cell = cells[new_index];
            memcpy( new_cell.local+new_cell.count, &bead, sizeof(Bead) );
            bead_locations[bead.bead_id] = new_index*MAX_BEADS_PER_CELL + new_cell.count;
            new_cell.count++;

            // We have created a gap at index i

            if(i+1 == cell.count){
                // The gap is the last thing in the array
                DEBUG_ASSERT(i+1 == cell.to_move && cell.to_move == cell.count);
                cell.count -= 1;
                // We are finished
                break;
            }
            
            unsigned src;
            // Either there is another bead todo, or an incoming bead which we should not do
            if(cell.to_move < cell.count){
                // A new bead has been moved in. Use it to fill the gap and advance i
                cell.count -= 1;
                
                memcpy(cell.local+i, cell.local+cell.count, sizeof(Bead));
                bead_locations[cell.local[i].bead_id] = cell.cell_index*MAX_BEADS_PER_CELL + i;

                // We need to move past this bead, as it is already processed in another cell
                ++i;
            }else{
                DEBUG_ASSERT(i+1 < cell.to_move && cell.to_move == cell.count); 
                // We still have at least one more bead, but nothing new
                // Use the last todo bead to fill the gap and don't advance i
                cell.to_move -= 1;
                cell.count -= 1;

                memcpy(cell.local+i, cell.local+cell.to_move, sizeof(Bead));
                bead_locations[cell.local[i].bead_id] = cell.cell_index*MAX_BEADS_PER_CELL + i;

                // Don't change i, as we still have to process this bead that now fills the slot
            }            
        }
    
        cell.to_move=0;
    }

    /*
    void update_mom_and_move_vec4(Cell &cell)
    {

        unsigned i=0;
        while(i<cell.to_move){
            Bead &bead=cell.local[i];
            
            __m128 pos = _mm_load_ps(bead.pos);
            __m128 mom = _mm_load_ps(bead.mom);
            __m128 force = _mm_load_ps(bead.force);

            __m128 delta_mom = _mm_set1_ps(halfdt) * force;
            mom = mom + delta_mom;
            pos = pos + _mm_set1_ps(dt) * mom + _mm_set1_ps(halfdt2) * force;
            mom = mom + delta_mom;
            force = _mm_set1_ps(0.0f);

            __m128 new_origin_f = _mm_floor_ps(pos);
            __m128 any_moved = _mm_cmpeq_ps(new_origin_f, _mm_load_ps(cell.origin));

            if(_mm_test_all_ones(_mm_castps_si128(any_moved))){
                _mm_store_ps(bead.pos, pos);
                _mm_store_ps(bead.mom, mom);
                _mm_store_ps(bead.force, force);
                ++i;
                continue;
            }

            __m128i metadata = _mm_load_si128( (__m128i*) &bead.bead_id );

            __m128 dims = _mm_load_ps(dims_float);
            __m128 off_low = _mm_cmplt_ps(pos, _mm_setzero_ps());
            pos = pos + _mm_and_ps( off_low, dims );

            __m128 off_high = _mm_cmpge_ps(pos, dims);
            pos = pos - _mm_and_ps( off_high, dims );

            new_origin_f = _mm_floor_ps(pos);
            __m128i new_origin_i = _mm_cvtps_epi32(new_origin_f);

            int32_t new_origin[4];
            _mm_store_si128((__m128i*)new_origin, new_origin_i);

            unsigned new_index = pos_to_cell_index(new_origin);
            Cell &new_cell = cells[new_index];
            Bead &new_bead=new_cell.local[new_cell.count];
            _mm_store_ps(new_bead.pos, pos);
            _mm_store_ps(new_bead.mom, mom);
            _mm_store_ps(new_bead.force, force);
            _mm_store_si128( (__m128i*)&new_bead.bead_id, metadata);
            bead_locations[bead.bead_id] = new_index*MAX_BEADS_PER_CELL + new_cell.count;
            new_cell.count++;

            // We have created a gap at index i

            // TODO: Below is probably wrong
            cell.count -= 1;
            if(i == cell.count){
                // The gap is the last thing in the array
                DEBUG_ASSERT(i+1 == cell.to_move && cell.to_move == cell.count+1);
                // We are finished
                break;
            }else{
                unsigned src;
                if(cell.to_move < cell.count){
                    // A new bead has been moved in. Use it to fill the gap
                    src = cell.count;
                }else{
                    DEBUG_ASSERT(i+1 < cell.to_move && cell.to_move == cell.count+1); 
                    // We still have at least one more bead, but nothing new
                    // Use the last todo bead to fill the gap
                    cell.to_move -= 1;
                    src = cell.to_move;
                }
                memcpy( cell.local+i,  cell.local+src, sizeof(Bead));
                bead_locations[cell.local[src].bead_id] = cell.cell_index * MAX_BEADS_PER_CELL + i;
                ++i;
            }
        }
    
        cell.to_move=0;
    }
    */
};

#endif
