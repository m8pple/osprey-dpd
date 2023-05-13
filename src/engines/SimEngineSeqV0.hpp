#ifndef SimEngineSeqV0_hpp
#define SimEngineSeqV0_hpp

#include <cstdint>
#include <vector>
#include <atomic>
#include <vector>
#include <mutex>
#include <cassert>
#include <cstdlib>
#include <set>
#include <cstring>

#include "AbstractBead.h"
#include "ISimBox.h"
#include "xxBase.h"
#include "CNTCell.h"
#include "Polymer.h"
#include "Bond.h"
#include "BondPair.h"

#include "ISimEngine.h"

#include "immintrin.h"
#include "pmmintrin.h"

#include "morton_codec.hpp"

#include "bond_info.hpp"

static bool require_fail_impl(const char *file, int line, const char *cond)
{
    fprintf(stderr, "%s:%d: requirement failed : %s\n", file, line, cond);
    assert(0);
    exit(1);
}

#define require(cond) \
    if(!(cond)){ require_fail_impl(__FILE__, __LINE__, #cond); }

static void declare(bool cond)
{
    if(!cond){
        assert(0);
        __builtin_unreachable();
    }
}


struct SimEngineSeqV0
    : public ISimEngine
{
public:
    std::string Name() const override
    {
        return "SimEngineSeqV0";
    }

    bool IsParallel() const override
    { return false; }



    std::string CanSupport(const ISimBox *box) const override
    {
        if(box->GetSimBoxXOrigin()!=0 || box->GetSimBoxXOrigin()!=0 || box->GetSimBoxZOrigin()!=0){
            return "Sim box origin is not at zero."; // Actually fairly easy to support
        }
        if( std::round(box->GetSimBoxXLength()) != box->GetSimBoxXLength() || std::round(box->GetSimBoxYLength()) != box->GetSimBoxYLength() || std::round(box->GetSimBoxZLength()) != box->GetSimBoxZLength()){
            return "Sim box does not have integer dimensions.";
        }
        if(box->GetLambda() != 0.5){
            return "Lambda != 0.5";
        }

        return {};
    }

    void Run(ISimBox *box, bool modified, unsigned num_steps) override
    {
        import_all(box);

        auto bead_source=[&](uint32_t bead_id) -> Bead &{
            return find(bead_id);
        };

        #ifndef NDEBUG
        validate_cells();
        #endif

        for(unsigned i=0; i<num_steps; i++){
            update_mom_and_move();
            m_generation += 1;

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
        }

        #ifndef NDEBUG
        validate_cells();
        #endif

        export_all(box);
    }

protected:

    static const int MAX_BEADS_PER_CELL = 32;

    using calc_t = double;

    struct Bead
    {
        calc_t pos[3];
        calc_t mom[3];
        calc_t force[3];
        uint32_t bead_id;
        uint8_t type;
        uint8_t frozen;
        int16_t unpbcWrap[3];
        uint32_t generation;
    };
    //static_assert(sizeof(Bead)==4*4*4);

    struct Neighbour
    {
        uint32_t index;
        std::array<uint8_t,4> edge_tags; // Tells us what boundary to apply to positions in this cell
    };

    std::array<uint8_t,4> make_edge_tag(const int32_t *home, const int32_t *neighbour, const int32_t *dims)
    {
        std::array<uint8_t,4> local;
        for(int d=0; d<3; d++){
            int delta=neighbour[d] - home[d];
            if(delta==0 || delta==1 || delta==-1){
                local[d]=0;
            }else if(delta==1-dims[d]){
                local[d]=1;
            }else if(delta==dims[d]-1){
                local[d]=2;
            }else{
                assert(0);
                fprintf(stderr, "Logic violation in make_edge_tag");
                exit(1);
            }
        }
        return local;
    }

    struct Cell
    {
        uint32_t count;
        uint32_t to_move; // to_move <= count

        calc_t origin[3];
        uint32_t cell_index;

        Neighbour nhood[13];

        Bead local[MAX_BEADS_PER_CELL];
    };

    calc_t dims_float[4];
    int32_t dims_int[4];
    calc_t half_dims_float[4];
    calc_t dt;
    calc_t halfdt;
    calc_t halfdt2;
    
    calc_t edge_adjustments[3][3];

    uint32_t num_bead_types;
    std::vector<std::pair<calc_t,calc_t>> coefficients;

    std::vector<Cell> cells;
    std::vector<int32_t> bead_locations; // The current location of each slot

    BondInfo<calc_t> bonds;

    AbstractBeadVector original_beads;

    uint32_t m_generation;

    uint64_t rng_state = 0;
    calc_t rng_scale;

    unsigned pos_to_cell_index(const calc_t *pos)
    {
        int parts[3];
        for(int d=0; d<3; d++){
            parts[d] = std::floor(pos[d]);
            assert( 0 <= parts[d] && parts[d] <dims_int[d]);
        }
        return pos_to_cell_index(parts);
    }

    unsigned pos_to_cell_index(const int *parts)
    {
        for(int d=0; d<3; d++){
            assert( 0 <= parts[d] && parts[d] <dims_int[d]);
        }
        unsigned index = parts[0] + dims_int[0]*parts[1] + dims_int[0]*dims_int[1]*parts[2];
        assert(index < cells.size());
        for(int d=0; d<3; d++){
            assert(cells[index].origin[d]==parts[d]);
        }
        return index;
    }

    template<class TF, class TP>
    void for_each_cell(TF &&f, TP &&prefetch)
    {
        unsigned n=cells.size();
        for(unsigned i=0; i<n; i++){
            if(i+1<n){
                prefetch( cells[i+1] );
            }
            f( cells[i] );
        }
    }

    std::vector<std::array<int,3>> make_forward_nhood()
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

    void import_all(ISimBox *box)
    {
        auto err=CanSupport(box);
        if(!err.empty()){
            fprintf(stderr, "%s\n", err.c_str());
            exit(1);
        }

        if(rng_state==0){
            rng_state=box->GetRNGSeed();
        }

        dims_float[0] = box->GetSimBoxXLength();
        dims_float[1] = box->GetSimBoxYLength();
        dims_float[2] = box->GetSimBoxZLength();
        dims_float[3] = 0;
        for(int d=0; d<4; d++){
            dims_int[d] = (int)dims_float[d];
            half_dims_float[d] = dims_float[d] / 2;
        }
        unsigned ncells=dims_int[0] * dims_int[1] * dims_int[2];
        cells.resize( ncells );
        for(int d=0; d<3; d++){
            edge_adjustments[d][0] = 0;
            edge_adjustments[d][1] = +dims_float[d];
            edge_adjustments[d][2] = -dims_float[d];
        }

        dt=box->GetStepSize();
        halfdt=dt * 0.5;
        halfdt2=halfdt * dt;

        double invrootdt = sqrt(24.0*CCNTCell::GetKt()/dt);
        rng_scale = ldexp(invrootdt, -32);

        const auto &bead_types=box->GetBeadTypes();
        num_bead_types = bead_types.size();
        coefficients.resize( num_bead_types * num_bead_types );
        for(unsigned i=0; i<num_bead_types; i++){
            for(unsigned j=0; j<num_bead_types; j++){
                coefficients[i*num_bead_types+j].first = box->GetDPDBeadConsInt(i, j);
                coefficients[i*num_bead_types+j].second = box->GetDPDBeadDissInt(i, j);
            }
        }

        std::vector<std::array<int,3>> forward_nhood=make_forward_nhood();

        unsigned cell_index=0;
        std::array<int,3> cell_pos;
        for(cell_pos[2]=0; cell_pos[2]< dims_int[2]; cell_pos[2]++){
            for(cell_pos[1]=0; cell_pos[1]< dims_int[1]; cell_pos[1]++){
                for(cell_pos[0]=0; cell_pos[0] < dims_int[0]; cell_pos[0]++){
                    auto &cell = cells[cell_index];
                    for(int d=0; d<3; d++){
                        cell.origin[d] = cell_pos[d];
                    }
                    cell.cell_index = cell_index;
                    assert( pos_to_cell_index(cell.origin) == cell_index );
                    cell_index += 1;
                }
            }
        }

        cell_index=0;
        for(cell_pos[2]=0; cell_pos[2]< dims_int[2]; cell_pos[2]++){
            for(cell_pos[1]=0; cell_pos[1]< dims_int[1]; cell_pos[1]++){
                for(cell_pos[0]=0; cell_pos[0] < dims_int[0]; cell_pos[0]++){
                    auto &cell = cells[cell_index];
                    unsigned neighbour_offset=0;
                    for(const auto &pos_delta : forward_nhood){
                        int neighbour_pos[3];
                        for(int d=0; d<3; d++){
                            neighbour_pos[d] = ( cell_pos[d] + pos_delta[d] + dims_int[d] ) % dims_int[d];
                        }

                        cell.nhood[neighbour_offset].index=pos_to_cell_index(neighbour_pos);
                        auto tags=make_edge_tag( cell_pos.data(), neighbour_pos, dims_int );
                        cell.nhood[neighbour_offset].edge_tags = tags;

                        neighbour_offset++;
                    }

                    cell.count = 0;
                    cell.to_move = 0;

                    ++cell_index;
                }
            }
        }

        const auto &beads=box->GetBeads();
        unsigned num_beads=beads.size();
        bead_locations.assign(num_beads, -1); // We are assuming bead ids are contiguous.

        original_beads=box->GetBeads(); // This method will create a new vector.
        for(CAbstractBead * b : original_beads){
            import_bead(*b);
        }

        m_generation=0;
        bonds.import_all(box);
    }

    void import_bead(CAbstractBead &b)
    {
        // Osprey bead ids start at 1
        require(1 <= b.GetId() && b.GetId() <= UINT32_MAX);
        require(0 <= b.GetType() && b.GetType() <= UINT8_MAX);

        Bead c;
        c.pos[0]=b.GetXPos();
        c.pos[1]=b.GetYPos();
        c.pos[2]=b.GetZPos();
        for(int d=0; d<3; d++){
            if(c.pos[d] == dims_float[d] ){ // Possible due to double->float rounding
                c.pos[d] = 0;  
            }
            require( 0 <= c.pos[d] && c.pos[d] < dims_float[d] );
        }
        b.SetunPBCXPos( b.GetunPBCXPos() - c.pos[0] );
        b.SetunPBCYPos( b.GetunPBCYPos() - c.pos[1] );
        b.SetunPBCZPos( b.GetunPBCZPos() - c.pos[2] );

        c.mom[0]=b.GetXMom();
        c.mom[1]=b.GetYMom();
        c.mom[2]=b.GetZMom();
        c.force[0]=b.GetXForce();
        c.force[1]=b.GetYForce();
        c.force[2]=b.GetZForce();
        c.bead_id=b.GetId() - 1; // Convert from 1-based to 0-based
        c.type=b.GetType();
        c.frozen=b.GetFrozen();
        for(int d=0; d<3; d++){
            c.unpbcWrap[d] = 0;
        }
        c.generation=0;

        for(int d=0; d<3; d++){
            c.mom[d] -= halfdt * c.force[d];
        }

        unsigned cell_index = pos_to_cell_index( c.pos );
        auto &cell = cells[cell_index];
        unsigned cell_offset = cell.count;
        require( cell_offset < MAX_BEADS_PER_CELL );
        cell.count += 1;
        cell.to_move += 1;
        cell.local[cell_offset] = c;

        unsigned bead_location=cell_index * MAX_BEADS_PER_CELL + cell_offset;
        require( c.bead_id < bead_locations.size() && bead_locations[c.bead_id] == -1);
        bead_locations[c.bead_id] = bead_location;
    }

    void export_all(ISimBox *box)
    { 
        auto simstate=box->GetSimBox();

        for(CAbstractBead * b : original_beads){
            unsigned id = b->GetId() - 1;
            const auto &ib = find(id);
            
            // Apply half mom correction to get out of loop skew
            b->SetXMom( ib.mom[0] + halfdt * ib.force[0] );
            b->SetYMom( ib.mom[1] + halfdt * ib.force[1] );
            b->SetZMom( ib.mom[2] + halfdt * ib.force[2] );

            b->SetXForce( ib.force[0] );
            b->SetYForce( ib.force[1] );
            b->SetZForce( ib.force[2] );

            b->SetunPBCXPos( b->GetunPBCXPos() + ib.pos[0] + ib.unpbcWrap[0] * dims_float[0] );
            b->SetunPBCYPos( b->GetunPBCYPos() + ib.pos[1] + ib.unpbcWrap[1] * dims_float[1] );
            b->SetunPBCZPos( b->GetunPBCZPos() + ib.pos[2] + ib.unpbcWrap[2] * dims_float[2] );

            simstate->MoveBeadBetweenCNTCells(b, ib.pos[0], ib.pos[1], ib.pos[2]);
        }
    }

    // Do sweep of all cells and validate beads. Pretty expensive.
    void validate_cells() const
    {
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
        }
    }

    const Bead &find(uint32_t bead_id) const
    {
        int32_t loc=bead_locations[bead_id];
        assert(0 <= loc);

        unsigned cell_index=loc/MAX_BEADS_PER_CELL;
        unsigned cell_offset=loc%MAX_BEADS_PER_CELL;
        auto &cell=cells[cell_index];
        auto &res=cell.local[cell_offset];
        assert(res.bead_id==bead_id);
        return res;
    }

    Bead &find(uint32_t bead_id)
    {
        int32_t loc=bead_locations[bead_id];
        assert(0 <= loc);

        unsigned cell_index=loc/MAX_BEADS_PER_CELL;
        unsigned cell_offset=loc%MAX_BEADS_PER_CELL;
        auto &cell=cells[cell_index];
        auto &res=cell.local[cell_offset];
        assert(res.bead_id==bead_id);
        return res;
    }

    calc_t RandUnifScaled()
    {
        uint32_t u32 = rng_state>>32;
        rng_state=6364136223846793005ull * rng_state + 1;
        
        int32_t i31;
    	memcpy(&i31, &u32, 4); // Avoid undefined behaviour. Gets number in range [-2^31,2^31)
	    return i31 * rng_scale; // rng_scale = ( CCNTCell::m_invrootdt * 2^-32  )
    }

    void calc_forces(
        Bead &home,
        Bead &other,
        calc_t other_x[4] // This includes any adjustment
    ){
        calc_t dx[4], dx2[4];

        for(int d=0; d<3; d++){
            dx[d] = home.pos[d] - other_x[d];
            assert( std::abs(dx[d]) <= 2);  // Something has gone wrong if they are more than two apart

            dx2[d] = dx[d] * dx[d];
        }
	
		calc_t dr2 = dx2[0] + dx2[1] + dx2[2];

        static const calc_t min_r = calc_t(0.000000001);
        static const calc_t min_r2 = min_r * min_r;

        if( dr2 < calc_t(1.0) && dr2 > min_r2)
        {		
            calc_t dv[4], dx_dv[4];

            for(int d=0; d<3; d++){
                dv[d] = home.mom[d] - other.mom[d];
                dx_dv[d] = dx[d] * dv[d];
            }

            declare(dr2 > 0);
            calc_t dr = std::sqrt(dr2);
            calc_t wr = (calc_t(1.0) - dr);
            calc_t wr2 = wr*wr;
            calc_t inv_dr = calc_t(1.0)/dr;

            const auto &coeffs = coefficients[ home.type * num_bead_types + other.type ];

            calc_t conForce = coeffs.first * wr;

            calc_t rdotv		= (dx_dv[0] + dx_dv[1] + dx_dv[2]) * inv_dr;
            calc_t gammap	= coeffs.second * wr2;

            calc_t dissForce	= -gammap*rdotv;
            declare(gammap > 0);
            calc_t rv = RandUnifScaled();

            /*
            static double rvSum=0;
            static double rvSumSqr=0;
            static unsigned rvNum=0;
            rvSum += rv;
            rvSumSqr += rv*rv;
            rvNum += 1;
            if(rvNum==1000000){
                fprintf(stderr, "Mean = %f, Stddev=%f\n", rvSum/rvNum, sqrt(rvSumSqr/rvNum));
            }
            */


            calc_t randForce	= std::sqrt(gammap)*rv;
            calc_t normTotalForce = (conForce + dissForce + randForce) * inv_dr;

            calc_t newForce[4];
            for(int d=0; d<3; d++){
                newForce[d] = normTotalForce * dx[d];
                home.force[d] += newForce[d];
                other.force[d] -= newForce[d];
            }
        }
    }

    void update_forces(Cell &home_cell)
    {
        home_cell.to_move = home_cell.count;

        if(home_cell.count==0){
            return;
        }

        for(unsigned i0=0; i0<home_cell.count-1; i0++){
            for(unsigned i1=i0+1; i1<home_cell.count; i1++){
                calc_forces( home_cell.local[i0], home_cell.local[i1], home_cell.local[i1].pos );
            }
        }

        for(unsigned nhood_index=0; nhood_index<13; nhood_index++){
            auto nlink=home_cell.nhood[nhood_index];
            Cell &other_cell=cells[nlink.index];

            calc_t other_delta[4];
            for(int d=0; d<3; d++){
                other_delta[d] = edge_adjustments[d][nlink.edge_tags[d]];
            }
            other_delta[3]=0;
            
            for(unsigned other_i=0; other_i<other_cell.count; other_i++){
                Bead &other_bead = other_cell.local[other_i];
                calc_t other_x[4];
                for(int d=0; d<3; d++){
                    other_x[d] = other_bead.pos[d] + other_delta[d];
                }
                
                for(unsigned home_i=0; home_i<home_cell.count; home_i++){
                    Bead &home_bead = home_cell.local[home_i];

                    calc_forces(home_bead, other_bead, other_x);
                }
            }
        }
    }

    virtual __attribute__((noinline)) void update_forces()
    {
        // TODO : branching based on any external may introduce data-dependent
        // control and bloat instruction size. Is minor saving worth it?
        for(Cell &cell : cells){
            update_forces(cell);
        }
    }

    virtual __attribute__((noinline)) void update_mom_and_move()
    {
        for(auto &cell : cells){
            update_mom_and_move(cell);
        }
    }

    void update_mom_and_move(Cell &cell)
    {
        for(unsigned i=0; i<cell.to_move; i++){
            assert(cell.local[i].generation==m_generation);
        }
        for(unsigned i=cell.to_move; i<cell.count; i++){
            assert(cell.local[i].generation==m_generation+1);
        }

        unsigned i=0;
        while(i<cell.to_move){
            Bead &bead=cell.local[i];

            assert(bead.generation==m_generation);
            bead.generation += 1;

            if(bead.frozen){
                for(int d=0; d<3; d++){
                    bead.force[d] = 0;
                }
                ++i;
                continue;
            }

            bool moved=false;

            calc_t orig_force[4];

            for(int d=0; d<3; d++){
                orig_force[d] = bead.force[d];

                // Mom from previous step (loop skewed)
                bead.mom[d] += halfdt * bead.force[d];
                bead.pos[d] += dt * bead.mom[d] + halfdt2 * bead.force[d];
		        // We know that lambda=0.5
                bead.mom[d] += halfdt * bead.force[d];
        		bead.force[d] = 0;

                calc_t new_origin_d = std::floor( bead.pos[d] );
                moved |= new_origin_d != cell.origin[d];
            }

            if(!moved){
                ++i;
                continue;
            }


            calc_t pre_wrap_pos[4];
            int new_origin[4];
            for(int d=0; d<3; d++){
                pre_wrap_pos[d] = bead.pos[d];
                if(bead.pos[d] < 0){
                    bead.pos[d] += dims_float[d];
                    bead.unpbcWrap[d] -= 1;
                }
                // The above could results in bead.pos[d] == dims_float[d] due to
                // rounding, so this is not an else if
                if(bead.pos[d] >= dims_float[d]){
                    bead.pos[d] -= dims_float[d];
                    bead.unpbcWrap[d] += 1;
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
                assert(i+1 == cell.to_move && cell.to_move == cell.count);
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
                assert(i+1 < cell.to_move && cell.to_move == cell.count); 
                // We still have at least one more bead, but nothing new
                // Use the last todo bead to fill the gap and don't advance i
                cell.to_move -= 1;
                cell.count -= 1;

                memcpy(cell.local+i, cell.local+cell.to_move, sizeof(Bead));
                bead_locations[cell.local[i].bead_id] = cell.cell_index*MAX_BEADS_PER_CELL + i;

                // Don't change i, as we still have to process this bead that now fills the slot
            }            
        }

        for(unsigned i=0; i<cell.count; i++){
            assert(cell.local[i].generation==m_generation+1);
        }
    
        cell.to_move=0;
    }

};

#endif
