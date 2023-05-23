#include "SimEngineSeq.hpp"


#if ( defined(_M_ARM64) || defined(__aarch64__) )

#include "arm_neon.h"

struct  SimEngineSeqNEONPolicy
    : EnginePolicyConcept
{
};

class SimEngineSeqNEON
    : public SimEngineSeq<SimEngineSeqNEONPolicy>
{
public:
    std::string Name() const override
    {
        return "SimEngineSeqNEON";
    }

private:

    using base_t = SimEngineSeq<SimEngineSeqNEONPolicy>;
    using rng_t = typename base_t::rng_t;

    uint64x2_t rng_state_neon;
    float rng_stddev_neon;

    support_result on_import(ISimBox *box) override 
    {
        uint64_t rng=SplitMix64(SplitMix64(box->GetCurrentTime()) + box->GetRNGSeed() + 21313112);

        rng_state_neon[0]= rng;
        rng_state_neon[1]= SplitMix64(rng);

        rng_stddev_neon=CalcRNGScaleForU31(base_t::rng_stddev);
        
        return {Supported};
    }

    float32x4_t RandUnifScaledNEON(uint64x2_t &state)
    {
        /*
        x ^= x << 13;
	    x ^= x >> 7;
	    x ^= x << 17;
        */

        uint64x2_t orig=state;

        state = state ^ (state << 13);
        state = state ^ (state >> 7);
        state = state ^ (state << 17);

        auto u = state; // orig + _mm256_shuffle_epi32( state, (2<<6)|(3<<4)|(0<<2)|(1<<0));

        auto ui32 = vreinterpretq_s32_u64(u);
        auto uf32 = vcvtq_f32_s32(ui32);
        return vmulq_f32( uf32 , vdupq_n_f32( rng_stddev_neon ));
    }

    bool any_nan(float32x4_t x)
    {
        auto is_non_nan = vceqq_f32(x,x);

        for(unsigned i=0; i<4; i++){
            if(!is_non_nan[i]){
                return true;
            }
        }
        return false;
    }

    virtual __attribute__((noinline)) void update_forces() override
    {
        rng_t rng(rng_stddev, global_seed, round_id, 0);

        // TODO : branching based on any external may introduce data-dependent
        // control and bloat instruction size. Is minor saving worth it?
        base_t::for_each_cell([&](Cell &cell){
            if(cell.any_external){
                update_forces<true>(rng, cell);
            }else{
                update_forces<false>(rng, cell);
            }
        },[](Cell &cell){

        });
    }

    void calc_force_neon(
        uint64x2_t &rng_state,

        unsigned home_count,
        const float32x4_t home_pos[3],
        const float32x4_t home_mom[3],
        float32x4_t home_force[3],
        uint32x4_t home_type,

        Bead &other,
        float other_x[4]
    ) {
        DEBUG_ASSERT(home_count <= 4);

        const float32x4_t ONE = vdupq_n_f32(1.0f);

        float32x4_t dx[3], dx2[3];

        for(int d=0; d<3; d++){
            dx[d] = home_pos[d] - vdupq_n_f32(other_x[d]);
            dx2[d] = dx[d] * dx[d];
        }
	
		float32x4_t dr2 = dx2[0] + dx2[1] + dx2[2];

        static const float min_r = 0.000000001f;
        static const float min_r2 = min_r * min_r;

        uint32x4_t active =  vcgtq_f32(dr2, vdupq_n_f32(min_r2)) & vcltq_f32(dr2, ONE );

        if(vmaxvq_u32(active)==0){
            return;
        }

        float32x4_t dv[3], dx_dv[3];

        for(int d=0; d<3; d++){
            dv[d] = home_mom[d] - vdupq_n_f32(other.mom[d]);
            dx_dv[d] = dx[d] * dv[d];
        }

        auto dr = vsqrtq_f32(dr2);
        auto wr = ONE - dr;
        auto wr2 = wr*wr;
        auto inv_dr = vrecpeq_f32(dr);
        inv_dr = vmulq_f32(vrecpsq_f32(dr, inv_dr), inv_dr);
        inv_dr = vmulq_f32(vrecpsq_f32(dr, inv_dr), inv_dr); // TODO: Can we get away with just one iteration?

        const float *coeff_row = (float*)&coefficients[other.type * num_bead_types]; // Note we are type-punning from std::pair<float> to float
        
        float32x4_t conservative_coeff=ONE;
        float32x4_t dissipative_coeff=ONE;
        for(int ci=0; ci<home_count; ci++){
            DEBUG_ASSERT(home_type[ci] < num_bead_types);
            conservative_coeff[ci] = coeff_row[home_type[ci]*2];
            dissipative_coeff[ci] = coeff_row[home_type[ci]*2+1];
        }

        auto conForce = conservative_coeff * wr;

        auto rdotv = (dx_dv[0] + dx_dv[1] + dx_dv[2]) * inv_dr;
        auto gammap	= dissipative_coeff * wr2;

        auto dissForce	= -gammap*rdotv;
        auto randForce	= vsqrtq_f32(gammap) * RandUnifScaledNEON(rng_state);
        auto normTotalForce = (conForce + dissForce + randForce) * inv_dr;

        normTotalForce = vreinterpretq_f32_u32( vandq_u32( vreinterpretq_u32_f32( normTotalForce) , active));

        float32x4_t newForce[3];
        for(int d=0; d<3; d++){
            newForce[d] = normTotalForce * dx[d];
            DEBUG_ASSERT(!any_nan(newForce[d]));

            home_force[d] += newForce[d];

            other.force[d] -= vaddvq_f32(newForce[d]);

            DEBUG_ASSERT(std::abs(other.force[d]) < 10000);
        }
    }

    template<bool AnyExternal>
    void update_forces(rng_t &base_rng, Cell &home_cell)
    {
        home_cell.to_move = home_cell.count;

        if(home_cell.count==0){
            return;
        }

        if(home_cell.count > 4){
            return SimEngineSeq::update_forces<AnyExternal>(base_rng, home_cell);
        }

        for(unsigned i0=0; i0<home_cell.count-1; i0++){
            for(unsigned i1=i0+1; i1<home_cell.count; i1++){
                 base_t::calc_force(
                    base_rng,
                    true,
                    home_cell.local[i0],
                    home_cell.local[i1],
                    home_cell.local[i1].pos
                );
            }
        }

        float home_pos_f[3][4] = {{0}};
        float home_mom_f[3][4] = {{0}};
        uint32_t home_type_u[4] = {0}; // Must be zero initialise to avoid invalid indexing for un-unused bead slots

        for(unsigned i=0; i<home_cell.count; i++){
            for(int d=0; d<3; d++){
                home_pos_f[d][i] = home_cell.local[i].pos[d];
                home_mom_f[d][i] = home_cell.local[i].mom[d];
            }
            home_type_u[i] = home_cell.local[i].type;
        }

        float32x4_t home_pos[3];
        float32x4_t home_mom[3];
        float32x4_t home_force[3];
        float32x4_t home_type = vld1q_u32(home_type_u);
        for(int d=0; d<3; d++){
            home_pos[d] = vld1q_f32(home_pos_f[d]);
            home_mom[d] = vld1q_f32(home_mom_f[d]);
            home_force[d] = vdupq_n_f32(0.0f);
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

                calc_force_neon(
                    rng_state_neon,
                    home_cell.count,
                    home_pos,
                    home_mom,
                    home_force,
                    home_type,
                    other_bead, AnyExternal ? other_x : other_bead.pos
                );
            }
        }

        float home_force_f[3][4];
        for(int d=0; d<3; d++){
            vst1q_f32(home_force_f[d], home_force[d]);
            for(int i=0; i<home_cell.count; i++){
                home_cell.local[i].force[d] += home_force_f[d][i];
            }
        }
    }


};

static bool reg_SimEngineSeqNEON = SimEngineBase<SimEngineSeqNEON>::Register();

#endif // __NEON__
