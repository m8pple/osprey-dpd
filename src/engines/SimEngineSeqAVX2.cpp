#include "SimEngineSeq.h"

#ifdef __AVX2__

#include "immintrin.h"
#include "pmmintrin.h"

struct  SimEngineSeqAVX2Policy
    : EnginePolicyConcept
{
};

class SimEngineSeqAVX2
    : public SimEngineSeq<SimEngineSeqAVX2Policy>
{
public:
    using IIntegrationEngine::support_result;
    using IIntegrationEngine::run_result;

    std::string Name() const override
    {
        return "SimEngineSeqAVX2";
    }

private:

    __m256i rng_state_avx2;
    float rng_scale;

    support_result on_import(ISimBox *box) override 
    {
        uint64_t rng_state = SplitMix64(SplitMix64( box->GetRNGSeed() ) + box->GetCurrentTime());

        rng_state_avx2 = _mm256_set_epi64x(
            rng_state=6364136223846793005ull * rng_state + 1,
            rng_state=6364136223846793005ull * rng_state + 1,
            rng_state=6364136223846793005ull * rng_state + 1,
            rng_state=6364136223846793005ull * rng_state + 1
          );

        rng_scale = ldexpf(rng_stddev,-32);

        return {Supported};
    }

    __m256 RandUnifScaledAVX2(__m256i &state)
    {
        /*
        x ^= x << 13;
	    x ^= x >> 7;
	    x ^= x << 17;
        */

       __m256i orig=state;

        state = state ^ _mm256_slli_epi64(state, 13);
        state = state ^ _mm256_srli_epi64(state, 7);
        state = state ^ _mm256_slli_epi64(state, 17);

        __m256i u = state; // orig + _mm256_shuffle_epi32( state, (2<<6)|(3<<4)|(0<<2)|(1<<0));

        __m256 uf = _mm256_cvtepi32_ps(u);
        return uf * _mm256_set1_ps(rng_scale);
    }

    bool any_nan(__m256 x)
    {
        float t[8];
        _mm256_storeu_ps(t, x);
        for(unsigned i=0; i<8; i++){
            if(isnanf(t[i])){
                return true;
            }
        }
        return false;
    }

    void calc_force_avx2(
        __m256i &rng_state,

        unsigned home_count,
        const __m256 home_pos[3],
        const __m256 home_mom[3],
        __m256 home_force[3],
        __m256i home_type,

        Bead &other,
        float other_x[4]
    ) {
        DEBUG_ASSERT(home_count <= 8);

        const __m256 ONE = _mm256_set1_ps(1.0f);

        __m256 dx[3], dx2[3];

        for(int d=0; d<3; d++){
            dx[d] = home_pos[d] - _mm256_set1_ps(other_x[d]);
            dx2[d] = dx[d] * dx[d];
        }
	
		__m256 dr2 = dx2[0] + dx2[1] + dx2[2];

        static const float min_r = 0.000000001f;
        static const float min_r2 = min_r * min_r;

        __m256 active =  _mm256_and_ps(  _mm256_cmp_ps(dr2, _mm256_set1_ps(min_r2), _CMP_GT_OQ), _mm256_cmp_ps( dr2, ONE, _CMP_LT_OQ  ) );

        if(_mm256_testz_si256(_mm256_castps_si256(active),_mm256_castps_si256(active))){
            return;
        }

        __m256 dv[3], dx_dv[3];

        for(int d=0; d<3; d++){
            dv[d] = home_mom[d] - _mm256_set1_ps(other.mom[d]);
            dx_dv[d] = dx[d] * dv[d];
        }

        __m256 dr = _mm256_sqrt_ps(dr2);
        __m256 wr = ONE - dr;
        __m256 wr2 = wr*wr;
        __m256 inv_dr = _mm256_rsqrt_ps(dr2); // = ONE / dr;

        const float *coeff_row = (float*)&coefficients[other.type * num_bead_types]; // Note we are type-punning from std::pair<float> to float
        
        __m256 conservative_coeff = _mm256_i32gather_ps(coeff_row, home_type, 8);    // Coeffs appear at 64-bit offsets
        __m256 dissipative_coeff = _mm256_i32gather_ps(coeff_row + 1, home_type, 8); // Then they are offset by 1 float into the array

        __m256 conForce = conservative_coeff * wr;

        __m256 rdotv = (dx_dv[0] + dx_dv[1] + dx_dv[2]) * inv_dr;
        __m256 gammap	= dissipative_coeff * wr2;

        __m256 dissForce	= -gammap*rdotv;
        __m256 randForce	= _mm256_sqrt_ps(gammap) * RandUnifScaledAVX2(rng_state);
        __m256 normTotalForce = (conForce + dissForce + randForce) * inv_dr;

        normTotalForce = _mm256_and_ps( normTotalForce, active );

        __m256 newForce[3];
        for(int d=0; d<3; d++){
            newForce[d] = normTotalForce * dx[d];
            DEBUG_ASSERT(!any_nan(newForce[d]));

            home_force[d] += newForce[d];

            if(1){
                // Hopefully these get interleaved with other instructions
                newForce[d] = _mm256_hadd_ps(newForce[d], newForce[d]);
                newForce[d] = _mm256_hadd_ps(newForce[d], newForce[d]);

                other.force[d] -= _mm256_cvtss_f32(newForce[d]);

                // A lot of the time this is wasted, as there are often 4 or fewer beads in
                // cell. Could be made conditional?
                if(home_count > 4){
                    other.force[d] -= _mm_cvtss_f32( _mm256_extractf128_ps( newForce[d], 1 ) );
                }
            }else{
                float tmp[8];
                _mm256_storeu_ps(tmp, newForce[d]);
                for(unsigned i=0; i<3; i++){
                    other.force[i] -= tmp[i];
                }
            }

            DEBUG_ASSERT(std::abs(other.force[d]) < 10000);
        }
    }

#if 0
    void calc_forces_avx2_dual(
        __m256i &rng_state,

        unsigned home_count,
        const __m256 home_pos[3],
        const __m256 home_mom[3],
        __m256 home_force[3],
        __m256i home_type,

        Bead &other_a,
        float other_x_a[4],
        Bead &other_b,
        float other_x_b[4]
    ) {
        DEBUG_ASSERT(home_count <= 8);

        const __m256 ONE = _mm256_set1_ps(1.0f);

        __m256 dx[3], dx2[3];

        for(int d=0; d<3; d++){
            __m256 other_x = _mm256_set_m128(_mm_set1_ps(other_x_b[d]),_mm_set1_ps(other_x_a[d]));
            dx[d] = home_pos[d] - other_x;
            dx2[d] = dx[d] * dx[d];
        }
	
		__m256 dr2 = dx2[0] + dx2[1] + dx2[2];

        static const float min_r = 0.000000001f;
        static const float min_r2 = min_r * min_r;

        __m256 active =  _mm256_and_ps(  _mm256_cmp_ps(dr2, _mm256_set1_ps(min_r2), _CMP_GT_OQ), _mm256_cmp_ps( dr2, ONE, _CMP_LT_OQ  ) );

        if(_mm256_testz_si256(_mm256_castps_si256(active),_mm256_castps_si256(active))){
            return;
        }

        __m256 dv[3], dx_dv[3];

        for(int d=0; d<3; d++){
            __m256 other_mom = _mm256_set_m128(_mm_set1_ps(other_b.mom[d]),_mm_set1_ps(other_a.mom[d]));
            dv[d] = home_mom[d] - other_mom;
            dx_dv[d] = dx[d] * dv[d];
        }

        __m256 dr = _mm256_sqrt_ps(dr2);
        __m256 wr = ONE - dr;
        __m256 wr2 = wr*wr;
        __m256 inv_dr = _mm256_rsqrt_ps(dr2); // ONE / dr;

        __m256i coeff_offset = _mm256_set_m128i(_mm_set1_epi32(other_b.type*num_bead_types),_mm_set1_epi32(other_a.type*num_bead_types));
        coeff_offset = _mm256_add_epi32(coeff_offset, home_type);

        const float *coeff_root = (float*)&coefficients[0]; // Note we are type-punning from std::pair<float> to float
        
        __m256 conservative_coeff = _mm256_i32gather_ps(coeff_root, coeff_offset, 8);    // Coeffs appear at 64-bit offsets
        __m256 dissipative_coeff = _mm256_i32gather_ps(coeff_root + 1, coeff_offset, 8); // Then they are offset by 1 float into the array

        __m256 conForce = conservative_coeff * wr;

        __m256 rdotv = (dx_dv[0] + dx_dv[1] + dx_dv[2]) * inv_dr;
        __m256 gammap	= dissipative_coeff * wr2;

        __m256 dissForce	= -gammap*rdotv;
        __m256 randForce	= _mm256_sqrt_ps(gammap) * RandUnifScaledAVX2(rng_state);
        __m256 normTotalForce = (conForce + dissForce + randForce) * inv_dr;

        normTotalForce = _mm256_and_ps( normTotalForce, active );

        __m256 newForce[3];
        for(int d=0; d<3; d++){
            newForce[d] = normTotalForce * dx[d];
            DEBUG_ASSERT(!any_nan(newForce[d]));

            home_force[d] += newForce[d];

            newForce[d] = _mm256_hadd_ps(newForce[d], newForce[d]);
            newForce[d] = _mm256_hadd_ps(newForce[d], newForce[d]);

            other_a.force[d] -= _mm256_cvtss_f32(newForce[d]);
            other_b.force[d] -= _mm_cvtss_f32( _mm256_extractf128_ps( newForce[d], 1 ) );

            DEBUG_ASSERT(std::abs(other_a.force[d]) < 10000);
        }
    }
    #endif

    template<bool AnyExternal>
    void update_forces(rng_t &rng, Cell &home_cell)
    {
        home_cell.to_move = home_cell.count;

        if(home_cell.count==0){
            return;
        }

        if(home_cell.count > 8){
            return SimEngineSeq::update_forces<AnyExternal>(rng, home_cell);
        }

        for(unsigned i0=0; i0<home_cell.count-1; i0++){
            for(unsigned i1=i0+1; i1<home_cell.count; i1++){
                calc_force(rng, /*reflect_force*/ true, home_cell.local[i0], home_cell.local[i1], home_cell.local[i1].pos );
            }
        }

        float home_pos_f[3][8] = {{0}};
        float home_mom_f[3][8] = {{0}};
        uint32_t home_type_u[8] = {0}; // Must be zero initialise to avoid invalid indexing for un-unused bead slots

        for(unsigned i=0; i<home_cell.count; i++){
            for(int d=0; d<3; d++){
                home_pos_f[d][i] = home_cell.local[i].pos[d];
                home_mom_f[d][i] = home_cell.local[i].mom[d];
            }
            home_type_u[i] = home_cell.local[i].type;
        }

        __m256 home_pos[3];
        __m256 home_mom[3];
        __m256 home_force[3];
        __m256i home_type = _mm256_loadu_si256((__m256i*)home_type_u);
        for(int d=0; d<3; d++){
            home_pos[d] = _mm256_loadu_ps(home_pos_f[d]);
            home_mom[d] = _mm256_loadu_ps(home_mom_f[d]);
            home_force[d] = _mm256_setzero_ps();
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

                calc_force_avx2(
                    rng_state_avx2,
                    home_cell.count,
                    home_pos,
                    home_mom,
                    home_force,
                    home_type,
                    other_bead, AnyExternal ? other_x : other_bead.pos
                );
            }
        }

        float home_force_f[3][8];
        for(int d=0; d<3; d++){
            _mm256_storeu_ps(home_force_f[d], home_force[d]);
            for(int i=0; i<home_cell.count; i++){
                home_cell.local[i].force[d] += home_force_f[d][i];
            }
        }
    }

#if 0
    template<bool AnyExternal>
    void update_forces_dual(Cell &home_cell)
    {
        home_cell.to_move = home_cell.count;

        DEBUG_ASSERT(home_cell.count<=4);

        if(home_cell.count==0){
            return;
        }

        for(unsigned i0=0; i0<home_cell.count-1; i0++){
            for(unsigned i1=i0+1; i1<home_cell.count; i1++){
                calc_forces( home_cell.local[i0], home_cell.local[i1], home_cell.local[i1].pos );
            }
        }

        float home_pos_f[3][8] = {{0}};
        float home_mom_f[3][8] = {{0}};
        uint32_t home_type_u[8] = {0}; // Must be zero initialise to avoid invalid indexing for un-unused bead slots

        for(unsigned i=0; i<home_cell.count; i++){
            for(int d=0; d<3; d++){
                home_pos_f[d][i] = home_cell.local[i].pos[d];
                home_mom_f[d][i] = home_cell.local[i].mom[d];
                home_pos_f[d][i+4] = home_cell.local[i].pos[d];
                home_mom_f[d][i+4] = home_cell.local[i].mom[d];
            }
            home_type_u[i] = home_cell.local[i].type;
            home_type_u[i+4] = home_cell.local[i].type;
        }

        __m256 home_pos[3];
        __m256 home_mom[3];
        __m256 home_force[3];
        __m256i home_type = _mm256_loadu_si256((__m256i*)home_type_u);
        for(int d=0; d<3; d++){
            home_pos[d] = _mm256_loadu_ps(home_pos_f[d]);
            home_mom[d] = _mm256_loadu_ps(home_mom_f[d]);
            home_force[d] = _mm256_setzero_ps();
        }

        float ready_x_a[3];
        Bead *ready_bead;

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

            unsigned other_i=0;
            while(other_i < other_cell.count){
                Bead &other_bead = other_cell.local[other_i];
                float other_x[4];
                if(AnyExternal){
                    for(int d=0; d<CALC_DIM; d++){
                        other_x[d] = other_bead.pos[d] + other_delta[d];
                    }
                }

                if(other_i+1 == other_cell.count){
                    static const __mm256 invalid_pos[3] = {_mm256_set1_ps(-1000),_mm256_set1_ps(-1000),_mm256_set1_ps(-1000)};

                    calc_forces_avx2(
                        rng_state_avx2,
                        home_cell.count,
                        home_pos,
                        home_mom,
                        home_force,
                        home_type,
                        other_bead, AnyExternal ? other_x : other_bead.pos,
                        other_bead, invalid_pos
                    );
                }else{
                    Bead &other_bead_b = other_cell.local[other_i];
                    float other_b_x[4];
                    if(AnyExternal){
                        for(int d=0; d<CALC_DIM; d++){
                            other_b_x[d] = other_bead_b.pos[d] + other_delta[d];
                        }
                    }

                    calc_forces_avx2_dual(
                        rng_state_avx2,
                        home_cell.count,
                        home_pos,
                        home_mom,
                        home_force,
                        home_type,
                        other_bead, AnyExternal ? other_x : other_bead.pos,
                        other_bead_b, AnyExternal ? other_b_x : other_bead_b.pos
                    );
                }
            }
        }

        float home_force_f[3][8];
        for(int d=0; d<3; d++){
            _mm256_storeu_ps(home_force_f[d], home_force[d]);
            for(int i=0; i<home_cell.count; i++){
                home_cell.local[i].force[d] += home_force_f[d][i] + home_force_f[d][i+4];
            }
        }
    }
    #endif


    __attribute__((noinline)) void update_forces()
    {
        rng_t rng(rng_stddev, global_seed, round_id, 0);

        // TODO : branching based on any external may introduce data-dependent
        // control and bloat instruction size. Is minor saving worth it?
        for(Cell &cell : cells){
            if(cell.any_external){
                update_forces<true>(rng, cell);
            }else{
                update_forces<false>(rng, cell);
            }
        }
    }

};

static bool reg_SimEngineSeqAVX2 = SimEngineBase<SimEngineSeqAVX2>::Register();

#endif