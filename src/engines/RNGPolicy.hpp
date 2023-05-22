
#ifndef RNGPolicy_hpp
#define RNGPolicy_hpp

#include <cstdint>
#include "DebugAssert.hpp"

#include "BeadIdHashRNG.hpp"

enum RNGPolicy
{
    RNGPolicy_Rng_LCG64,
    RNGPolicy_Hash_BeadTagXorShiftAdd,  // r = (a.hash_tag ^ (b.hash_tag<<16)) + (b.hash_tag ^ (a.hash_tag<<16))
    RNGPolicy_Hash_PositionXorMulSum,    // r = sum( (a.x[d] + b.x[d]) * round_mul[d], d=0..2 ) 
    RNGPolicy_Hash_BeadIdHash,           

    RNGPolicy_AsymmetricRng    = RNGPolicy_Rng_LCG64,                   // Requires use of newton's 3rd law
    RNGPolicy_HashBeadTag      = RNGPolicy_Hash_BeadTagXorShiftAdd,     // Can be used anywhere, requires hash_tag in bead
    RNGPolicy_HashPosition     = RNGPolicy_Hash_PositionXorMulSum          // Can be used anywhere, requires round_mul[3] and bit-exact position sharing
};



template<RNGPolicy TPolicy>
struct RNGImpl;

inline uint64_t SplitMix64( uint64_t z ) {
	z = z + 0x9E3779B97F4A7C15ull;
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
	return z ^ (z >> 31);
}

inline float CalcRNGScaleForU31(double stddev)
{
    // StdDev(U(-0.5,0.5)) = 0.28867513459481287 = sqrt( (b-a)**2 / 12 ), where a=-0.5, b=0.5
    const double u31_to_uHalf = 1.0 / 4294967296.0;                  // U(-2^31,2^31) -> U(-0.5,0.5)
    const double uHalf_to_uNorm = 1.0 / 0.28867513459481287;         // U(-0.5,0.5) -> Uniform with stddev=1
    const double u31_to_uNorm = u31_to_uHalf * uHalf_to_uNorm;       // U(-2^31,2^32) -> Uniform with stddev=1
    return u31_to_uNorm * stddev ;                              // U(-2^31,2^32) -> Uniform with given        
}

inline float CalcRNGScaleForSymmetricUniform(double stddev)
{
    // StdDev(U(-0.5,0.5)) = 0.28867513459481287 = sqrt( (b-a)**2 / 12 ), where a=-0.5, b=0.5
    const double uHalf_to_uNorm = 1.0 / 0.28867513459481287;         // U(-0.5,0.5) -> Uniform with stddev=1
    return uHalf_to_uNorm * stddev ;                              // U(-2^31,2^32) -> Uniform with given        
}

template<>
struct RNGImpl<RNGPolicy_Rng_LCG64>
{
    static const int RNG_POLICY = RNGPolicy_Rng_LCG64;

    static const bool USES_BEAD_ROUND_TAG = false;

    uint64_t state;
    float rng_scale;

    static const char *Name() 
    { return "LCG64"; }

    RNGImpl(double stddev, uint64_t global_seed, uint64_t round_id, uint64_t cluster_unq)
    {
        state = SplitMix64( SplitMix64(global_seed + cluster_unq) + round_id);
        rng_scale=CalcRNGScaleForU31(stddev);
    }


    template<class TBead>
    float NextUnifScaled(const TBead &, const TBead &)
    {
        uint32_t u32 = state>>32;
        state=6364136223846793005ull * state + 1;
        
        int32_t i31;
    	memcpy(&i31, &u32, 4); // Avoid undefined behaviour. Gets number in range [-2^31,2^31)
	    return i31 * rng_scale; // rng_scale = ( CCNTCell::m_invrootdt * 2^-32  )
    }

    uint32_t MakeBeadTag(uint32_t /*bead_id*/)
    {
        DEBUG_ASSERT(0); // Not expected for this rng
        return 0;
    }
};

template<>
struct RNGImpl<RNGPolicy_Hash_BeadTagXorShiftAdd>
{
    static const int RNG_POLICY = RNGPolicy_Hash_BeadTagXorShiftAdd;

    static const bool USES_BEAD_ROUND_TAG = true;

    float rng_scale;
    uint64_t state;

    static const char *Name() 
    { return "BeadTagXorShiftAdd"; }

    RNGImpl(double stddev, uint64_t global_seed, uint64_t round_id, uint64_t cluster_unq)
    {
        state = SplitMix64( SplitMix64(global_seed + cluster_unq) + round_id);
        rng_scale = CalcRNGScaleForU31(stddev);
    }


    template<class TBead>
    float NextUnifScaled(const TBead &a, const TBead &b)
    {
        uint32_t u32 = (a.round_tag ^ (b.round_tag<<16) ) + (b.round_tag ^ (a.round_tag<<16) );
        
        int32_t i31;
    	memcpy(&i31, &u32, 4); // Avoid undefined behaviour. Gets number in range [-2^31,2^31)
	    return i31 * rng_scale; // rng_scale = ( CCNTCell::m_invrootdt * 2^-32  )
    }

    uint32_t MakeBeadTag(uint32_t /*bead_id*/)
    {
        uint32_t res=state>>32;
        state=6364136223846793005ull * state + 1;
        return res;
    }

};

template<>
struct RNGImpl<RNGPolicy_Hash_PositionXorMulSum>
{
    static const int RNG_POLICY = RNGPolicy_Hash_PositionXorMulSum;

    static const bool USES_BEAD_ROUND_TAG = false;

    static const char *Name() 
    { return "BeadTagXorMulSum"; }

    uint32_t dim_scales[3];
    float rng_scale;

    uint32_t to_bits(const float &x)
    {
        uint32_t tmp;
        memcpy(&tmp, &x, 4);
        return tmp;
    }

    unsigned weight(uint32_t x)
    {
        return __builtin_popcount(x);
    }

    RNGImpl(double stddev, uint64_t global_seed, uint64_t round_id, uint64_t cluster_unq)
    {
        rng_scale = CalcRNGScaleForU31(stddev);
        uint64_t lcg=SplitMix64( SplitMix64(global_seed + cluster_unq) + round_id);
        for(int d=0; d<3; d++){
    retry:
            lcg = 6364136223846793005ull * lcg + 1;
            dim_scales[d] = (lcg >> 32) | 1;
            uint32_t w=weight(dim_scales[d]);
            if( (w < 8) || (24 < 8) ){
                goto retry;
            }
        }
    }

    template<class TBead>
    float NextUnifScaled(const TBead &a, const TBead &b)
    {
        uint32_t u32 = 0;
        for(unsigned d=0; d<3; d++){
            u32 += dim_scales[d] * ( to_bits(a.pos[d]) + to_bits(b.pos[d]) );
        }
        
        int32_t i31;
    	memcpy(&i31, &u32, 4); // Avoid undefined behaviour. Gets number in range [-2^31,2^31)
	    return i31 * rng_scale; // rng_scale = ( CCNTCell::m_invrootdt * 2^-32  )
    }

    uint32_t MakeBeadTag(uint32_t /*bead_id*/)
    {
        DEBUG_ASSERT(0); // Not expected for this rng
        return 0;
    }
};


template<>
struct RNGImpl<RNGPolicy_Hash_BeadIdHash>
{
    static const int RNG_POLICY = RNGPolicy_Hash_BeadIdHash;

    static const bool USES_BEAD_ROUND_TAG = false;

    uint64_t round_base;
    float rng_scale;

    static const char *Name() 
    { return "BeadIdHash"; }

    RNGImpl(double stddev, uint64_t global_seed, uint64_t round_id, uint64_t cluster_unq)
    {
        rng_scale = CalcRNGScaleForU31(stddev);
        round_base = bead_id_hash_rng__round_hash(global_seed, round_id);
    }

    template<class TBead>
    float NextUnifScaled(const TBead &a, const TBead &b)
    {
        return rng_scale * bead_id_hash_rng__random_symmetric_uniform(round_base, a.bead_id, b.bead_id);
    }
};





#endif
