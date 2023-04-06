
#ifndef RNGPolicy_hpp
#define RNGPolicy_hpp

#include <cstdint>

enum RNGPolicy
{
    RNGPolicy_Rng_LCG64,
    RNGPolicy_Hash_RoundTagAddShiftAdd,  // r = a.round_tag + b.round_tag; r = r + (r<<16)

    RNGPolicy_AsymmetricRng    = RNGPolicy_Rng_LCG64,                // Requires use of newton's 3rd law
    RNGPolicy_HashBeadTag      = RNGPolicy_Hash_RoundTagAddShiftAdd,          // Can be used anywhere, requires round_tag in bead

};

template<RNGPolicy TPolicy>
struct RNGImpl;

template<>
struct RNGImpl<RNGPolicy_Rng_LCG64>
{
    static const int RNG_POLICY = RNGPolicy_Rng_LCG64;

    uint64_t &rng_state;
    float rng_scale;

    template<class TContext>
    RNGImpl(TContext &ctxt)
        : rng_state(ctxt.rng_state)
        , rng_scale(ctxt.rng_scale)
    {
    }

    template<class TBead>
    float NextUnifScaled(const TBead &, const TBead &)
    {
        uint32_t u32 = rng_state>>32;
        rng_state=6364136223846793005ull * rng_state + 1;
        
        int32_t i31;
    	memcpy(&i31, &u32, 4); // Avoid undefined behaviour. Gets number in range [-2^31,2^31)
	    return i31 * rng_scale; // rng_scale = ( CCNTCell::m_invrootdt * 2^-32  )
    }
};

template<>
struct RNGImpl<RNGPolicy_Hash_RoundTagAddShiftAdd>
{
    static const int RNG_POLICY = RNGPolicy_Hash_RoundTagAddShiftAdd;

    float rng_scale;

    template<class TContext>
    RNGImpl(TContext &ctxt)
        : rng_scale(ctxt.rng_scale)
    {}

    template<class TBead>
    float NextUnifScaled(const TBead &a, const TBead &b)
    {
        uint32_t u32 = (a.round_tag + b.round_tag ) * 19937;
        
        int32_t i31;
    	memcpy(&i31, &u32, 4); // Avoid undefined behaviour. Gets number in range [-2^31,2^31)
	    return i31 * rng_scale; // rng_scale = ( CCNTCell::m_invrootdt * 2^-32  )
    }
};
/*
template<>
struct RNGImpl<RNGPolicy_Hash_PositionHash>
{
    static const int RNG_POLICY = RNGPolicy_Hash_PositionHash;

    float rng_scale;
    const uint32_t *round_position_hash_biases;

    template<class TContext>
    RNGImpl(TContext &ctxt)
        : rng_scale(ctxt.rng_scale)
        , round_position_hash_biases(ctxt.round_position_hash_biases)
    {}

    template<class TBead>
    float NextUnifScaled(const TBead &a, const TBead &b)
    {
        uint32_t u32 = (a.round_tag + b.round_tag ) * 19937;
        
        int32_t i31;
    	memcpy(&i31, &u32, 4); // Avoid undefined behaviour. Gets number in range [-2^31,2^31)
	    return i31 * rng_scale; // rng_scale = ( CCNTCell::m_invrootdt * 2^-32  )
    }
};
*/

#endif
