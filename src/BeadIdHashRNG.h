#ifndef bead_id_hash_rng_hpp
#define bead_id_hash_rng_hpp

#include <cstdint>
#include <memory>

inline uint64_t bead_id_hash_rng__SplitMix64(uint64_t z)
{
	z = z + 0x9E3779B97F4A7C15ull;
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
	return z ^ (z >> 31);
}

inline uintptr_t bead_id_hash_rng__round_hash(uint64_t global_seed, uint64_t round_id)
{
    uint64_t a = bead_id_hash_rng__SplitMix64(global_seed);
    return bead_id_hash_rng__SplitMix64(a + round_id);
}

inline uint32_t bead_id_hash_rng__bead_hash(uint64_t round_hash, uint32_t ida, uint32_t idb)
{
    // This is a completely arbitrary constant. It hasn't been optimised
    // at all. I don't really recommend people use this in production,
    // it is mainly to support debugging.
    /*const uint32_t C1 = 1459366445; // e/8
    const uint32_t C2 = 1686629713; // pi/8

    uint32_t a = std::min(ida,idb) + uint32_t(round_hash);
    uint32_t b = std::max(ida,idb) + uint32_t(round_hash>>32);

    uint32_t c = (a + (b>>16)) & 0x7FFFFFFFul;
    uint32_t d = (b + (a>>16)) & 0x7FFFFFFFul;

    c *= C1;
    d *= C2;

    return c ^ d;*/


    uint64_t ab=(uint64_t(std::max(ida,idb))<<32) | std::min(ida,idb);
    uint64_t xx=bead_id_hash_rng__SplitMix64(round_hash+ab);
    return uint32_t(xx) ^ uint32_t(xx>>32);
}

// Return random float in [-0.5,+0.5)
inline float bead_id_hash_rng__random_symmetric_uniform(uintptr_t round_hash, uint32_t ida, uint32_t idb)
{
    uint32_t u=bead_id_hash_rng__bead_hash(round_hash, ida, idb);
    int32_t s;
    memcpy(&s, &u, 4);
    return s * 0.000000000232831f; // s * 2^-32:  [-2^31...,2^31)] -> [-0.5,0.5)
}


#endif
