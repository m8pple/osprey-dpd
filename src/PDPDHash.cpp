#include "PDPDHash.h"

#include <cassert>
#include <cmath>

static uint64_t riscv_mix64_m2(uint64_t x)
{
    // From custom hash prospector. Not as good as splitmix64 in quality, but much faster/smaller on riscv
    const uint32_t C=0xd392d2a7; // 2
    x = x ^ (x>>32);     // 1
    x = x * C;           // 4
    x = x ^ (x>>32);     // 1
    x = x * C;           // 4
    x = x ^ (x>>32);     // 1
    return x;
}

uint32_t PDPDHashCalcBeadHash(uint32_t bead_type, bool is_monomer, uint32_t polymer_id, uint32_t polymer_offset)
{
    assert( is_monomer ? (polymer_offset==0 && polymer_id<(1u<<27))
                                    : (polymer_offset<128 && polymer_id<(1u<<20)));
    uint32_t base=(uint32_t(polymer_offset)<<20) | polymer_id;
    base |= uint32_t(is_monomer)<<27;
    base |= bead_type << 28;
    return base;
}

uint64_t PDPDHashGetTimeStepConstant(uint64_t seed, uint32_t time_step)
{
    uint64_t base=seed + riscv_mix64_m2(time_step);
    return riscv_mix64_m2(base);
}

uint32_t PDPDHashBits32(uint64_t time_step_constant, uint32_t hash_a, uint32_t hash_b)
{
    if(hash_a>hash_b){
        std::swap(hash_a,hash_b);
    }
    // m3 is safer from randomness testing perspective, but m2 is ok in practise for many beads.
    return riscv_mix64_m2( (hash_a|(uint64_t(hash_b)<<32)) ^ time_step_constant);
}

double PDPDHashU01(uint64_t time_step_constant, uint32_t hash_a, uint32_t hash_b)
{
    return PDPDHashBits32(time_step_constant, hash_a, hash_b) * 0.00000000023283064365386962890625;
}

static int32_t uint32_to_int32(uint32_t x)
{
    // Sigh, Avoid undefined behaviour. Compiler should optimise it out.
    // https://stackoverflow.com/a/13208789
    int32_t res;
    if (x <= INT32_MAX) {
        res=static_cast<uint32_t>(x);
    }else{
        assert(x >= (uint32_t)INT32_MIN);
        res= static_cast<uint32_t>(x - INT32_MIN) + INT32_MIN;
    }
    assert(x == uint32_t(res) );  // int32_t -> uint32_t is well-defined
    return res;
}

double PDPDHashUSym(uint64_t time_step_constant, uint32_t hash_a, uint32_t hash_b)
{
    uint32_t bits= PDPDHashBits32(time_step_constant, hash_a, hash_b);
    int32_t ibits=uint32_to_int32(bits);

    const double scale=std::ldexp(1.0, -32); // Gives range of [-0.5,0.5]  (same as Osprey-DPD)
    return ibits * scale; 
}