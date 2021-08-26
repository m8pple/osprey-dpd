#ifndef PDPDHash_hpp
#define PDPDHash_hpp

#include "AbstractBead.h"

uint64_t PDPDHashGetTimeStepConstant(uint64_t seed, uint32_t time_step);


uint32_t PDPDHashCalcBeadHash(uint32_t bead_type, bool is_monomer, uint32_t polymer_id, uint32_t polymer_offset);

uint32_t PDPDHashBits32(uint64_t time_step_constant, uint32_t hash_a, uint32_t hash_b);

double PDPDHashU01(uint64_t time_step_constant, uint32_t hash_a, uint32_t hash_b);

double PDPDHashUSym(uint64_t time_step_constant, uint32_t hash_a, uint32_t hash_b);

#endif