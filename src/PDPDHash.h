#ifndef PDPDHash_hpp
#define PDPDHash_hpp

#include "AbstractBead.h"

uint64_t PDPDHashGetTimeStepConstant(uint64_t seed, unsigned time_step);

double PDPDHash(uint64_t time_step_constant_, const CAbstractBead *a, const CAbstractBead *b);

#endif