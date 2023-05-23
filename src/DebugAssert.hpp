#ifndef Osprey_DebugAssert_hpp
#define Osprey_DebugAssert_hpp

#include <cassert>

#if defined(OSPREY_ENABLE_DEBUG_ASSERT)
#define DEBUG_ASSERT(msg) assert(msg)
#else
#define DEBUG_ASSERT(msg) ((void)0)
#endif

#endif