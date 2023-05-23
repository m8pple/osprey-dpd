#ifndef Osprey_DebugAssert_hpp
#define Osprey_DebugAssert_hpp

#if defined(OSPREY_ENABLE_DEBUG_ASSERT)
#define DEBUG_ASSERT(msg) DEBUG_ASSERT(msg)
#else
#define DEBUG_ASSERT(msg) ((void)0)
#endif

#endif