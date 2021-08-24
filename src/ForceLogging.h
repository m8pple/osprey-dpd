#ifndef ForceLogging_h
#define ForceLogging_h

#include "AbstractBead.h"

enum ForceLoggingFlags
{
    Asymmetric,       // need to record property seperately for (b0,b1) and (b1,b0)
    Symmetric,        // Same value for (b0,b1) and (b1,b0)
    SymmetricFlipped  // Value for (b0,b1) is negative of value for (b1,b0)
};

class ForceLogging
{
public:
    void SetTime(long t);
    void LogBeadProperty(long bead_id, const char *name, int dims, const double *x);
    void LogBeadPairProperty(long bead_id0,long bead_id1, const char *name, ForceLoggingFlags flags, int dims, const double *x);
    void LogBeadTripleProperty(long bead_id0, long bead_id1, long bead_id2, const char *name, int dims, const double *x);

    static ForceLogging *logger;
private:
    long m_t;
};

inline void LogBeadProperty(const CAbstractBead *b, const char *name, int dims, const double *x)
{
    if(ForceLogging::logger){
        ForceLogging::logger->LogBeadProperty(b->GetId(),name,dims,x);
    }
}

inline void LogBeadPairProperty(const CAbstractBead *b0,const CAbstractBead *b1, const char *name, ForceLoggingFlags flags, int dims, const double *x)
{
    if(ForceLogging::logger){
        ForceLogging::logger->LogBeadPairProperty(b0->GetId(),b1->GetId(),name,flags,dims,x);
    }
}

inline void LogBeadPairProperty(const CAbstractBead *b0,const CAbstractBead *b1, const char *name, ForceLoggingFlags flags, double x)
{
    if(ForceLogging::logger){
        ForceLogging::logger->LogBeadPairProperty(b0->GetId(),b1->GetId(),name,flags,1,&x);
    }
}

inline void LogBeadTripleProperty(const CAbstractBead *b0,const CAbstractBead *b1,const CAbstractBead *b2, const char *name, int dims, const double *x)
{
    if(ForceLogging::logger){
        ForceLogging::logger->LogBeadTripleProperty(b0->GetId(), b1->GetId(), b2->GetId(),name,dims,x);
    }
}

#endif