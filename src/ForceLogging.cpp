#include "ForceLogging.h"

#include <iostream>

ForceLogging *ForceLogging::logger=0;

void ForceLogging::SetTime(long t)
{
    m_t=t;
}

void ForceLogging::LogBeadProperty(long bead_id, const char *name, int dims, const double *x)
{
    std::cout<<"Prop,"<<m_t<<","<<bead_id<<",,,"<<name;
    for(int i=0; i<dims; i++){
        std::cout<<","<<x[i];
    }
    for(int i=dims;i<3;i++){
        std::cout<<",";
    }
    std::cout<<"\n";
}

void ForceLogging::LogBeadPairProperty(long bead_id0,long bead_id1, const char *name, ForceLoggingFlags flags, int dims, const double *x)
{
    bool flip=flags!=Asymmetric && (bead_id0 > bead_id1);
    if(flip){
        std::swap(bead_id0, bead_id1);
    }
    std::cout<<"Prop,"<<m_t<<","<<bead_id0<<","<<bead_id1<<",,"<<name;
    for(int i=0; i<dims; i++){
        if(flip && (flags==SymmetricFlipped)){
            std::cout<<","<<-x[i];
        }else{
            std::cout<<","<<x[i];
        }
    }
    for(int i=dims;i<3;i++){
        std::cout<<",";
    }
    std::cout<<"\n";
}

void ForceLogging::LogBeadTripleProperty(long bead_id0, long bead_id1, long bead_id2, const char *name, int dims, const double *x)
{
    std::cout<<"Prop,"<<m_t<<","<<bead_id0<<","<<bead_id1<<","<<bead_id2<<","<<name;
    for(int i=0; i<dims; i++){
        std::cout<<","<<x[i];
    }
    for(int i=dims;i<3;i++){
        std::cout<<",";
    }
    std::cout<<"\n";
}
