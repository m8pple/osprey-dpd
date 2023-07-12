#ifndef ParallelContext_NaiveThreads_h
#define ParallelContext_NaiveThreads_h

#include "ParallelContext_Sequential.h"

#include "xxBase.h"

#include <thread>
#include <atomic>
#include "DebugAssert.h"

class ParallelContext_NaiveThreads
    : public ParallelContext_Sequential
{
private:
    mutable std::atomic<int> m_curr_workers;
    int m_max_workers = -1;

    bool OnRegister() override
    {
        DEBUG_ASSERT(m_max_workers==-1);
       
        // TODO: uniform way of controlling this
        m_max_workers=std::thread::hardware_concurrency()-1;
        m_curr_workers = 0;

        return true;
    }

    template<class T1, class T2>
    void ParTasksImpl(T1 &&t1, T2 &&t2) const
    {
        int curr=m_curr_workers.load();
        if(curr < m_max_workers){
            m_curr_workers.fetch_add(1);
            std::thread worker([&](){
                t2();
                m_curr_workers.fetch_sub(1);
            });
            t1();
            worker.join();
        }else{
            t1();
            t2();
        }
    }
public:
    virtual const char *Name() const override
    { return "NaiveThreads"; }

    virtual int Merit() const override
    { return 10; }

    virtual void ParTasks(
        std::function<void()> f1,
        std::function<void()> f2
    ) const override
    {
        ParTasksImpl(f1, f2);
    }

    virtual void ParFor(
        const range_1d &range,
        unsigned grainSize,
        std::function<void(const range_1d &)> f
    ) const override
    {
        if(range.volume() <= grainSize ){
            f(range);
        }else{
            auto parts = range.split();
            ParTasksImpl(
                [&](){ ParFor(parts.first, grainSize, f); },
                [&](){ ParFor(parts.second, grainSize, f); }
            );
        }
    }

    virtual void ParFor(
        const range_3d &range,
        unsigned grainSize,
        std::function<void(const range_3d &)> f
    ) const override
    {
        auto longest_dim_and_length=range.longest_dimension_and_length();

        if( longest_dim_and_length.second <= grainSize ){
            f(range);
        }else{
            auto parts = range.split(longest_dim_and_length.first);
            ParTasksImpl(
                [&](){ ParFor(parts.first, grainSize, f); },
                [&](){ ParFor(parts.second, grainSize, f); }
            );
        }
    }


    virtual void ParForWithSafeHalo(
        const range_3d &range,
        unsigned grainDim,
        std::function<void(const range_3d &)> f
    ) const override {
        unsigned G=std::max<int64_t>(2, grainDim);

        ///////////////////////////////////////////
        // Check range is big enough to allow parallelism.
        // For convenience we only go parallel if there is parallelism in every axis,
        // as it avoids special casing for degenenerate axis with only one block
        int parDims=0;
        for(int d=0; d<3; d++){
            parDims += range.axes[d].volume() >= 4*G; // Only have parallelism if there dimension is more than 8
        }

        // Bail if there is less than two parallelism in each dimension
        // TODO: This could be more sophisticated if we have very long skinny spaces.
        if(parDims!=3){
            f(range);
            return;
        }

        std::vector<range_1d> ranges[3];
        for(unsigned d=0; d<3; d++){
            ranges[d] = SplitForSafeHalo(range.axes[d], G);
            DEBUG_ASSERT( (ranges[d].size()%2)==0 );
        }
        range_3d coarse_dims(
            { 0, int64_t(ranges[0].size()/2) },
            { 0, int64_t(ranges[1].size()/2) },
            { 0, int64_t(ranges[2].size()/2) }
        );

        for(int sector = 0; sector < 8; sector++){
            auto fh=[&](const range_3d &fine){
                int ox=sector%2, oy=(sector/2)%2, oz=(sector/4)%2;
                fine.ForEach([&](int x, int y, int z){
                    // Select one of the 8 sub-cells for this cell at [2*x,2*x+2)x[2*y,2*y+2)x[2*z,2*z+2)
                    range_3d rr( ranges[0][x*2+ox], ranges[1][y*2+oy], ranges[2][z*2+oz] );
                    f(rr);
                });
            };

            //fh(coarse_dims);
            ParFor(coarse_dims, grainDim, fh);
        }
    }
};

#endif
