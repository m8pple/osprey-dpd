#ifndef ParallelContext_TBB_h
#define ParallelContext_TBB_h

#include "ParallelContext_Sequential.h"

#include "xxBase.h"

#include <thread>
#include "DebugAssert.h"

#include "tbb/parallel_invoke.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range3d.h"

class ParallelContext_TBB
    : public ParallelContext_Sequential
{
public:
    virtual const char *Name() const override
    { return "TBB"; }

    virtual int Merit() const override
    { return 100; }

    virtual void ParTasks(
        std::function<void()> f1,
        std::function<void()> f2
    ) const override
    {
        tbb::parallel_invoke(
            f1,
            f2
        );
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
            tbb::parallel_for(
                tbb::blocked_range<int64_t>(range.lower, range.upper, grainSize),
                [&]( const tbb::blocked_range<int64_t> &fine ){
                    f( {fine.begin(), fine.end()} );
                }
            );
        }
    }

    virtual void ParFor(
        const range_3d &range,
        unsigned grainSize,
        std::function<void(const range_3d &)> f
    ) const override
    {
        if(range.volume() <= grainSize ){
            f(range);
        }else{
            tbb::blocked_range3d<int64_t> grid(
                range.axes[0].lower, range.axes[0].upper, grainSize,
                range.axes[1].lower, range.axes[1].upper, grainSize,
                range.axes[2].lower, range.axes[2].upper, grainSize
            );
            tbb::parallel_for(
                grid,
                [&]( const tbb::blocked_range3d<int64_t> &fine ){
                    range_3d block{
                        { fine.pages().begin(), fine.pages().end() },
                        { fine.rows().begin(), fine.rows().end() },
                        { fine.cols().begin(), fine.cols().end() }
                    };
                    
                    f( block );
                }
            );
        }
    }


    virtual void ParForWithSafeHalo(
        const range_3d &range,
        unsigned grainSize,
        std::function<void(const range_3d &)> f
    ) const override {
        int64_t G=std::max<int64_t>(2, grainSize);

        ///////////////////////////////////////////
        // Check range is big enough to allow parallelism.
        // For convenience we only go parallel if there is parallelism in every axis,
        // as it avoids special casing for degenenerate axis with only one block
        int parDims=0;
        for(int d=0; d<3; d++){
            parDims += range.axes[d].volume() >= 4*G; // Only have parallelism if there dimension is at least 4 blocks
        }

        // Bail if there is less than two parallelism in each dimension
        // TODO: This could be more sophisticated if we have very long skinny spaces or flat spaces
        if(parDims!=3){
            f(range);
            return;
        }

        std::vector<range_1d> ranges[3];
        for(unsigned d=0; d<3; d++){
            ranges[d] = SplitForSafeHalo(range.axes[d], grainSize);
            DEBUG_ASSERT( (ranges[d].size()%2)==0 );
        }
        tbb::blocked_range3d<int64_t> coarse(
            0, int64_t(ranges[0].size()/2), 1,
            0, int64_t(ranges[1].size()/2), 1,
            0, int64_t(ranges[2].size()/2), 1
        );

        for(int sector = 0; sector < 8; sector++){
            auto fh=[&](const tbb::blocked_range3d<int64_t> &fine){
                int ox=sector%2, oy=(sector/2)%2, oz=(sector/4)%2;

                range_3d local{{},{},{}};                
                for(int64_t iz=fine.cols().begin(); iz < fine.cols().end(); iz++){
                    local.axes[2]=ranges[2][iz*2+oz];
                    for(int64_t iy=fine.rows().begin(); iy < fine.rows().end(); iy++){
                        local.axes[1]=ranges[1][iy*2+oy];
                        for(int64_t ix=fine.pages().begin(); ix < fine.pages().end(); ix++){
                            local.axes[0]=ranges[0][ix*2+ox];
                            f(local);
                        }                    
                    }
                }
            };

            tbb::parallel_for(coarse, fh, tbb::simple_partitioner());
        }
    }
};

#endif
