#ifndef ParallelContext_Sequential_h
#define ParallelContext_Sequential_h

#include "ParallelContext.h"

#include "xxBase.h"

/*
This is a stateless default implementation.

Because it is stateless, other implementations are explicitly allowed
to derive from this one.
*/
class ParallelContext_Sequential
    : public ParallelContext
{
protected:
    struct Range1DContextImpl
        : public Range1DContext
    {
        Range1DContextImpl(const range_1d &_r, unsigned _grainSize)
            : range(_r)
            , grainSize(_grainSize)
        {}

        range_1d range;
        unsigned grainSize;
    };

    struct Range3DContextImpl
        : public Range3DContext
    {
        Range3DContextImpl(const range_3d &_r, unsigned _grainVolume)
            : range(_r)
            , grainVolume(_grainVolume)
        {}

        range_3d range;
        unsigned grainVolume;
    };

public:
    virtual const char *Name() const override
    { return "Sequential"; }

    virtual int Merit() const override
    { return 1; }

    virtual Range1DToken CreateParForToken(const range_1d& range, unsigned grainSize) const override
    {
        return Range1DToken( new Range1DContextImpl(range,grainSize) );
    }
    
    virtual Range3DToken CreateParForToken(const range_3d& range, unsigned grainVolume) const override
    {
        return Range3DToken( new Range3DContextImpl(range, grainVolume ) );
    }

    virtual Range3DToken CreateParForWithSafeHaloToken(const range_3d& range, unsigned grainVolume) const override
    {
        return Range3DToken( new Range3DContextImpl(range, grainVolume) );
    }

    virtual void ParTasks(
        std::function<void()> f1,
        std::function<void()> f2
    ) const override
    {
        f1();
        f2();
    }

    virtual void ParFor(
        const range_1d &range,
        unsigned grainSize,
        std::function<void(const range_1d &)> f
    ) const override
    {
        f(range);
    }

    virtual void ParFor(
        const Range1DToken &ctxt,
        std::function<void(const range_1d &)> f
    ) const override {
        auto ic=dynamic_cast<const Range1DContextImpl *>(ctxt.get());
        if(!ic){
            xxBase::FatalTraceGlobal("ParFor called with incorrect context.");
        }
        ParFor(ic->range, ic->grainSize, f);
    }

    virtual void ParFor(
        const range_3d &range,
        unsigned grainVolume,
        std::function<void(const range_3d &)> f
    ) const override {
        f(range);
    }

    virtual void ParFor(
        const Range3DToken &ctxt,
        std::function<void(const range_3d &)> f
    ) const override {
        auto ic=dynamic_cast<const Range3DContextImpl *>(ctxt.get());
        if(!ic){
            xxBase::FatalTraceGlobal("ParFor called with incorrect context.");
        }
        ParFor(ic->range, ic->grainVolume, f);
    }

    /*
    Execute everything in the given range, but ensure each executing volume
    has an exclusion halo of one in every dimension.
    Simplest method is to split the range into blocks of 2x2x2, then run
    them as 8 groups, but with strange size volumes it is more complex.
    */
    virtual void ParForWithSafeHalo(
        const range_3d &range,
        unsigned grainVolume,
        std::function<void(const range_3d &)> f
    ) const override {
        f(range);
    }

    virtual void ParForWithSafeHalo(
        const Range3DToken &ctxt,
        std::function<void(const range_3d &)> f
    ) const override {
        auto ic=dynamic_cast<const Range3DContextImpl *>(ctxt.get());
        if(!ic){
            xxBase::FatalTraceGlobal("ParFor called with incorrect context.");
        }
        ParFor(ic->range, ic->grainVolume, f);
    }
};

#endif
