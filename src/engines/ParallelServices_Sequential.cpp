#include "ParallelServices.hpp"

#include "xxBase.h"

class ParallelServices_Sequential
    : public ParallelServices
{
public:
    class RangeContext
    {
    public:
        virtual ~RangeContext()
        {}
    };

    struct Range1DContextImpl
        : public Range1DContext
    {
        Range1DContextImpl(const range_1d &_r)
            : range(_r)
        {}

        range_1d range;
    };

    struct Range3DContextImpl
        : public Range3DContext
    {
        Range3DContextImpl(const range_3d &_r)
            : range(_r)
        {}

        range_3d range;
    };

    virtual Range1DToken CreateParForToken(const range_1d& range) override
    {
        return Range1DToken( new Range1DContextImpl(range) );
    }
    
    virtual Range3DToken CreateParForToken(const range_3d& range) override
    {
        return Range3DToken( new Range3DContextImpl(range) );
    }

    virtual Range3DToken CreateParForWithSafeHaloToken(const range_3d& range) override
    {
        return Range3DToken( new Range3DContextImpl(range) );
    }

    virtual void Par(
        std::function<void()> f1,
        std::function<void()> f2
    ) override
    {
        f1();
        f2();
    }

    virtual void ParFor(
        const range_1d &range,
        std::function<void(const range_1d &)> f
    ) override
    {
        f(range);
    }

    virtual void ParFor(
        const Range1DToken &ctxt,
        std::function<void(const range_1d &)> f
    ) override {
        auto ic=dynamic_cast<const Range1DContextImpl *>(ctxt.get());
        if(!ic){
            xxBase::FatalTraceGlobal("ParFor called with incorrect context.");
        }
        ParFor(ic->range, f);
    }

    virtual void ParFor(
        const range_3d &range,
        std::function<void(const range_3d &)> f
    ) override {
        f(range);
    }

    virtual void ParFor(
        const Range3DToken &ctxt,
        std::function<void(const range_3d &)> f
    ) override {
        auto ic=dynamic_cast<const Range3DContextImpl *>(ctxt.get());
        if(!ic){
            xxBase::FatalTraceGlobal("ParFor called with incorrect context.");
        }
        ParFor(ic->range, f);
    }

    /*
    Execute everything in the given range, but ensure each executing volume
    has an exclusion halo of one in every dimension.
    Simplest method is to split the range into blocks of 2x2x2, then run
    them as 8 groups, but with strange size volumes it is more complex.
    */
    virtual void ParForWithSafeHalo(
        const range_3d &range,
        std::function<void(const range_3d &)> f
    ) override {
        f(range);
    }

    virtual void ParForWithSafeHalo(
        const Range3DToken &ctxt,
        std::function<void(const range_3d &)> f
    ) override {
        auto ic=dynamic_cast<const Range3DContextImpl *>(ctxt.get());
        if(!ic){
            xxBase::FatalTraceGlobal("ParFor called with incorrect context.");
        }
        ParFor(ic->range, f);
    }
};