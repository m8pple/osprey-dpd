#include "ParallelServices.hpp"

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
        range_1d range;
    };

    struct Range3DContextImpl
        : public Range3DContext
    {
        range_3d range;
    };

    virtual Range1DToken CreateParForToken(const range_1d& range) override
    {
        return std::make_unique<Range1DContextImpl>(range);
    }
    
    virtual Range3DToken CreateParForToken(const range_3d& range) override
    {
        return std::make_unique<Range3DContextImpl>(range);
    }

    virtual Range3DToken CreateParForWithSafeHaloToken(const range_3d& range) override
    {
        return std::make_unique<Range3DContextImpl>(range);
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
    ) {
        auto ic=dynamic_cast<const Range1DContextImpl *>(&ctxt);
        if(!ic){

        }
    }

    virtual void ParFor(
        const range_1d &range,
        std::function<void(const range_3d &)> f
    ) = 0;

    virtual void ParFor(
        const Range3DToken &range,
        std::function<void(const range_3d &)> f
    ) = 0;

    /*
    Execute everything in the given range, but ensure each executing volume
    has an exclusion halo of one in every dimension.
    Simplest method is to split the range into blocks of 2x2x2, then run
    them as 8 groups, but with strange size volumes it is more complex.
    */
    virtual void ParForWithSafeHalo(
        const range_3d &range,
        std::function<void(const range_3d &)> f
    ) = 0;

    virtual void ParForWithSafeHalo(
        const Range3DContext &range,
        std::function<void(const range_3d &)> f
    ) = 0;
};