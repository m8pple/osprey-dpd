#ifndef ParallelServices_hpp
#define ParallelServices_hpp

#include <functional>
#include <cstdint>
#include <memory>

/*
These provide wrapper services around possible parallel providers, with
a base-case for sequential.
Because the back-ends are variable, they might not be that efficient, so
don't assume TBB levels of optimisation. It is recommended that tasks
be chunky enough to be efficient on almost any possible threading implementation.
Usually TBB wants at least 1us per task, but here we probably want to
aim for 10-100 us.
*/
class ParallelServices
{
public:
    struct range_1d
    {
        int64_t begin;
        int64_t end;

        template<class TFunc>
        void ForEach(TFunc &f)
        {
            for(int64_t x=begin; x<end; x++){
                f(x);
            }
        }
    };

    struct range_3d
    {
        int64_t begin[3];
        int64_t end[3];

        template<class TFunc>
        void ForEach(TFunc &f)
        {
            for(int64_t z=begin[2]; z<end[2]; z++){
                for(int64_t y=begin[1]; y<end[1]; y++){
                    for(int64_t x=begin[0]; x<end[0]; x++){
                        f(x,y,z);
                    }
                }   
            }
        }
    };

    virtual ~ParallelServices()
    {}

    static ParallelServices &GetServices();

    class RangeContext
    {
    public:
        virtual ~RangeContext()
        {}
    };

    class Range1DContext
        : public RangeContext
    {};

    class Range3DContext
        : public RangeContext
    {};

    using Range1DToken = std::unique_ptr<Range1DContext>;
    using Range3DToken = std::unique_ptr<Range3DContext>;

    virtual Range1DToken CreateParForToken(const range_1d& range)=0;
    virtual Range3DToken CreateParForToken(const range_3d& range)=0;
    virtual Range3DToken CreateParForWithSafeHaloToken(const range_3d& range)=0;

    virtual void Par(
        std::function<void()> f1,
        std::function<void()> f2
    ) =0;

    virtual void ParFor(
        const range_1d &range,
        std::function<void(const range_1d &)> f
    ) = 0;

    virtual void ParFor(
        const Range1DToken &ctxt,
        std::function<void(const range_1d &)> f
    ) =0;

    virtual void ParFor(
        const range_3d &range,
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
    them as 8 groups. With non-even size volumes we need to split the
    range into groups of 3 or 2 to preserve halo properties.
    */
    virtual void ParForWithSafeHalo(
        const range_3d &range,
        std::function<void(const range_3d &)> f
    ) = 0;

    virtual void ParForWithSafeHalo(
        const Range3DToken &range,
        std::function<void(const range_3d &)> f
    ) = 0;
};

#endif
