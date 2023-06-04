#ifndef ParallelContext_hpp
#define ParallelContext_hpp

#include <functional>
#include <cstdint>
#include <memory>
#include <array>
#include <string>

#include "DebugAssert.h"

struct range_1d
{
    int64_t lower;
    int64_t upper;

    struct const_iterator
    {
        int64_t x;

        bool operator!=(const const_iterator &o) const
        { return x!=o.x; }

        const_iterator &operator++()
        {
            ++x;
            return *this;
        }

        int64_t operator*() const
        { return x; }
    };

    const_iterator begin() const
    { return {lower}; }

    const_iterator end() const
    { return {upper}; }

    template<class TFunc>
    void ForEach(TFunc &f)
    {
        for(int64_t x : *this){
            f(x);
        }
    }

    int64_t volume() const
    {
        assert(lower <= upper);
        return upper-lower;
    }

    std::pair<range_1d,range_1d> split() const
    {
        return { {lower, lower+(upper-lower)/2}, {lower+(upper-lower)/2, upper}};
    }
};

struct point_3d
{
    int64_t x, y, z;

    bool operator==(const point_3d &o)  const
    { return !memcmp(this, &o, sizeof(*this)); }

    bool operator!=(const point_3d &o)  const
    { return memcmp(this, &o, sizeof(*this)); }
        
    int64_t operator[](unsigned d) const
    { return d==0 ? x : d==1 ? y : z; }

    int64_t &operator[](unsigned d)
    { return d==0 ? x : d==1 ? y : z; }
};
static_assert(sizeof(point_3d) == sizeof(int64_t)*3, "Assuming that array in union will alias over points...");

struct range_3d
{
    range_1d axes[3];

    range_3d(const range_1d &x, const range_1d &y, const range_1d &z)
        : axes{x,y,z}
    {}

    struct const_iterator
    {
        const range_3d &parent;
        point_3d point;

        void assert_valid() const
        {
            DEBUG_ASSERT(parent.contains(point) || point == parent.end().point );
        }

        bool operator!=(const_iterator &o) const {
            assert_valid();
            return point != o.point;
        }

        const_iterator &operator++() {
            assert_valid();
            if(++point[0] == parent.axes[0].upper){
                point[0]=parent.axes[0].lower;
                if(++point[1] == parent.axes[1].upper){
                    point[1]=parent.axes[1].lower;
                    ++point[2]; 
                }
            }
            assert_valid();
            return *this;
        }

        const point_3d &operator*() const
        {
            assert_valid();
            return point;
        }
    };

    bool contains(const point_3d &p) const
    {
        for(int d=0; d<3; d++){
            if(p[d] < axes[d].lower || axes[d].upper <= p[d]) return false;
        }
        return true;
    }

    uint64_t to_linear_index(const point_3d &p) const
    {
        assert(contains(p));
        //uint64_t res= x.x + x.y * w + x.z * w * h;
        int64_t x=p.x-axes[0].lower, y=p.y-axes[1].lower, z=p.z-axes[2].lower;
        DEBUG_ASSERT( x < axes[0].volume());
        DEBUG_ASSERT( y < axes[1].volume());
        DEBUG_ASSERT( z < axes[2].volume());
        auto res = x + axes[0].volume() * (y + axes[1].volume() * (z) );       
        DEBUG_ASSERT( res < volume() );
        return res; 
    }

    point_3d lower_corner() const
    { return {axes[0].lower, axes[1].lower, axes[2].lower}; }

    const_iterator begin() const
    { return {*this, lower_corner()}; }

    point_3d upper_corner() const
    { return {axes[0].lower, axes[1].lower, axes[2].upper}; }

    const_iterator end() const
    { return {*this, upper_corner()}; }

    template<class TFunc>
    void ForEach(TFunc &&f) const
    {
        for(int64_t z : axes[2]){
            for(int64_t y : axes[1]){
                for(int64_t x : axes[0]){
                    f(x,y,z);
                }
            }   
        }
    }

    int64_t volume() const
    { return axes[0].volume() * axes[1].volume() * axes[2].volume(); }

    std::pair<int,int64_t> longest_dimension_and_length() const
    {
        int64_t len=-1;
        int dim=-1;
        for(int d=0; d<3; d++){
            auto curr=axes[d].volume();
            if(curr > len){
                dim=d;
                len=curr;
            }
        }
        return {dim,len};
    }

    int longest_dimension() const
    { return longest_dimension_and_length().first; }

    std::pair<range_3d,range_3d> split(int dimension) const
    {
        range_3d a{*this}, b{*this};
        a.axes[dimension].upper=axes[dimension].lower + (axes[dimension].upper-axes[dimension].lower)/2;
        b.axes[dimension].lower=a.axes[dimension].upper;
        return {a,b};
    }

    point_3d wrap(const point_3d &p) const
    {
        point_3d res(p);
        for(int d=0; d<3; d++){
            auto &v=res[d];
            if(v < axes[d].lower) v += (axes[d].upper - axes[d].lower);
            else if(v >= axes[d].upper) v -= (axes[d].upper - axes[d].lower);
        }        
        DEBUG_ASSERT(contains(res));
        return res;
    }
};

/*
These provide wrapper services around possible parallel providers, with
a base-case for sequential.
Because the back-ends are variable, they might not be that efficient, so
don't assume TBB levels of optimisation. It is recommended that tasks
be chunky enough to be efficient on almost any possible threading implementation.
Usually TBB wants at least 1us per task, but here we probably want to
aim for 10-100 us.

Parallel contexts are always expected to be singletons, so only one
instance of each context ever exists.

Parallel contexts are registered before main(), at which point some
contexts may or may not get loaded (depending on available features)
*/
class ParallelContext
{
public:

    virtual ~ParallelContext()
    {}

    // Highest numerical merit is the preferred services
    static const ParallelContext &GetDefaultContext();

    static std::shared_ptr<ParallelContext> GetContext(const std::string &name);
    static const std::vector<std::shared_ptr<ParallelContext> > &GetContexts();
    static bool RegisterContext(std::shared_ptr<ParallelContext> services);

    // Gives the parallel context a chance to do any one-off registration when
    // it is registered.
    // This avoids needing to mess around with too many static initialisation tricks.
    // Note that parallel contexts are supposed to be singletons, so a given
    // concrete class only gets registered once.
    // Return false to indicate the class failed to initialise
    virtual bool OnRegister()
    {
        return true;
    }

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


    virtual const char *Name() const=0;

    /* Merit is how efficient this method is. Higher is better.
    Suggested ranges are:
    - -1 : This is some sort of debug engine or trial implementation 
    - 0 : This is not parallel, but is ok sequential
    - 10 : Parallel, but not that efficient at it
    - 100 : Parallel and pretty efficient
    - 1000 : Parallel, efficient, and with extra features (e.g. optimisation or locality)
    */
    virtual int Merit() const=0;

    using Range1DToken = std::unique_ptr<Range1DContext>;
    using Range3DToken = std::unique_ptr<Range3DContext>;

    virtual Range1DToken CreateParForToken(const range_1d& range, unsigned grainSize=1) const=0;
    virtual Range3DToken CreateParForToken(const range_3d& range, unsigned grainVolume=1 ) const=0;
    virtual Range3DToken CreateParForWithSafeHaloToken(const range_3d& range, unsigned grainVolume) const=0;

    virtual void ParTasks(
        std::function<void()> f1,
        std::function<void()> f2
    ) const =0;

    virtual void ParFor(
        const range_1d &range,
        unsigned grainSize,
        std::function<void(const range_1d &)> f
    ) const = 0;

    void ParFor(
        const range_1d &range,
        std::function<void(const range_1d &)> f
    ) const {
        ParFor(range, 1, f);
    }

    virtual void ParFor(
        const Range1DToken &ctxt,
        std::function<void(const range_1d &)> f
    ) const =0;

    /*
    For simplicit we define just one grain size which applies to all dimensions.
    So for example if we say grainDim=2, we are saying don't split below a volume of 8.
    */
    virtual void ParFor(
        const range_3d &range,
        unsigned grainDim,
        std::function<void(const range_3d &)> f
    ) const = 0;

    void ParFor(
        const range_3d &range,
        std::function<void(const range_3d &)> f
    ) const {
        ParFor(range, 1, f);
    }

    virtual void ParFor(
        const Range3DToken &range,
        std::function<void(const range_3d &)> f
    ) const = 0;

    /*
    Execute everything in the given range, but ensure each executing volume
    has an exclusion halo of one in every dimension.
    grainDim defines the smallest size in each dimension for a parallel group,
    and for safe halos this must be at least 2. So if grainDim=1 is passed,
    then it will either be bumped up to 2, or there will be no parallelism
    (depending on the implementation). Let G = max(2, grainDim)

    For each dimension we can split it into an even number of blocks of G+1 or G.
    We then execute 8 rounds of parallel 3d execution:
       0 : odd x, odd y, odd z
       1 : odd x, odd y, even z
       2 : odd x, evey y, odd z
       ...
       7 : even x, even y, even z

    A dimension is degenerate (one block), if length() < 2*G
    To provide two blocks in each dimension, we need length() >= 2*G
    To provide parallel blocks in each dimension we need length() >= 4*G

    If one dimension is degenerate (only one block) or non-parallel (two blocks), an implementation
    might give up and run things sequentially. This might happen with extremely long thin
    tubes, or very flat grids. If it becomes an issue then implementations could be
    upgraded to deal with it.
    */
    virtual void ParForWithSafeHalo(
        const range_3d &range,
        unsigned grainDim,
        std::function<void(const range_3d &)> f
    ) const = 0;

    void ParForWithSafeHalo(
        const range_3d &range,
        std::function<void(const range_3d &)> f
    ) const {
        ParForWithSafeHalo(range, 2, f);
    }

    virtual void ParForWithSafeHalo(
        const Range3DToken &range,
        std::function<void(const range_3d &)> f
    ) const = 0;

    /* Split the range into contiguous blocks such that:
    - The number of blocks is even, or there is just one block
    - Each block has size at least grainSize, or there is just one block
    - Ideally: there are at least four blocks, to expose parallelism 
    - Ideally: all blocks are within size 1 of each other, for load balancing
    - Ideally: Larger blocks outnumber smaller blocks, for load balancing

    Larger blocks outnumbering smaller blocks is to deal with the case where
    we split a cube using a grainSize of 2, and we end up with
    each dimension having a single block of size 3. That means there will
    be one task of volume 27, a small number with volume 18, many with volume 10, and most with size 8
    If we are unlucky, then the big tasks go last, and hold up the group.
    If we do mostly size 3, and one of size 2, we end up with
    one task of 8, a few with volume 10, many with volume 18, and most with size 27.

    This needs to be balanced against creating the number of parallel tasks though...
    Might be over-thinking it, but this is a case where dynamic splitting from
    TBB doesn't help as we have already manually split it.
    */
    static std::vector<range_1d> SplitForSafeHalo(const range_1d &r, unsigned grainSize);
};

#endif
