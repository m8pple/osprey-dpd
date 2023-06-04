
#include "ParallelContext.h"
#include "ParallelContext_Sequential.h"
#include "ParallelContext_NaiveThreads.h"

#ifdef OSPREY_DPD_ENABLE_TBB
#include "ParallelContext_TBB.h"
#endif

static const ParallelContext *sm_default = 0;

static std::vector<std::shared_ptr<ParallelContext> > sm_services;

const ParallelContext &ParallelContext::GetDefaultContext()
{
    assert(sm_default); // We should always have the sequential one
    return *sm_default;
}

const std::vector<std::shared_ptr<ParallelContext> > &ParallelContext::GetContexts()
{
    return sm_services;
}


std::shared_ptr<ParallelContext> ParallelContext::GetContext(const std::string &name)
{
    for(const auto &x : sm_services){
        if(name == x->Name()){
            return x;
        }
    }
    return {};
}

bool ParallelContext::RegisterContext(std::shared_ptr<ParallelContext> ctxt)
{
    for(auto &x : sm_services){
        if( !strcmp(x->Name(), ctxt->Name())){
            fprintf(stderr, "Attempt to register services twice for '%s'\n", ctxt->Name());
            return false;
        }
    }

    if(!ctxt->OnRegister()){
        return false;
    }

    sm_services.push_back(ctxt);
    if(sm_default==0 || ctxt->Merit() > sm_default->Merit()){
        sm_default=ctxt.get();
    }

    return true;
}

/* Split the range into an even number of blocks such that:
    - The number of blocks is even, or there is just one block
    - Each block has size at least two, or there is just one block
    - Ideally: there are at least four blocks, to expose parallelism 
    - Ideally: all blocks are within size 1 of each other, for load balancing
    - Ideally: Larger blocks outnumber smaller blocks, for load balancing

    Assuming grainSize=2:
    n=1, A
    n=2, AA
    n=3, AAA
    n=4, AABB  <- First true split, though it is pointless for parallel implementations
    n=5, AAABB
    n=6, AAABBB
    n=7, AAAABBB
    n=8, AABBAABB <- Next increase in block bount, at this point we have parallelism of two
    n=9, AAABBAABB
    n=10, AAABBBAABB
    n=11, AAABBBAAABB
    n=12, AABBAABBAABB <- Next increase in block count, parallelism of three

    So:
    - If n < 4, there is just one block
    - If n is divisible by 4, we just split into blocks of size 2
    - Otherwise, we get as many blocks of size 3 as possible, then the remainer of size 2

Let nThree = number of block size 3, nTwo = number of blocks size 2

n = nThree * 3 + nTwo * 2,  where (nThree+nTwo) is even


Warning: this has got way over-complicated, and hasn't been properly performance tested...

*/
std::vector<range_1d> ParallelContext::SplitForSafeHalo(const range_1d &r, unsigned grainSize)
{
    std::vector<range_1d> res;

    int64_t G=std::max<int64_t>(2, grainSize);

    auto n=r.volume();
    if(n<4*G){
        // Degenerate or no parallelism - just return one block
        res.push_back({r.lower, r.upper});
    }else if( (n % (2*G)) == 0 ) {
        // Simple but fairly common case, e.g. G=2 and ( n is a binary power or divisible by 4
        int64_t base=r.lower;
        while(base < r.upper){
            res.push_back({base, base+G});
            base += G;
        }
    }else{
        int nSmall=0; // Number of size G
        int nBig=0; // Number of size G+1

        // TODO : If G is large this could get slow?
        // But then if G is large then the things getting scheduled are very
        // big...
        // There is probably a way of jumping straight to the best point, but
        // it would require some maths. All we are looking for is:
        //   min nSmall  such that  G*nSmall + (G+1)*nBig ==n  and  (nSmall+nBig) is even
        //          where nBig = (n - G*nSmall)/(G+1)
        //  
        //  Subst  nBig = (n - G*nSmall)/(G+1)  into  G*nSmall + (G+1)*nBig ==n
        //        G*nSmall + (G+1)*(n - G*nSmall)/(G+1) == n
        //        G*nSmall + (n - G*nSmall) == n
        //        G*nSmall + n - G*nSmall == n
        //        n == n
        //        Huzzah!
        //
        // Condition for this is that nBig is an integer, so a constraint
        // is that (G+1) divides (n - G*nSmall)
        // 
        // Constraint that (nSmall+nBig) is even means that
        //                 (nSmall + (n - G*nSmall)/(G+1)) is even
        //
        // So there is a way to do this more intelligently, I'm just too lazy right now

        bool fail=false;
        while(1){
            nBig = (n - G*nSmall)/(G+1);
            DEBUG_ASSERT(nBig>0);
            if(nBig <= 0){
                fail=true; // I don't _think_ this can happen, but quit gracefully just in case
                break;
            }

            if( G*nSmall + (G+1)*nBig == n && (((nSmall+nBig)&1)==0) ){
                break; // good exit
            }
            ++nSmall;

            // It is not expected to take this long
            if(nSmall > G*G || nSmall > 100){
                DEBUG_ASSERT(false); // This is not an expected code path
                fail=true;
                break;
            }
        }

        if(!fail){
            int64_t base=r.lower;
            for(int i=0; i<nBig; i++){
                res.push_back( { base, base+G+1 } );
                base += G+1;
            }
            for(int i=0; i<nSmall; i++){
                res.push_back( { base, base+G } );
                base += G;
            }
            DEBUG_ASSERT(base==r.upper);
        }else{
            // Do something sensible in case of failure

            std::function<void(const range_1d &)> recurse=[&](const range_1d &rr){
                DEBUG_ASSERT(rr.volume() >= 2*G);
                auto rrr=rr.split();
                if( rr.volume() < 4 * G ){
                    res.push_back(rrr.first);
                    res.push_back(rrr.second);
                }else{
                    recurse(rrr.first);
                    recurse(rrr.second);
                }
            };

            DEBUG_ASSERT(r.volume() >= 2*G);
            recurse(r);
        }
    }

    DEBUG_ASSERT( (res.size()%2)==0 );
    DEBUG_ASSERT( res[0].lower==r.lower );
    for(int i=0; i<res.size(); i++){
        DEBUG_ASSERT(res[i].volume() >= 2);
        if(i>0){
            DEBUG_ASSERT(res[i-1].upper == res[i].lower); 
        }
    }
    DEBUG_ASSERT( res.back().upper == r.upper);
    return res;
}


static bool regSequential=ParallelContext::RegisterContext(std::make_shared<ParallelContext_Sequential>());
static bool regNaiveThreads=ParallelContext::RegisterContext(std::make_shared<ParallelContext_NaiveThreads>());
#ifdef OSPREY_DPD_ENABLE_TBB
static bool regTBB=ParallelContext::RegisterContext(std::make_shared<ParallelContext_TBB>());
#endif