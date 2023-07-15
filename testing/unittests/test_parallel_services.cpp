#include "ParallelContext.h"

#include <vector>
#include <algorithm>
#include <numeric>

#include <random>
#include <atomic>

unsigned test_num=0;

void test_parallel_context_par_tasks(ParallelContext &svc)
{
    std::atomic<int> count;
    count=0;

    std::function<void(unsigned)> rec = [&](unsigned n)
    {
        assert(n>0);
        if(n==1){
            count.fetch_add(1);
        }else{
            svc.ParTasks(
                [&](){
                    rec(n/2);
                },
                [&](){
                    rec(n-n/2);
                }
            );
        }
    };

    unsigned n=1000000;
    unsigned reps=10;
    for(int i=0; i<reps; i++){
        rec(n);
    }
    
    int sum=count.load();

    bool pass=(sum==n*reps);
    fprintf(stdout, "%s %d - test_parallel_context_par_tasks-%s\n", pass?"ok":"not ok", ++test_num, svc.Name());

    fflush(stdout);
}


void test_parallel_context_parfor_1d(ParallelContext &svc)
{
    std::atomic<int64_t> count;
    count=0;

    int64_t n=1000000;

    unsigned reps=10;
    int64_t want=0;
    for(int i=0; i<reps; i++){
        svc.ParFor({i,n}, [&](const range_1d &r){
            for(int64_t ii : r){
                count.fetch_add(ii);
            }
        });
        want += n*(n-1)/2  - i*(i-1)/2;
    }
    
    int64_t got=count.load();

    bool pass= (got==want);
    fprintf(stdout, "%s %d - test_parallel_context_parfor_1d-%s # got=%lld, want=%lld\n", pass?"ok":"not ok", ++test_num, svc.Name(), got, want);

    fflush(stdout);
}

void test_parallel_context_parfor_1d_token(ParallelContext &svc)
{
    std::atomic<int64_t> count;
    count=0;

    int64_t n=1000000;

    auto token = svc.CreateParForToken({10,n}, 3);

    unsigned reps=10;
    for(int i=0; i<reps; i++){
        svc.ParFor(token, [&](const range_1d &r){
            for(int64_t ii : r){
                count.fetch_add(ii);
            }
        });
    }
    
    int64_t got=count.load();
    int64_t want = ( n*(n-1)/2 - (10*(10-1)/2) ) * reps;

    bool pass= (got==want);
    fprintf(stdout, "%s %d - test_parallel_context_parfor_1d_token-%s # got=%lld, want=%lld\n", pass?"ok":"not ok", ++test_num, svc.Name(), got, want);

    fflush(stdout);
}

std::string test_parallel_context_parfor_3d(ParallelContext &svc, range_3d &range, unsigned grainSize)
{
    std::atomic<int64_t> count;
    count=0;

    std::atomic<bool> grain_failed;
    grain_failed=false;

    std::vector<unsigned> hits(range.volume(), 0);

    unsigned reps=10;
    int64_t want=0;
    for(int i=0; i<reps; i++){
        svc.ParFor(range, [&](const range_3d &r){
            for(int d=0; d<3; d++){
                if( r.axes[d].volume() < std::min<int64_t>( range.axes[d].volume(), grainSize ) ){
                    grain_failed=true;
                }
            }

            for(const auto p : r){
                count.fetch_add(1);
                hits[ range.to_linear_index(p) ] += 1;
            }
        });
        want += range.volume();
    }
    
    int64_t got=count.load();

    if(got!=want){
        return "count is wrong";
    }

    if(hits != std::vector<unsigned>(range.volume(), reps)){
        return "hits is wrong.";
    }

    return {};
}

void test_parallel_context_parfor_3d(ParallelContext &svc, std::mt19937_64 &rng)
{
    std::vector<int64_t> sizes{1,2,3,4,5,6,7,8,9,10};
    for(int i=0; i<10; i++){
        sizes.push_back( (rng() % 100)+10);
    }

    for(unsigned i=0; i<100; i++){
        range_3d r{{},{},{}};
        for(int d=0; d<3; d++){
            r.axes[d].lower = (rng()%100)-50;
            r.axes[d].upper = r.axes[d].lower + sizes[rng()%sizes.size()];
        }

        int grainSize=1+(rng()%4); // Allow grainSize of 1 to happen
        std::string res=test_parallel_context_parfor_3d(svc, r, grainSize);
        if(!res.empty()){
            fprintf(stdout, "not ok %d - test_parallel_context_parfor_3d-%s # %s\n", ++test_num,  svc.Name(), res.c_str());    
            return;
        }
    }

    fprintf(stdout, "ok %d - test_parallel_context_parfor_3d-%s\n", ++test_num, svc.Name());    
    fflush(stdout);
}

std::string test_parallel_context_parfor_3d_halo(ParallelContext &svc, const range_3d &range, unsigned grainSize)
{
    std::atomic<int64_t> count, count2;
    count=0;
    count2=0;

    std::vector<unsigned> hits(range.volume(), 0);

    std::atomic<bool> grainFailed;
    grainFailed=false;

    int64_t G=std::max<int64_t>(2, grainSize);

    unsigned reps=3;
    int64_t want=0;
    for(int i=0; i<reps; i++){
        std::function<void(const range_3d &range)> ff = [&](const range_3d &r){
            for(int d=0; d<3; d++){
                if(r.axes[d].volume() < std::min<int64_t>(range.axes[d].volume(), G)){
                    grainFailed.store(true);
                }
            }
            
            for(const auto p : r){
                count.fetch_add(1);
                range_3d halo{ {p.x-1, p.x+2}, {p.y-1, p.y+2}, {p.z-1,p.z+2} };
                for(const auto h : halo){
                    auto pr=range.wrap( h );
                    hits[ range.to_linear_index(pr) ] += 1;
                    count2.fetch_add(1);
                }
            }
        };

        svc.ParForWithSafeHalo(range, grainSize, ff);
        want += range.volume();
    }
    
    int64_t got=count.load();
    if(got!=want){
        return "count, got="+std::to_string(got)+", want="+std::to_string(want);
    }

    got=count2.load();
    want *= 27;
    if(got!=want){
        return "count2, got="+std::to_string(got)+", want="+std::to_string(want);
    }

    auto sum_hits=std::accumulate( hits.begin(), hits.end(), 0ll );
    if( hits != std::vector<unsigned>(range.volume(), reps*27) ){
        return "hits wrong";
    }

    if(grainFailed.load()){
        return "grain size failed.";
    }

    return {};
}


void test_parallel_context_parfor_3d_halo(ParallelContext &svc, std::mt19937_64 &rng)
{
    std::vector<int64_t> sizes{1,2,3,4,5,6,7,8,9,10};
    for(int i=0; i<10; i++){
        sizes.push_back( (rng() % 100)+10);
    }

    for(unsigned i=0; i<100; i++){
        range_3d r{{},{},{}};
        for(int d=0; d<3; d++){
            r.axes[d].lower = (rng()%100)-50;
            r.axes[d].upper = r.axes[d].lower + sizes[rng()%sizes.size()];
        }

        int grainSize=1+(rng()%4); // Allow grainSize of 1 to happen
        std::string res=test_parallel_context_parfor_3d_halo(svc, r, grainSize);
        if(!res.empty()){
            fprintf(stdout, "not ok %d - test_parallel_context_parfor_3d_halo-%s # %s\n", ++test_num,  svc.Name(), res.c_str());    
            return;
        }
    }

    fprintf(stdout, "ok %d - test_parallel_context_parfor_3d_halo-%s\n", ++test_num, svc.Name());    
    fflush(stdout);
    
}

void test_parallel_context(ParallelContext &svc)
{
    std::mt19937_64 rng(1);

    test_parallel_context_par_tasks(svc);
    test_parallel_context_parfor_1d(svc);
    test_parallel_context_parfor_1d_token(svc);
    test_parallel_context_parfor_3d(svc, rng);
    test_parallel_context_parfor_3d_halo(svc, rng);
}

int main()
{
    for(const auto &dut : ParallelContext::GetContexts()){
        test_parallel_context(*dut);
    }
}
