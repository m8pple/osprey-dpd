// This file is not intended to be platform portable. It only works under linux...
#include <unistd.h>
#include <fcntl.h>

// ... and it libtestu01 is installed
extern "C"{
    #include "testu01.h"
}

#include "BeadIdHashRNG.h"


#include <vector>
#include <set>
#include <utility>
#include <string>
#include <iostream>

uint32_t SplitMix32(uint32_t  s) {
	s += 2654435769UL;
	s ^= s >> 16;
	s *= 0x85ebca6b;
	s ^= s >> 13;
	s *= 0xc2b2ae35;
	s ^= s >> 16;
	return s;
}

struct bead_hash_t
{
    uint32_t hash;
};

class HashRNG
{
protected:
    HashRNG()
    {}
public:
    ~HashRNG()
    {}

    virtual void BeginRound(uint64_t seed, uint64_t time) =0;

    virtual bead_hash_t BeadHash(uint32_t x) =0;

    virtual uint32_t BeadHashPair(bead_hash_t hash1, bead_hash_t hash2) =0 ;
};




class BeadIdHashRNG final
    : public HashRNG
{
private:
    uint64_t round_seed;

public:
    BeadIdHashRNG()
    {}

    void BeginRound(uint64_t seed, uint64_t time)
    {
        round_seed = bead_id_hash_rng__round_hash(seed, time);
    }

    bead_hash_t BeadHash(uint32_t x)
    {
        return {x};
    }

    uint32_t BeadHashPair(bead_hash_t id1, bead_hash_t id2)
    {
        return bead_id_hash_rng__bead_hash(round_seed, id1.hash, id2.hash);
    }
};


class HashA final
    : public HashRNG
{
private:
    uint64_t round_seed;

public:
    HashA()
    {}

    void BeginRound(uint64_t seed, uint64_t time)
    {
        round_seed = bead_id_hash_rng__round_hash(seed, time);
    }

    bead_hash_t BeadHash(uint32_t x)
    {
        return {x};
    }

    uint32_t BeadHashPair(bead_hash_t id1, bead_hash_t id2)
    {
        uint32_t ma=std::min(id1.hash,id2.hash);
        uint32_t mb=std::max(id1.hash,id2.hash);

        ma += round_seed>>32;
        mb += round_seed&0xFFFFFFFFul;

        ma *= 0x85ebca6b;
        mb *= 0x85ebca6b;

        ma += (mb >> 13);
        mb += (ma >> 13);

        uint32_t m = ma ^ mb;

        m *= 0xc2b2ae35;

        return m;
    }
};


class HashB final
    : public HashRNG
{
private:
    uint64_t round_seed;

public:
    HashB()
    {}

    void BeginRound(uint64_t seed, uint64_t time)
    {
        round_seed = bead_id_hash_rng__round_hash(seed, time);
    }

    bead_hash_t BeadHash(uint32_t x)
    {
        x += (round_seed>>32);
        x *= 0x85ebca6b;
        x += (x>>16);
        x *= 0x85ebca6b;
        x += (x>>16);
        return {x};
    }

    uint32_t BeadHashPair(bead_hash_t id1, bead_hash_t id2)
    {
        uint32_t ma= std::min(id1.hash, id2.hash);
        uint32_t mb= std::max(id1.hash, id2.hash);

        ma *= 0x85ebca6b;
        mb *= 0x85ebca6b;

        ma += (mb >> 13) + (round_seed&0xFFFFFFFFul);
        mb += (ma >> 13) + (round_seed&0xFFFFFFFFul);

        ma *= 0xc2b2ae35;
        mb *= 0xc2b2ae35;

        uint32_t m = ma ^ mb;

        m *= 0xc2b2ae35;

        return m;
    }
};


class HashC final
    : public HashRNG
{
private:
    uint64_t round_seed;

public:
    HashC()
    {}

    void BeginRound(uint64_t seed, uint64_t time)
    {
        round_seed = bead_id_hash_rng__round_hash(seed, time);
    }

    bead_hash_t BeadHash(uint32_t x)
    {
        x += (round_seed>>32);
        x *= 0x85ebca6b;
        x += (x>>16);
        x *= 0x85ebca6b;
        x += (x>>16);
        return {x};
    }

    uint32_t BeadHashPair(bead_hash_t id1, bead_hash_t id2)
    {
        uint32_t ma= std::min(id1.hash, id2.hash);
        uint32_t mb= std::max(id1.hash, id2.hash);

        ma += (mb>>16) + (round_seed&0xFFFFFFFFul);
        mb += (ma>>16);
        mb *= 0x85ebca6b;

        uint32_t m = ma ^ mb;

        m *= 0x85ebca6b;
        return m;
    }
};


class HashD final
    : public HashRNG
{
private:
    uint64_t round_seed;

public:
    HashD()
    {}

    void BeginRound(uint64_t seed, uint64_t time)
    {
        round_seed = bead_id_hash_rng__round_hash(seed, time);
    }

    bead_hash_t BeadHash(uint32_t x)
    {
        x += (round_seed>>32);
        x *= 0x85ebca6b;
        x += (x>>16);
        x *= 0x85ebca6b;
        x += (x>>16);
        return {x};
    }

    uint32_t BeadHashPair(bead_hash_t id1, bead_hash_t id2)
    {
        uint32_t ma= std::min(id1.hash, id2.hash);
        uint32_t mb= std::max(id1.hash, id2.hash);

        const uint32_t C = 0b11111000011111000111011101;

        ma += (mb>>16) + (round_seed&0xFFFFFFFFul);
        ma *= C;

        mb += (ma>>16);
        mb *= C;
        
        uint32_t m = ma ^ mb;

        m *= C;
        return m;
    }
};


/*
class HashE final
    : public HashRNG
{
private:
    uint64_t round_seed;

public:
    HashE()
    {}

    void BeginRound(uint64_t seed, uint64_t time)
    {
        round_seed = bead_id_hash_rng__round_hash(seed, time);
    }

    bead_hash_t BeadHash(uint32_t x)
    {
        return { SplitMix32(x+(round_seed>>32)) };
    }

    uint32_t BeadHashPair(bead_hash_t id1, bead_hash_t id2)
    {
        uint32_t ma= id1.hash;
        uint32_t mb= id2.hash;

        //const uint32_t C = 0b11111000011111000111011101;
        const uint32_t C = 0x85ebca6b;

        uint32_t x=ma + (mb>>2);
        uint32_t y=mb + (ma>>2);

        ma=x + (y>>16);
        mb=y + (x>>16);

        ma *= C;
        mb *= C;

        uint32_t res=ma^mb;

        return res ^ (res<<16);
    }
};
*/

class HashF final
    : public HashRNG
{
private:
    uint64_t round_seed;

public:
    HashF()
    {}

    void BeginRound(uint64_t seed, uint64_t time)
    {
        round_seed = bead_id_hash_rng__round_hash(seed, time);
    }

    bead_hash_t BeadHash(uint32_t x)
    {
        return {SplitMix32(x+(round_seed>>32))};
    }

    uint32_t BeadHashPair(bead_hash_t id1, bead_hash_t id2)
    {
        uint32_t ma= std::min(id1.hash, id2.hash);
        uint32_t mb= std::max(id1.hash, id2.hash);

        uint32_t C=0x85ebca6b;

        /*
        Multiplier is A = (c<<32)+1
       
        Overall op is 

        (hi':lo') = A*(hi:lo)
        (hi':lo') ^= (hi':lo')>>48
        */

        uint32_t hi0=ma;
        uint32_t lo0=mb;

        uint32_t hi1=hi0 + lo0 * C;
        uint32_t lo1=lo0;

        lo1 ^= (hi1>>16);

        uint32_t hi2=hi1 + lo1 * C;
        uint32_t lo2=lo1;

        lo2 ^= (hi2>>16);

        uint32_t hi3=hi2 + lo2 * C;
        uint32_t lo3=lo2;
        
        return hi3;
    }
};


class ExhaustiveWalker
{
    HashRNG *rng;
    uint64_t global_seed;
    unsigned max_bead_id;
    unsigned round_id;
    uint32_t ida;
    uint32_t idb;
    bead_hash_t hash_b;


public:
    ExhaustiveWalker &operator=(const ExhaustiveWalker &) = default;
    ExhaustiveWalker(const ExhaustiveWalker &) = default;

    ExhaustiveWalker(HashRNG *_rng, uint64_t _global_seed, unsigned _max_bead_id)
        : rng(_rng)
        , global_seed(_global_seed)
        , max_bead_id(_max_bead_id)
        , round_id(0)
        , ida(0)
        , idb(1)
    {
        rng->BeginRound(global_seed, round_id);
        hash_b = rng->BeadHash(idb);
    }

    std::string Name()
    {
        return "Exhaustive[seed="+std::to_string(global_seed)+";max_beed="+std::to_string(max_bead_id)+"]";
    }

    uint32_t next()
    {
        uint32_t res=rng->BeadHashPair(rng->BeadHash(ida), hash_b);
        ++ida;
        if(ida == idb)
        {
            ida=0;
            ++idb;
            if(idb==max_bead_id)
            {
                idb=1;
                ++round_id;
                rng->BeginRound(global_seed, round_id);
            }
            hash_b = rng->BeadHash(idb);
        }
        return res;
    }
};


class SlidingWalker
{
    HashRNG *rng;
    uint64_t global_seed;
    unsigned max_bead_id;
    unsigned round_id;
    uint32_t ida;
    bead_hash_t hash_a;
public:
    SlidingWalker &operator=(const SlidingWalker &) = default;
    SlidingWalker(const SlidingWalker &) = default;

    SlidingWalker(HashRNG *_rng,uint64_t _global_seed, unsigned _max_bead_id)
        : rng(_rng)
        , global_seed(_global_seed)
        , max_bead_id(_max_bead_id)
        , round_id(0)
        , ida(0)
    {
        rng->BeginRound(global_seed, round_id);
        hash_a=rng->BeadHash(ida);
    }

    std::string Name()
    {
        return "Sliding[seed="+std::to_string(global_seed)+";max_beed="+std::to_string(max_bead_id)+"]";
    }

    uint32_t next()
    {
        bead_hash_t hash_b=rng->BeadHash(ida+1);
        uint32_t res=rng->BeadHashPair(hash_a, hash_b);
        ++ida;
        hash_a=hash_b;
        if(ida+1 == max_bead_id)
        {
            ida=0;
            ++round_id;
            rng->BeginRound(global_seed, round_id);  
            hash_a=rng->BeadHash(ida);
        }
        return res;
    }
};

/*
Points are positioned at integer points on a grid. Each bead's neighbourhood
is all beads within distance 2.
This results in a forwards neighbourhood of 16, and an implicit full neighbourhood of 32.
*/
class NeighbourWalker
{
    HashRNG *rng;
    uint64_t global_seed;
    unsigned width;
    unsigned height;
    unsigned depth;

    unsigned round_id;
    
    std::vector<std::pair<uint32_t,uint32_t> > pairs;
    unsigned pair_offset;

    std::vector<bead_hash_t> bead_hashes;

    std::vector<std::tuple<int,int,int>> make_forwards_nhood()
    {
        std::set<std::tuple<int,int,int>> nhood;

        for(int dx=-2; dx<=+2; dx++){
            for(int dy=-2; dy<=+2; dy++){
                for(int dz=-2; dz<=+2; dz++){
                    if(dx==0 && dy==0 && dz==0){
                        continue;
                    }
                    if( dx*dx+dy*dy+dz*dz > 2*2 ){
                        continue;
                    }
                    auto rev=std::tuple{-dx,-dy,-dz};
                    if(nhood.find(rev) == nhood.end()){
                        nhood.insert({dx,dy,dz});
                    }
                }
            }
        }
        return {nhood.begin(),nhood.end()};
    }

    uint32_t get_id(const std::tuple<int,int,int> &p)
    {
        return std::get<0>(p) + std::get<1>(p)*width + std::get<2>(p)*width*height;
    }

    void make_pairs()
    {
        auto nhood=make_forwards_nhood();
        //std::cerr<<"Nhood size="<<nhood.size()<<"\n";

        pairs.clear();
        for(unsigned x=0; x<width; x++){
            for(unsigned y=0; y<height; y++){
                for(unsigned z=0; z<height; z++){
                    uint32_t home=get_id({x,y,z});
                    for(auto [dx,dy,dz] : nhood){
                        uint32_t other = get_id({(x+width+dx)%width, (y+height+dy)%height, (z+depth+dz)%depth});
                        pairs.push_back({home,other});
                    }
                }
            }
        }
        //std::cerr<<"Pairs size = "<<pairs.size()<<"\n";
    }

    void BeginStep()
    {
        rng->BeginRound(global_seed, round_id);

        bead_hashes.resize( width*height*depth );
        for(unsigned i=0; i<bead_hashes.size(); i++){
            bead_hashes[i] = rng->BeadHash(i);
        }

        pair_offset=0;
        round_id += 1;
    }
public:
    NeighbourWalker &operator=(const NeighbourWalker &) = default;
    NeighbourWalker(const NeighbourWalker &) = default;

    NeighbourWalker(HashRNG *_rng, uint64_t _global_seed, unsigned _width, unsigned _height, unsigned _depth)
        : rng(_rng)
        , global_seed(_global_seed)
        , width(_width)
        , height(_height)
        , depth(_depth)
    {
        make_pairs();

        BeginStep();
    }

    std::string Name()
    {
        return "Neighbour[seed="+std::to_string(global_seed)+";w="+std::to_string(width)+";h="+std::to_string(height)+";d="+std::to_string(depth)+"]";
    }

    uint32_t next()
    {
        const auto &pair=pairs[pair_offset];
        uint32_t res = rng->BeadHashPair( bead_hashes[pair.first], bead_hashes[pair.second] );

        ++pair_offset;
        if(pair_offset==pairs.size()){
            BeginStep();
        }

        return res;
    }
};

template<class TSrc>
struct rng_src
    : unif01_Gen
{
    TSrc src;
    std::shared_ptr<std::string> pname;
    uint64_t lsb_rng = 12345662423242;

    static unsigned long GetBitsImpl(void *param, void *state)
    {
        auto me=(rng_src*)state;
        me->lsb_rng = me->lsb_rng * 6364136223846793005 + 1;
        return (me->src.next() & 0xFFFFFF00) | (me->lsb_rng>>56);
    }

    static double GetU01Impl(void *param, void *state)
    {
        return GetBitsImpl(param,state) * 0.00000000023283064365386962890625;
    }

    static void WriteImpl(void *state)
    {
        auto me=(rng_src*)state;
        fprintf(stdout, "%s", me->name);
    }

    rng_src(const TSrc &_src, const std::string &_name)
        : src(_src)
        , pname(std::make_shared<std::string>(_name))
    {
        name=(char*)pname->c_str();
    
        GetBits=GetBitsImpl;
        GetU01=GetU01Impl;
        state = this;
        param = 0;
        Write = WriteImpl;
    }

};

struct SupressStdout
{
    // https://stackoverflow.com/a/4832902

    int prev, curr;

    SupressStdout()
    {
        fflush(stdout);
        prev = dup(1);
        curr = open("/dev/null", O_WRONLY);
        dup2(curr, 1);
        close(curr);
    }

    ~SupressStdout()
    {
        fflush(stdout);
        dup2(prev, 1);
        close(prev);
    }

};

struct battery_res_t
{
    int nTest;
    int nWarn;
    int nFail;
};

battery_res_t RunSmallCrush(unif01_Gen *gen)
{
    {
        SupressStdout s;
        bbattery_SmallCrush(gen);
    }

    int fails=0, warns=0;
    for(int z=0; z<bbattery_NTests; z++){
        double p=bbattery_pVal[z];
        p=std::min(p, 1-p);
        if(p < 1e-8){
            fails += 1;
            warns += 1;
            fprintf(stdout, "  # %s : %.12f\n", bbattery_TestNames[z], p);
        }else if(p < 1e-5){
            warns += 1;
        }
    }

    return {bbattery_NTests, warns, fails};
}

void TestRNG(HashRNG *rng)
{
    std::vector<uint64_t> global_seeds({0,1,2,3, 1ull<<32, 1ull<<63});
    for(auto global_seed : global_seeds){
        for(int w=8; w<=32; w+=2){
            {
                NeighbourWalker walker(rng, global_seed, w,w,w);
                rng_src<NeighbourWalker> src( walker, "walker" );

                auto res=RunSmallCrush(&src);

                fprintf(stdout, " %s, %llu, %d, %d, %d\n", walker.Name().c_str(), (unsigned long long)global_seed, res.nTest, res.nFail, res.nWarn);
            }
        }

        for(int i=4; i<=(1<<28); i=(i*8)/5){
            {
                ExhaustiveWalker walker(rng, global_seed, i);
                rng_src<ExhaustiveWalker> src( walker, "walker" );

                auto res=RunSmallCrush(&src);

                fprintf(stdout, " %s, %llu, %d, %d, %d\n", walker.Name().c_str(), (unsigned long long)global_seed, res.nTest, res.nFail, res.nWarn);
            }

            {
                SlidingWalker walker(rng, global_seed, i);
                rng_src<SlidingWalker> src( walker, "walker" );

                auto res=RunSmallCrush(&src);

                fprintf(stdout, " %s, %llu, %d, %d, %d\n", walker.Name().c_str(), (unsigned long long)global_seed, res.nTest, res.nFail, res.nWarn);
            }
        }
    }
}

int main()
{
    swrite_Basic=0;

    HashF rng;
    TestRNG(&rng);
}
