
#include <cstdint>

/*
    The world is split into cells.
    Each cell has:
    - A unique linear id
        - These are contiguous, regardless of world dimensions
    - A unique 3D integer location

    Cells are grouped into SequentialClusters
    All cells within an SequentialCluster are processed by the same thread.

    SequentialClusters are grouped into NearClusters
    All cells within a NearCluster are assumed to be "local" in some sense, so that:
    - Atomic memory accesses are "close"
    - Cached memory accesses are shared, e.g. at L3

    within SequentialCluster -> cells running within one TBB task or a thread, no locks needed
    within NearCluster -> clusters running in one CPU group, cheap locks (e.g. atomics) appropriate
    

*/
class ClusterPolicy
{
public:
    struct CellInfo
    {
        uint32_t cell_x : 10;
        uint32_t cell_y : 10;
        uint32_t cell_z : 10;
        uint32_t cell_id;
        uint32_t seq_cluster_id;
        uint32_t local_cluster_id;
    };

    virtual const CellInfo &GetCellInfo(uint32_t cell_id) const =0;
    virtual const CellInfo &GetCellInfo(uint32_t pos[3]) const =0;


};