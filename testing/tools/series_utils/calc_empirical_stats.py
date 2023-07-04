import numpy as np
import time
import sys

def calc_uniform_ad_stat(paths:np.ndarray) -> np.ndarray:
    """
    Calculates anderson darling statistics
    - paths.shape[0]: number of points in paths
    - paths.shape[1]: number of paths
    """
    n=paths.shape[0]
    reps=paths.shape[1]
    inc=np.array( range(0,n), dtype=np.float32 ).reshape((n,1))
    inc=(2*inc+1) / n

    paths=np.sort(paths, axis=0)
    
    parts=inc * np.log(paths) * np.log(1-np.flip(paths,axis=0))
    stats=-n-np.sum(parts, axis=0)
    return stats

_ad_distribution_cache={}

def calc_uniform_ad_pvalue(x:np.ndarray) -> np.ndarray:
    reps=1000000

    assert len(x.shape)==1
    n=x.shape[0]
    points = _ad_distribution_cache.get(n, None)
    if points is None:
        points=calc_uniform_ad_distribution(n, reps)
        points.sort()
        _ad_distribution_cache[n]=points

    return np.digitize(x, points, right=True) / reps


rng=np.random.default_rng()

def calc_uniform_ad_distribution(n:int, m:int) -> np.ndarray:
    start=time.time_ns()
    u=rng.random(size=(n,m))
    s=calc_uniform_ad_stat(u)
    s.sort()
    finish=time.time_ns()
    sys.stderr(f"Took {((finish-start)*1e-9)} secs for {n} * {m}\n")

