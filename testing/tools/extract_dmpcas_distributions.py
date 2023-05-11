#!/usr/bin/env python3 

import sys
import os
import shutil
import subprocess
import re
import glob
from typing import *
import numpy as np
import json
import math

from series_utils.time_series_bundle import TimeSeriesBundle

from series_utils.parse_dmpcas import parse_dmpcas
from run_replicates import make_dir_name, hash_text

def calculate_quantised_walk_statistic(nbuckets:int, ntimes:int, nreps:int, boundaries : np.ndarray, transition_probs : np.ndarray, paths : np.ndarray):
    assert transition_probs.shape==(ntimes,nbuckets,nbuckets)    
    reshaped_probs=np.zeros(shape=(ntimes,nbuckets*nbuckets))
    for b in range(nbuckets):
        reshaped_probs[:,b*nbuckets:(b+1)*nbuckets]=nbuckets / transition_probs[:,b,:]
    
    assert paths.shape==(ntimes,nreps) or paths.shape==(ntimes,)
    indices=np.zeros(shape=paths.shape,dtype=int)
    for t in range(ntimes):
        indices[t,:] = np.digitize(paths[t,:], boundaries[t,:], right=True)
    
    hop_indices = indices[:-1,:] * nbuckets + indices[1:,:] #type: np.ndarray
    assert hop_indices.shape == (ntimes-1,nreps)

    walk_probs = np.zeros(shape=(ntimes-1,nreps))
    for t in range(ntimes-1):
        pp=reshaped_probs[t+1][hop_indices[t]]
        assert not np.any(pp==0), f"t={t}"
        walk_probs[t,:] = pp

    assert not np.any(walk_probs==0)

    log_probs = np.log(walk_probs)
    res=np.sum(log_probs, axis=0)
    return res

def calculate_extreme_value_statistic(ntimes:int, nreps:int, ul_boundaries : np.ndarray, paths : np.ndarray):
    assert ul_boundaries.shape == (ntimes,2)
    assert paths.shape==(ntimes,nreps)
    assert np.sum(ul_boundaries[:,0] < ul_boundaries[:,1]) > ntimes/2, str(ul_boundaries)+"\n"+str(paths)
    extrema = (paths < ul_boundaries[:,0].reshape((ntimes,1))) | (ul_boundaries[:,1].reshape((ntimes,1)) < paths)
    statistics=np.sum(extrema, axis=0) / ntimes
    assert statistics.shape==(nreps,)
    return statistics

if __name__ == "__main__":

    name_to_slices={} #type: Mapping[int,Dict[str,List[float]]]
    rep_to_series={} #type: Mapping[int,Dict[str,List[float]]]

    dmpci_file=sys.argv[1]
    with open(dmpci_file,"r") as src:
        dmpci_contents=src.read()
    dmpci_hash=hash_text(dmpci_contents)
    m=re.match("dmpci[.]([-_a-zA-Z0-9]+)", dmpci_file)
    assert m, "Couldn't parse '{dmpci_file}' into 'dmpci.NAME' "
    dmpci_name=m.group(1)
    sys.stderr.write(f"dmpci {dmpci_name}, hash {dmpci_hash}\n")


    if len(sys.argv) < 3:
        working_dir=f"./working"
        os.makedirs(working_dir,exist_ok=True)
    else:
        working_dir=sys.argv[2]

    if len(sys.argv) < 4:
        output_dir=f"./data"
        os.makedirs(working_dir,exist_ok=True)
    else:
        output_dir=sys.argv[3]

    if len(sys.argv) < 5:
        start=0
    else:
        start=int(sys.argv[4])

    sys.stderr.write(f"Looking for results for {dmpci_name}-{dmpci_hash} in directory {working_dir}, starting at index {start}\n")
    sys.stderr.write(f"Writing to output dir {output_dir}")

    all_slices=[]

    r=start
    while True:
        replicate_dir=f"{working_dir}/{make_dir_name(r)}"
        if not os.path.isdir(replicate_dir):
            sys.stderr.write(f"Missing directory {replicate_dir} finished search\n")
            break
        
        if not os.path.exists(f"{replicate_dir}/finished.{dmpci_name}-{dmpci_hash}"):
            sys.stderr.write(f"Skipping {replicate_dir}\n")
        else:
            slice=parse_dmpcas(replicate_dir, f"{dmpci_name}-{dmpci_hash}")
            all_slices.append(slice)
        r=r+1

    assert(len(all_slices) > 0)

    times=[t for t in all_slices[0].keys()]
    names=list(all_slices[0][times[0]].keys())
    reps=range(0,len(all_slices))


    nnames=len(names)
    nreps=len(reps)
    ntimes=len(times)

    res={}

    res["dmpci_name"]=dmpci_name
    res["dmpci_hash"]=dmpci_hash
    res["names"]=names
    res["nreps"]=nreps
    res["times"]=times

    res["slices"]={ name:{ "times":[] } for name in names }

    data_cube=np.zeros(shape=[nnames,ntimes,nreps])
    for (rindex,slice) in enumerate(all_slices):
        for (tindex,(t,vals)) in enumerate(slice.items()):
            for (n,v) in vals.items():
                nindex=names.index(n)
                data_cube[nindex,tindex,rindex]=v

    stats_cube=np.zeros(shape=[nnames,ntimes,2])

    means = np.zeros([nnames,ntimes])
    stddevs = np.zeros([nnames,ntimes])

    nbuckets=5

    transition_probs=np.zeros(shape=(nnames,ntimes,nbuckets,nbuckets))
    transition_boundaries=np.zeros(shape=(nnames,ntimes,nbuckets))
    for (nindex,name) in enumerate(names):
        prev_bucket=None
        for (tindex,t) in enumerate(times):
            x=data_cube[nindex,tindex,:]
            mean=x.mean()
            stddev=x.std()
            means[nindex,tindex]=mean
            stddevs[nindex,tindex]=stddev
            if stddev==0:
                continue

            skewness=np.average( (x-mean)**3 ) / stddev**3
            kurtosis=np.average( (x-mean)**4 ) / stddev**4
            minval=x.min()
            q01=np.quantile(x,0.01)
            q10=np.quantile(x,0.10)
            q50=np.median(x)
            q90=np.quantile(x,0.90)
            q99=np.quantile(x,0.99)
            maxval=x.max()
            stats_cube[nindex,tindex,0]=q10
            stats_cube[nindex,tindex,1]=q90
            stats={"t":t,"name":name,"mean":mean,"stddev":stddev,"skewness":skewness,"kurtosis":kurtosis,"minval":minval,"q01":q01,"q10":q10, "q50":q50, "q90":q90, "q99":q99, "maxval":maxval}

    # remove any metrics that have collapsed or are constant
    if False:
        for nindex in range(nnames-1,-1,-1):
            name=names[nindex]
            assert name in res["slices"]
            if any(stddevs[nindex,:]==0):
                sys.stderr.write(f"Removing {name}\n")
                slices=res["slices"]
                del slices[name]
                del names[nindex]
                stddevs=np.delete(stddevs, [nindex], axis=0)
                means=np.delete(means, [nindex], axis=0)
                # Not needed as res["names"] is the same object as names
                # res["names"].remove(name)
                data_cube=np.delete(data_cube, [nindex], axis=0)
                stats_cube=np.delete(stats_cube, [nindex], axis=0)
                transition_probs=np.delete(transition_probs,[nindex],axis=0)
                transition_boundaries=np.delete(transition_boundaries,[nindex],axis=0)
                
                nnames -= 1

                assert transition_probs.shape == (nnames,ntimes,nbuckets,nbuckets)

    assert data_cube.shape == (nnames,ntimes,nreps)
    assert len(names) == nnames
    assert stats_cube.shape == (nnames,ntimes,2)
    assert transition_probs.shape == (nnames,ntimes,nbuckets,nbuckets)
    assert transition_boundaries.shape == (nnames,ntimes,nbuckets)


    os.makedirs("data", exist_ok=True)
    with open(f"{output_dir}/{dmpci_name}-{dmpci_hash}.json", "w") as dst:
        json.dump(res,dst,indent=" ")

    tsb=TimeSeriesBundle(dmpci_name, dmpci_hash, names, times, [str(i) for i in range(nreps)], data_cube)
    tsb.get_summary_stats()
    tsb.save_as_hdf5(f"{output_dir}/{dmpci_name}-{dmpci_hash}.h5")
