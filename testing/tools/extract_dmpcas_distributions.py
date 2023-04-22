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

from parse_dmpcas import parse_dmpcas
from run_replicates import make_dir_name, hash_text

def calculate_quantised_walk_statistic(nbuckets:int, ntimes:int, nreps:int, boundaries : np.ndarray, transition_probs : np.ndarray, paths : np.ndarray):
    assert transition_probs.shape==(ntimes,nbuckets,nbuckets)    
    reshaped_probs=np.zeros(shape=(ntimes,nbuckets*nbuckets))
    for b in range(nbuckets):
        reshaped_probs[:,b*nbuckets:(b+1)*nbuckets]=nbuckets / transition_probs[:,b,:]
    
    assert paths.shape==(ntimes,nreps) or paths.shape==(ntimes,)
    sys.stderr.write(f"paths={paths.shape}, boundaries={boundaries.shape}\n")
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
    sys.stderr.write(str(log_probs.shape)+"\n")
    res=np.sum(log_probs, axis=0)
    sys.stderr.write(str(res.shape)+"\n")
    return res

def calculate_extreme_value_statistic(ntimes:int, nreps:int, ul_boundaries : np.ndarray, paths : np.ndarray):
    assert ul_boundaries.shape == (ntimes,2)
    assert paths.shape==(ntimes,nreps)
    assert np.sum(ul_boundaries[:,0] < ul_boundaries[:,1]) > ntimes/2, str(ul_boundaries)+"\n"+str(paths)
    extrema = (paths < ul_boundaries[:,0].reshape((ntimes,1))) | (ul_boundaries[:,1].reshape((ntimes,1)) < paths)
    statistics=np.sum(extrema, axis=0) / ntimes
    assert statistics.shape==(nreps,)
    return statistics


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
    working_dir=f"./working/{dmpci_name}-{dmpci_hash}"
    os.makedirs(working_dir,exist_ok=True)
else:
    working_dir=sys.argv[2]

sys.stderr.write(f"Looking for results for {dmpci_name}-{dmpci_hash} in directory {working_dir}\n")

all_slices=[]

r=0
while True:
    replicate_dir=f"{working_dir}/{make_dir_name(r)}"
    if not os.path.isdir(replicate_dir):
        sys.stderr.write(f"Missing directory {replicate_dir} finished search\n")
        break

    assert os.path.isfile(f"{replicate_dir}/dmpci.{dmpci_name}"), f"Missing {replicate_dir}/dmpci.{dmpci_name}"
    
    if not os.path.exists(f"{replicate_dir}/finished"):
        sys.stderr.write("Skipping {replicate_dir}\n")
    else:
        slice=parse_dmpcas(replicate_dir, dmpci_name)
        all_slices.append(slice)
    r=r+1

times=[t for t in all_slices[0].keys()]
names=list(all_slices[0][times[0]].keys())
reps=range(0,len(all_slices))

nnames=len(names)
nreps=len(reps)
ntimes=len(times)

res={}

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

sys.stderr.write(str(res["slices"].keys())+"\n")

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
        
        ecdf_boundaries=[ np.quantile(x,(i+1)/nbuckets) for i in range(nbuckets) ]
        ecdfs={"buckets":ecdf_boundaries}
        ecdf_boundaries=np.array(ecdf_boundaries)
        transition_boundaries[nindex,tindex,:]=ecdf_boundaries
        
        curr_bucket=np.digitize(x, ecdf_boundaries,right=True)
        (ecdf_count,bb)=np.histogram(curr_bucket, nbuckets, range=(0, nbuckets))
        assert(bb[0]==0 and bb[-1]==nbuckets)
        assert np.sum(ecdf_count)==nreps
        prob=ecdf_count / nreps
        ecdfs["probs"]=list(prob)


        if prev_bucket is None:
            for prev_index in range(nbuckets):
                transition_probs[nindex,tindex,prev_index,:]=prob
        else:
            rows=[]
            for prev_index in range(nbuckets):
                in_prev_bucket = prev_bucket==prev_index
                cols=[]
                for curr_index in range(nbuckets):
                    in_curr_bucket = curr_bucket==curr_index
                    in_prev_and_curr_bucket = in_curr_bucket & in_prev_bucket
                    prob=np.count_nonzero(in_prev_and_curr_bucket) / np.count_nonzero(in_prev_bucket)
                    cols.append(prob)
                    transition_probs[nindex,tindex,prev_index,curr_index]=prob
                assert( abs( sum(cols) - 1 ) < 0.000001 )
                rows.append(cols)
            ecdfs["transition_probs"]=rows

        prev_bucket=curr_bucket
        
        stats["ecdf"]=ecdfs
        res["slices"][name]["times"].append(stats)
        #print(f"{t:3}, {name:11}, {len(x)}, {mean:.8}, {stddev:.6}, {skewness:.5}, {kurtosis:.5}, {minval:.8}, {q01:.8}, {q10:.8}, {q50:.8}, {q90:.8}, {q99:.8}, {maxval:.8}")

assert transition_probs.shape == (nnames,ntimes,nbuckets,nbuckets)

# remove any metrics that have collapsed or are constant
for nindex in range(nnames-1,-1,-1):
    name=names[nindex]
    assert name in res["slices"]
    if any(stddevs[nindex,:]==0):
        sys.stderr.write(f"Removing {name}, means={means[nindex,:]} stddev={stddevs[nindex,:]}\n")
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

for (nindex,name) in enumerate(names):
    s=calculate_quantised_walk_statistic(nbuckets, ntimes, nreps, transition_boundaries[nindex,:,:], transition_probs[nindex,:,:,:], data_cube[nindex,:,:])
    s.sort()
    sys.stderr.write(str(s.shape)+"\n")
    sys.stderr.write(name+"\n")
    res["slices"][name]["walk_stat_ecdf"]=[float(x) for x in s]

    ul_boundaries=stats_cube[nindex,:,:]
    s=calculate_extreme_value_statistic(ntimes, nreps, ul_boundaries, data_cube[nindex,:,:])
    s.sort()
    sys.stderr.write(str(s))

sys.stderr.write(str(res["slices"].keys())+"\n")


os.makedirs("data", exist_ok=True)
with open(f"data/{dmpci_name}-{dmpci_hash}.json", "w") as dst:
    json.dump(res,dst,indent=" ")
