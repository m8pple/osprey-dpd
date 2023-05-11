#!/usr/bin/env python3
import json
import sys
import re
from typing import *
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from run_replicates import hash_text
from series_utils.time_series_bundle import TimeSeriesBundle, StatisticId
from series_utils.walk_statistics import TimeSeriesBucketTransitions

def plot_summary_stats(tsb:TimeSeriesBundle, name, walks:Optional[TimeSeriesBundle]):
    stats=tsb.get_summary_stats(name)
    sel=[StatisticId.MIN,StatisticId.Q01,StatisticId.Q10,StatisticId.Q50,StatisticId.Q90,StatisticId.Q99,StatisticId.MAX]
    alpha=1.0
    style="-"
    if walks:
        alpha=0.5
        style=":"
    for s in sel:
        plt.plot(tsb.times, stats[:,s.value], style, label=s.name, alpha=alpha)
    if walks:
        data=walks.get_data_for_name(name)
        for i in range(0,data.shape[1]):
            plt.plot(tsb.times, data[:,i], label=f"Repeat {i}")
    
    plt.legend()
    plt.title(name)

def plot_walk_quantiles(tsb:TimeSeriesBundle, name, walks:TimeSeriesBundle):
    ref_data=tsb.get_data_for_name(name)
    walk_data=walks.get_data_for_name(name)

    ref_data_sorted=ref_data.copy()
    ref_data_sorted.sort(axis=1)

    walk_indices=np.zeros(shape=walk_data.shape, dtype=np.float64)
    for t in range(tsb.ntimes):
        walk_indices[t,:]=np.digitize(walk_data[t,:],ref_data_sorted[t,:], right=True)
    walk_quantiles=walk_indices / tsb.nreps
    
    for i in range(0, walk_quantiles.shape[1]):
        plt.plot(tsb.times, walk_quantiles[:,i], label=f"Repeat {i}")
    plt.xlabel(f"Empirical quantile (from {tsb.nreps} reference samples)")
    plt.legend()
    plt.title(name) 

def to_colour(p,nbuckets):
    e=1.0/nbuckets
    if p < e:
        r = p / e
        return (r,0,0)
    if p >= e:
        r = ( p - e ) / (1-e)
        return (0,r,0)
    
def plot_transitions_heat_v2(tsb:TimeSeriesBundle, name:str, walk:Optional[TimeSeriesBundle] = None):
    nbuckets=5
    bb=TimeSeriesBucketTransitions(tsb, name, 5)

    nsubplots=nbuckets+1
    if walk:
        nsubplots+=1

    (fix,axes)=plt.subplots(nsubplots,1,sharex=True)

    axes[0].plot(bb.bucket_centres)

    for from_bucket in range(nbuckets):
        probs=bb.transition_probs[:,from_bucket*nbuckets:(from_bucket+1)*nbuckets]
        axes[from_bucket+1].imshow(probs.transpose(), vmin=0, vmax=1)

    if walk:
        aa=axes[nbuckets+1]
        data=walk.get_data_for_name(name)
        indices=bb.calc_buckets(data)
        for r in range(data.shape[1]):
            aa.plot(indices[:,r])

    plt.title(name)


walks=None

dest_pdf_file=None

spec=sys.argv[1]
if spec.startswith("dmpci."):
    dmpci_file=spec
    with open(dmpci_file,"r") as src:
        dmpci_contents=src.read()
    dmpci_hash=hash_text(dmpci_contents)
    m=re.match("dmpci[.]([-_a-zA-Z0-9]+)", dmpci_file)
    assert m, "Couldn't parse '{dmpci_file}' into 'dmpci.NAME' "
    dmpci_name=m.group(1)
    sys.stderr.write(f"dmpci {dmpci_name}, hash {dmpci_hash}\n")

    tsb=TimeSeriesBundle.load_from_hdf5(f"data/{dmpci_name}-{dmpci_hash}.h5")

    dest_pdf_file=f"data/{dmpci_name}-{dmpci_hash}.pdf"
elif spec.endswith(".h5"):
    data_file=spec
    walks=TimeSeriesBundle.load_from_hdf5(data_file)
    dmpci_name=walks.dmpci_name
    dmpci_hash=walks.dmpci_hash
    
    tsb=TimeSeriesBundle.load_from_hdf5(f"data/{dmpci_name}-{dmpci_hash}.h5")
    assert tsb.dmpci_hash==walks.dmpci_hash

    dest_pdf_file=data_file.replace(".h5",".pdf")
else:
    sys.stderr.write("Didn't udnerstand source spec.")
    sys.exit(1)


names=tsb.names
times=tsb.times
reps=tsb.rep_names
nreps=len(reps)

with PdfPages(dest_pdf_file) as pdf:
    for name in names:
        tsb_stats=tsb.get_summary_stats(name)
        if np.any(tsb_stats[:,StatisticId.STDDEV]==0):
            continue

        if True:
            plt.figure(figsize=(8, 6))
            plot_transitions_heat_v2(tsb,name, walks)
            pdf.savefig()
            plt.close()
        
        plt.figure(figsize=(8, 6))
        plot_summary_stats(tsb, name, walks)
        pdf.savefig()
        plt.close()

        if walks:
            plt.figure(figsize=(8,6))
            plot_walk_quantiles(tsb, name, walks)
            pdf.savefig()
            plt.close()
