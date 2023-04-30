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

def to_colour(p,nbuckets):
    e=1.0/nbuckets
    if p < e:
        r = p / e
        return (r,0,0)
    if p >= e:
        r = ( p - e ) / (1-e)
        return (0,r,0)
    
def plot_transitions_heat_v2(tsb:TimeSeriesBundle, name:str):
    nbuckets=5
    bb=TimeSeriesBucketTransitions(tsb, name, 5)

    (fix,axes)=plt.subplots(nbuckets+1,1,sharex=True)

    axes[0].plot(bb.bucket_centres)

    for from_bucket in range(nbuckets):
        probs=bb.transition_probs[:,from_bucket*nbuckets:(from_bucket+1)*nbuckets]
        axes[from_bucket+1].imshow(probs.transpose(), vmin=0, vmax=1)

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
            plot_transitions_heat_v2(tsb,name)
            pdf.savefig()
            plt.close()
        
        plt.figure(figsize=(8, 6))
        plot_summary_stats(tsb, name, walks)
        pdf.savefig()
        plt.close()

