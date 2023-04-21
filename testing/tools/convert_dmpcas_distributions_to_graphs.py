#!/usr/bin/env python3
import json
import sys
import re
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from run_replicates import hash_text

def plot_slice(name):
    print(name)
    slice=data["slices"][name]["times"]
    x=np.zeros(shape=(len(times),5))
    for (tindex,time) in enumerate(times):
        print(time)
        stats=slice[tindex]
        assert stats["t"]==time
        x[tindex,:]=( stats["minval"],stats["q10"],stats["q50"],stats["q90"],stats["maxval"] )
    plt.plot(x[:,0],"-",label="q1")
    plt.plot(x[:,1],"-",label="q10")
    plt.plot(x[:,2],"-",label="median")
    plt.plot(x[:,3],"-",label="q90")
    plt.plot(x[:,4],"-",label="q99")
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

def plot_transitions(name):
    slice=data["slices"][name]["times"]
    nbuckets=len(slice[0]["ecdf"]["buckets"])

    buckets=np.zeros(shape=(len(times),1+nbuckets))
    for (tindex,time) in enumerate(times):
        stats=slice[tindex]
        assert stats["t"]==time
        ecdf=stats["ecdf"]
        buckets[tindex,0]=stats["minval"]
        buckets[tindex,1:]=ecdf["buckets"]
    
    centres=(buckets[:,:-1]+buckets[:,1:])/2
    for (tindex,time) in enumerate(times):
        if tindex == 0:
            continue
        stats=slice[tindex]
        assert stats["t"]==time
        ecdf=stats["ecdf"]
        xs=[ times[tindex-1],times[tindex] ]
        for fromb in range(nbuckets):
            print(ecdf)
            transition_probs=ecdf["transition_probs"][fromb]
            for tob in range(nbuckets):                
                #ys=[centres[tindex-1,fromb],centres[tindex,tob]]
                ys=[fromb,tob]
                plt.plot(xs,ys,alpha=0.5,color=to_colour(transition_probs[tob],nbuckets))

            
    #for b in range(nbuckets):
    #    plt.scatter(times, centres[:,b], label=f"bucket {b}" )


    plt.legend()
    plt.title(name)

def get_transition_probs(slice,from_bucket):
    nbuckets=len(slice[0]["ecdf"]["buckets"])

    probs=np.zeros(shape=(len(times),nbuckets))
    for (tindex,time) in enumerate(times):
        if tindex==0:
            probs[tindex,:]=1/nbuckets
        else:
            stats=slice[tindex]
            assert stats["t"]==time
            ecdf=stats["ecdf"]
            probs[tindex,:]=ecdf["transition_probs"][from_bucket]
    return probs

def plot_transitions_heat(name):
    slice=data["slices"][name]["times"]
    nbuckets=len(slice[0]["ecdf"]["buckets"])

    (fix,axes)=plt.subplots(nbuckets+1,1,sharex=True)

    buckets=np.zeros(shape=(len(times),1+nbuckets))
    for (tindex,time) in enumerate(times):
        stats=slice[tindex]
        assert stats["t"]==time
        ecdf=stats["ecdf"]
        buckets[tindex,0]=stats["minval"]
        buckets[tindex,1:]=ecdf["buckets"]
    centres=(buckets[:,:-1]+buckets[:,1:])/2
    axes[0].plot(centres)

    for from_bucket in range(nbuckets):
        probs=get_transition_probs(slice,from_bucket)
        axes[from_bucket+1].imshow(probs.transpose(),vmin=0, vmax=1)

    plt.title(name)


dmpci_file=sys.argv[1]
with open(dmpci_file,"r") as src:
    dmpci_contents=src.read()
dmpci_hash=hash_text(dmpci_contents)
m=re.match("dmpci[.]([-_a-zA-Z0-9]+)", dmpci_file)
assert m, "Couldn't parse '{dmpci_file}' into 'dmpci.NAME' "
dmpci_name=m.group(1)
sys.stderr.write(f"dmpci {dmpci_name}, hash {dmpci_hash}\n")


with open(f"data/{dmpci_name}-{dmpci_hash}.json", 'r') as src:
    data=json.load(src )

names=data["names"]
times=data["times"]
nreps=data["nreps"]

with PdfPages(f"data/{dmpci_name}-{dmpci_hash}.pdf") as pdf:
    for name in names:
        plt.figure(figsize=(8, 6))
        plot_transitions_heat(name)
        pdf.savefig()
        plt.close()
        
        plt.figure(figsize=(8, 6))
        plot_slice(name)
        pdf.savefig()
        plt.close()

