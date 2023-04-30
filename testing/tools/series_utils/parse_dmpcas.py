import re
import os
import sys
from typing import *
import numpy as np
from series_utils.time_series_bundle import TimeSeriesBundle
from run_replicates import calc_dmpci_file_hash, make_dir_name

dmpcas_scalar_names=[
    "Temperature",
    "Pressure",
    "[a-zA-Z]+ EE distance",
    "Bond Length",
    "[a-zA-Z]+ bond length"
]
dmpcas_scalar_names_re="|".join(dmpcas_scalar_names)

dmpcas_vector_names=[
    "CM Mom",
    "CM Pos"
]
dmpcas_vector_names_re="|".join(dmpcas_vector_names)

def parse_dmpci_types(directory, dmpci_name):
    bead_types=list()
    polymer_types=list()

    with open(f"{directory}/dmpci.{dmpci_name}") as src:
        for line in src:
            line=line.strip()
            parts=line.split(" ")
            if parts[0]=="Bead":
                bead_types.append(parts[1])
            if parts[0]=="Polymer":
                polymer_types.append(parts[1])

    return (bead_types,polymer_types)

def parse_and_add_dmpchs(directory, dmpci_name, res):
    (bead_types,polymer_types)=parse_dmpci_types(directory, dmpci_name)

    times=res.keys()

    with open(f"{directory}/dmpchs.{dmpci_name}") as src:
        for line in src:
            line=line.strip()
            while True:
                pline=line
                line=pline.replace("  "," ")
                if line==pline:
                    break

            cols=line.split(" ")

            time=int(cols[0])
            if time not in times:
                continue

            for bt in range(len(bead_types)):
                v=float(cols[6+bt].strip())
                res[time][f"Bead {bt} drift"]=v

def parse_dmpcas(directory, dmpci_name):
    res={}

    with open(f"{directory}/dmpcas.{dmpci_name}") as src:
        for line in src:
            line=line.strip()

            if line=="":
                continue

            m=re.match("Time = ([0-9]+)", line)
            if m:
                time=int(m.group(1))
                res[time]={}
                continue
            assert time is not None


            if re.match(dmpcas_scalar_names_re,line):
                values=next(src) # TODO : is this legal, while iterating over it?
                values=values.split()
                assert len(values)==2
                mean=float(values[0])
                
                res[time][line]=mean

            if re.match(dmpcas_vector_names_re,line):
                _ = next(src) # Skip x y and z
                _ = next(src)
                _ = next(src)
                values = next(src) # Get magnitude
                values=values.split()
                assert len(values)==2
                mean=float(values[0])
                
                res[time][line]=mean

    parse_and_add_dmpchs(directory, dmpci_name, res)

    return res

def parse_dmpcas_to_bundle(directory, dmpci_name, dmpci_hash):
    all_slices={}
    r=0
    while True:
        replicate_dir=f"{directory}/{make_dir_name(r)}"
        if not os.path.isdir(replicate_dir):
            sys.stderr.write(f"Missing directory {replicate_dir} finished search\n")
            break

        assert os.path.isfile(f"{replicate_dir}/dmpci.{dmpci_name}-{dmpci_hash}"), f"Missing {replicate_dir}/dmpci.{dmpci_name}-{dmpci_hash}"
        
        if not os.path.exists(f"{replicate_dir}/finished.{dmpci_name}-{dmpci_hash}"):
            sys.stderr.write(f"Skipping {replicate_dir}\n")
        else:
            slice=parse_dmpcas(replicate_dir, f"{dmpci_name}-{dmpci_hash}")

            all_slices[r]=slice
            got_times=list(slice.keys())
            got_times.sort()
            for (k,v) in slice.items():
                got_names=list(v.keys())
            got_names.sort()
            if len(all_slices)==1:
                times=got_times
                names=got_names
            else:
                assert times==got_times
                assert names==got_names
        r=r+1

    assert len(all_slices)>0

    reps=list(all_slices.keys())
    reps.sort()
    nreps=len(reps)


    cube=np.zeros(shape=(len(names), len(times), nreps))
    for (rindex,rep) in enumerate(reps):
        data=all_slices[rep]
        for (tindex,time) in enumerate(times):
            slice=data[time]
            for (name,value) in slice.items():
                assert name in names, f"{names}, {name}"
                nindex=names.index(name)
                cube[nindex,tindex,rindex]=value

    res=TimeSeriesBundle(
        dmpci_name=dmpci_name,
        dmpci_hash=dmpci_hash,
        names=names,
        rep_names=reps,
        times=times,
        data=cube
    )
    return res