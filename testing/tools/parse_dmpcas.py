import re
from typing import *

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
