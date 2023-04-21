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

def parse_dmpcas(file_name):
    res={} # type: Dict[int,Dict[str,float]]

    with open(file_name) as src:
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

    return res
