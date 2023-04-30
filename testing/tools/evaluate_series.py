#!/usr/bin/env python3

import os
import sys
import json
import numpy as np
from series_utils.parse_dmpcas import parse_dmpcas_to_bundle, calc_dmpci_file_hash
from series_utils.walk_statistics import TimeSeriesBucketTransitions, TimeSeriesExcessionStatistics, TimeSeriesOutlierStatistics
from series_utils.time_series_bundle import TimeSeriesBundle, StatisticId

results_dir = sys.argv[1]
dmpci_file = sys.argv[2]
assert os.path.exists(dmpci_file)

dmpci_name = dmpci_file.removeprefix("dmpci.")

dmpci_hash = calc_dmpci_file_hash(".", dmpci_name)

sys.stderr.write(f"Loading data/{dmpci_name}-{dmpci_hash}.h5\n")
ref_data=TimeSeriesBundle.load_from_hdf5(f"data/{dmpci_name}-{dmpci_hash}.h5")

series=parse_dmpcas_to_bundle(results_dir, f"{ref_data.dmpci_name}", dmpci_hash=dmpci_hash)


"""
We run N reps, each of which is expected to be IID uniform

We calculate p'=min(p[i], i=1..N)

We want a thresh such that P[ p' < thresh' ] = thresh

1-(1-thresh')**N = thresh
1-thresh = (1-thresh')**N
(1-thresh)**(1/N)=1-thresh'
thresh'=1-(1-thresh)**(1/N)

"""
def to_message_and_code(p):
    assert(len(p.shape)==1)
    reps=p.shape[0]
    threshfail=1-(1-0.01)**(1.0/reps)
    threshwarn=1-(1-0.05)**(1.0/reps)
    sys.stderr.write(f"threshwarn={threshwarn}\n")
    res=np.ndarray(shape=p.shape, dtype=object)
    code=np.zeros(shape=p.shape)
    res[:]="FAIL"
    code[p<threshfail]=1
    res[p>threshfail] = "warn"
    res[p>threshwarn] = "pass"
    return (res,code)

if ref_data.nreps < 10*series.nreps:
    sys.stderr.write(f"Warning : ref data has {ref_data.nreps} reps, series has {series.nreps}. Stats will be dubious.")

ref_stats=ref_data.get_summary_stats()
fails=0
for (nindex,name) in enumerate(ref_data.names):
    sys.stderr.write(f"{nindex}, {name}\n")
    if np.any(ref_stats[nindex,:,StatisticId.STDDEV]==0):
        sys.stderr.write("  Skipping\n")
        continue

    data=series.get_data_for_name(name)

    tt=TimeSeriesBucketTransitions(ref_data, nindex, 5)
    pvalues=tt.calc_pvalues(data)
    (message,code)=to_message_and_code(pvalues)
    for (rindex,pvalue) in enumerate(pvalues):
        print(f"{ref_data.dmpci_name}, {nindex}, {name}, {rindex}, bucket_walk, {message[rindex]}, {pvalue:.5f}")
    fails+=np.sum(code)

    if False:
        tt=TimeSeriesExcessionStatistics(ref_data, nindex)
        pvalues=tt.calc_pvalues(data)
        (message,code)=to_message_and_code(pvalues)
        for (rindex,pvalue) in enumerate(pvalues):
            print(f"{ref_data.dmpci_name}, {nindex}, {name}, {rindex}, excessions, {message[rindex]}, {pvalue:.5f}")
        fails+=np.sum(code)

    tt=TimeSeriesOutlierStatistics(ref_data, nindex)
    pvalues=tt.calc_pvalues(data)
    (message,code)=to_message_and_code(pvalues)
    for (rindex,pvalue) in enumerate(pvalues):
        print(f"{ref_data.dmpci_name}, {nindex}, {name}, {rindex}, outliers, {message[rindex]}, {pvalue:.5f}")
    fails+=np.sum(code)


if ref_data.nreps < 10*series.nreps:
    sys.stderr.write(f"Warning : ref data has {ref_data.nreps} reps, series has {series.nreps}. Stats will be dubious.")

if fails > 0:
    sys.exit(1)
