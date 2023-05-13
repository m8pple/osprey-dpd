from .time_series_bundle import TimeSeriesBundle, StatisticId
import numpy as np
import math
from typing import *
import sys


class TimeSeriesExcessionStatistics:
    def __init__(self, bundle:TimeSeriesBundle, name_index:Union[int,str]):
        ntimes=len(bundle.times)
        nreps=bundle.data.shape[2]

        self.ntimes=ntimes
        self.nreps=nreps

        data=bundle.get_data_for_name(name_index)
        summary_stats=bundle.get_summary_stats(name_index)
        self.thresh_high=summary_stats[:,StatisticId.Q90].reshape((ntimes,1))
        self.thresh_low=summary_stats[:,StatisticId.Q10].reshape((ntimes,1))
        self.ref_stats=self.calc_statistics(data)
        self.ref_stats.sort()

    def calc_statistics(self, walks:np.ndarray):
        assert walks.shape[0]==self.ntimes
        nreps=walks.shape[1]
        excessions=(walks > self.thresh_high) | (walks<self.thresh_low)
        assert excessions.shape==(self.ntimes,nreps)
        counts=excessions.sum(axis=0)
        counts=self.ntimes-counts # We want lots of excessions to appear first, so that they correspond to low p-values
        assert len(counts.shape)==1 and counts.shape[0]==nreps
        return counts
    
    def calc_pvalues(self, walks):
        stats=self.calc_statistics(walks)
        place = np.digitize(stats, self.ref_stats, right=True)
        return place / self.nreps
    
class TimeSeriesOutlierStatistics:
    def __init__(self, bundle:TimeSeriesBundle, name_index:Union[int,str]):
        ntimes=len(bundle.times)
        nreps=bundle.data.shape[2]

        self.ntimes=ntimes
        self.nreps=nreps

        data=bundle.get_data_for_name(name_index)
        self.data=data.copy()        
        self.data.sort( axis=1 ) # Sort along nreps

        self.ref_stats=self.calc_statistics(data)
        self.ref_stats.sort()

    def calc_statistics(self, walks:np.ndarray):
        assert walks.shape[0]==self.ntimes
        nreps=walks.shape[1]

        probs=np.zeros( shape=(self.ntimes, nreps), dtype=np.float64)
        for t in range(self.ntimes):
            probs[t] = np.digitize( walks[t,:], self.data[t,:], right=True )
        probs=(probs+0.5) / (1+self.nreps) # Move into (0,1)
        assert np.all( 0<probs) and np.all(probs <1)
        probs=np.minimum(probs, 1-probs) * 2

        log_probs=np.log(probs)
        statistics=np.sum(log_probs, axis=0)

        assert len(statistics.shape)==1 and statistics.shape[0]==nreps
        return statistics
    
    def calc_pvalues(self, walks):
        stats=self.calc_statistics(walks)
        place = np.digitize(stats, self.ref_stats, right=True)
        return place / self.nreps
    

# Bucket transitions for one time-series
class TimeSeriesBucketTransitions:
    def __init__(self, bundle:TimeSeriesBundle, name_index:Union[int,str], nbuckets:int):
        nindex=bundle.get_name_index(name_index)
        ntimes=len(bundle.times)
        nreps=bundle.data.shape[2]

        self.nbuckets=nbuckets
        self.ntimes=ntimes
        self.nreps=nreps

        data=bundle.get_data_for_name(name_index)

        bucket_probs=np.zeros(shape=(ntimes,nbuckets))
        transition_probs=np.zeros(shape=(ntimes,nbuckets*nbuckets))  # encoded as prev_index*nbuckets+curr_index
        bucket_boundaries=np.zeros(shape=(ntimes,nbuckets))
        bucket_centres=np.zeros(shape=(ntimes,nbuckets))

        for i in range(nbuckets):
            bucket_boundaries[:,i]=np.quantile(data, (i+1)/nbuckets, axis=1 )
            bucket_centres[:,i]=np.quantile(data, (i+0.5)/nbuckets, axis=1)

        prev_bucket=None
        for tindex in range(ntimes):
            x=data[tindex,:]
            stddev=x.std()
            assert stddev!=0
            
            curr_bucket=np.digitize(x, bucket_boundaries[tindex,:],right=True)
            (ecdf_count,bb)=np.histogram(curr_bucket, nbuckets, range=(0, nbuckets))
            assert(bb[0]==0 and bb[-1]==nbuckets)
            assert np.sum(ecdf_count)==nreps
            prob=ecdf_count / nreps

            if prev_bucket is None:
                for prev_index in range(nbuckets):
                    transition_probs[tindex,prev_index*nbuckets:(prev_index+1)*nbuckets]=prob
            else:
                for prev_index in range(nbuckets):
                    in_prev_bucket = prev_bucket==prev_index
                    for curr_index in range(nbuckets):
                        in_curr_bucket = curr_bucket==curr_index
                        in_prev_and_curr_bucket = in_curr_bucket & in_prev_bucket
                        tprob=np.count_nonzero(in_prev_and_curr_bucket) / np.count_nonzero(in_prev_bucket)
                        transition_probs[tindex,prev_index*nbuckets+curr_index] = tprob
                            
            prev_bucket=curr_bucket

        missing=transition_probs==0
        nmissing=np.count_nonzero(missing)
        if nmissing>0:
            sys.stderr.write(f"Warning: empty transition for {bundle.names[nindex]}. Found {nmissing} out of {ntimes*nbuckets*nbuckets} ({nmissing/(ntimes*nbuckets*nbuckets)})\n")
            transition_probs[missing] = 0.5/nreps
            #TODO: The probabilities are no longer true probabilities as they don't sum to one
            # However, the statistics process doesn't actually care about that, as it just multiplies them.

        self.transition_probs=transition_probs
        self.bucket_boundaries=bucket_boundaries
        self.bucket_centres=bucket_centres

        self.ref_stats=self.calc_statistics(data)
        self.ref_stats.sort()

    def calc_buckets(self, walks:np.ndarray):
        ntimes=self.ntimes
        nbuckets=self.nbuckets

        assert walks.shape[0]==ntimes
        nreps=walks.shape[1]
        
        indices=np.zeros(shape=walks.shape,dtype=int)

        for t in range(ntimes):
            indices[t,:] = np.digitize(walks[t,:], self.bucket_boundaries[t,:], right=True)
            indices[t,:] = np.minimum(indices[t,:], nbuckets-1)
            assert (indices[t,:] < nbuckets).all(), indices[t,:]

        assert ( 0<=indices ).all() and (indices < nbuckets).all()

        return indices

    def calc_statistics(self, walks:np.ndarray):
        ntimes=self.ntimes
        nbuckets=self.nbuckets
        nreps=walks.shape[1]

        indices=self.calc_buckets(walks)
        
        hop_indices = indices[:-1,:] * nbuckets + indices[1:,:] #type: np.ndarray
        assert hop_indices.shape == (ntimes-1,nreps)

        walk_probs = np.zeros(shape=(ntimes-1,nreps))
        for t in range(ntimes-1):
            pp=self.transition_probs[t+1,:][hop_indices[t]]
            walk_probs[t,:] = pp

        #assert not np.any(walk_probs==0)

        log_probs = np.log(walk_probs)
        res=np.sum(log_probs, axis=0)

        unique = len(set(res)) / len(res)

        return res
    
    def calc_pvalues(self, walks):
        stats=self.calc_statistics(walks)
        place = np.digitize(stats, self.ref_stats, right=True)
        return place / self.nreps