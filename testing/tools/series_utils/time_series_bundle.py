import numpy as np
import scipy
from typing import * 
from enum import IntEnum
import h5py

class StatisticId(IntEnum):
    MIN = 0
    Q01 = 1
    Q10 = 2
    Q50 = 3
    Q90 = 4
    Q99 = 5
    MAX = 6
    MEAN = 7
    STDDEV = 8
    SKEWNESS = 9
    KURTOSIS = 10

class TimeSeriesBundle:

    def __init__(self, dmpci_name:str, dmpci_hash:str, names:List[str], times:List[int], rep_names:List[str], data : Optional[np.ndarray] = None, stats_cube : Optional[np.ndarray] = None):
        assert len(names)>0
        
        self.dmpci_name=dmpci_name
        self.dmpci_hash=dmpci_hash
        self.names=np.array(names, dtype=str)
        assert len(set(rep_names)) == len(rep_names), rep_names
        self.rep_names=np.array(rep_names, dtype=str)
        self.nreps=len(self.rep_names)
        self.times=np.array(times, dtype=float)
        self.ntimes=len(self.times)
        self.data = data
        if self.data is None:
            self.data = np.zeros(shape=(len(names), len(times), len(rep_names)), dtype=float)
        else:
            self.data = np.array(data)

        if stats_cube is not None:
            assert stats_cube.shape == (len(names),len(times),len(StatisticId))
        self.stats_cube=stats_cube

    def save_as_hdf5(self, file_name):
        with h5py.File(file_name, 'w') as dest:
            dest.attrs["dmpci_name"] = self.dmpci_name
            dest.attrs["dmpci_hash"] = self.dmpci_hash
            dest.create_dataset("names", data= [n.encode("ascii","ignore") for n in self.names])
            dest.create_dataset("times", data=self.times)
            dest.create_dataset("rep_names", data=[n.encode("ascii", "ignore") for n in self.rep_names])
            dest.create_dataset("data", data=self.data)

            #if self.stats_cube is not None:
            #    dest["stats_cube"] = self.stats_cube

    @staticmethod
    def load_from_hdf5(file_name):
        with h5py.File(file_name, "r") as src:
            stats_cube=None
            if "stats_cube" in src:
                stats_cube=src["stats_cube"]

            return TimeSeriesBundle(
                dmpci_name = src.attrs["dmpci_name"],
                dmpci_hash = src.attrs["dmpci_hash"],
                names = [n.decode("ascii") for n in src["names"]],
                times = src.get("times"),
                rep_names = [n.decode("ascii") for n in src["rep_names"]],
                data=src["data"]
            )

    def set_data_point(self, nindex, tindex, rindex, value):
        self.data[nindex, tindex, rindex] = value

    def get_name_index(self, index:Union[int,str]):
        if isinstance(index, str):
            index=list(self.names).index(index)
        return index
    
    def get_data_for_name(self, index:Union[int,str]):
        return self.data[self.get_name_index(index), :, :]
    
    def get_rep_index(self, index:Union[int,str]):
        if isinstance(index, str):
            index=self.rep_names.index(index)
        return index

    def remove_series(self, nindex):
        self.names.pop(nindex)
        np.erase( self.names, [nindex], axis=0) 
        if self.stats_cube:
            np.erase(self.stats_cube, [nindex], axis=0)
    
    def get_summary_stats(self, nindex=None):
        if nindex is not None:
            return self.get_summary_stats()[self.get_name_index(nindex),:,:]
        
        if self.stats_cube is not None:
            return self.stats_cube

        self.stats_cube=np.zeros(shape=(len(self.names), len(self.times), 11))

        self.stats_cube[:,:,StatisticId.MIN] = np.min(self.data, axis = 2)
        self.stats_cube[:,:,StatisticId.Q01] = np.quantile(self.data, 0.01, axis = 2)
        self.stats_cube[:,:,StatisticId.Q10] = np.quantile(self.data, 0.10, axis = 2)
        self.stats_cube[:,:,StatisticId.Q50] = np.quantile(self.data, 0.50, axis = 2)
        self.stats_cube[:,:,StatisticId.Q90] = np.quantile(self.data, 0.90, axis = 2)
        self.stats_cube[:,:,StatisticId.Q99] = np.quantile(self.data, 0.99, axis = 2)
        self.stats_cube[:,:,StatisticId.MAX] = np.max(self.data, axis = 2)
        self.stats_cube[:,:,StatisticId.MEAN] = np.mean(self.data, axis = 2)
        self.stats_cube[:,:,StatisticId.STDDEV] = np.std(self.data, axis = 2)
        self.stats_cube[:,:,StatisticId.SKEWNESS] = scipy.stats.skew(self.data, axis = 2)
        self.stats_cube[:,:,StatisticId.KURTOSIS] = scipy.stats.kurtosis(self.data, axis = 2)

        return self.stats_cube