import awkward as ak
import numpy as np
from coffea.nanoevents import NanoEventsFactory
from ecalphisym import EcalPhiSymSchema
from coffea.processor import AccumulatorABC, ProcessorABC, Runner, FuturesExecutor


class EcalPhiSymAccumulator(AccumulatorABC):
    def __init__(self, ebhits=None, eehits=None, info=None, run=None):
        self._ebhits = ebhits
        self._eehits = eehits
        self._info = info
        self._run = run 

    def add(self, other):
        if self._ebhits is not None or self._eehits is not None or self._info is not None:
            self._ebhits = ak.concatenate([self._ebhits, other._ebhits], behavior=self._ebhits.behavior)
            self._eehits = ak.concatenate([self._eehits, other._eehits], behavior=self._eehits.behavior)
            self._info = ak.concatenate([self._info, other._info], behavior=self._info.behavior)  
            self._run = ak.concatenate([self._run, other._run])  
        else:
            self._ebhits = other._ebhits
            self._eehits = other._eehits
            self._info = other._info            
            self._run = other._run        
            
    def identity(self):
        return EcalPhiSymAccumulator()
    
    @property
    def ebhits(self):
        """Return the ebhits array"""
        return self._ebhits

    @property
    def eehits(self):
        """Return the eehits array"""
        return self._eehits
    
    @property
    def info(self):
        """Return the info array"""
        return self._info

    @property
    def run(self):
        """Return the runs array"""
        return self._run
    

class Processor(ProcessorABC):   
    def process(self, runs):
        unique_runs = np.unique(runs.run, return_index=True)[0]
        counts = np.unique(runs.run, return_index=True)[1]
        splits = np.diff(np.concatenate([counts, [len(runs.run)]]))
        info = ak.unflatten(runs.EcalPhiSymInfo, splits, axis=0, behavior=runs.behavior).sum(axis=1)
        ebhits = ak.unflatten(runs.EcalPhiSymEB, splits, axis=0, behavior=runs.behavior).sum(axis=1)
        eehits = ak.unflatten(runs.EcalPhiSymEE, splits, axis=0, behavior=runs.behavior).sum(axis=1)
        return EcalPhiSymAccumulator(ebhits=ebhits, eehits=eehits, info=info, run = unique_runs )
    
    def postprocess(self, accumulator):
        return accumulator

