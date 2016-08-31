# -*- coding: utf-8 -*-
"""
Class and methods for analysis variance and covariance

TRX team

"""

#%%

import numpy as np
import os




#%%

class Raw1D_dataset:
    """
    
    
    
    
    """
    def __init__(self, datapath=None, file_extension = '.chi', header_number = 24):
        self.datapath = datapath
        self.file_extension = file_extension
        self.header_number = header_number        
        
        
    def load_data(self):
        self.files = os.listdir(self.datapath)
        
        IQ = []
        for file in self.files:
            if file.endswith(self.file_extension):
                curve = np.genfromtxt(self.datapath+file, skip_header=self.header_number)
                IQ.append(curve[:,1])
        self.Q = curve[:,0]
        self.IQ = np.array(IQ).T 
        
        
        
        
    def sum_curves(self):
        self.summed_intensity = np.sum(self.IQ,axis=0)
        return self.summed_intensity
        
       
                
        
        
    
    