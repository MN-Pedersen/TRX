# -*- coding: utf-8 -*-
"""
Class and methods for analysis variance and covariance

TRX team

"""

#%%

import numpy as np
import os
import re
import matplotlib.pyplot as plt


#%%

class Raw1D_dataset:
    """
    
    
    
    
    """
    def __init__(self, datapath=None, file_extension='.chi', header_number=24,
                 sort=True):
        self.datapath = datapath
        self.file_extension = file_extension
        self.header_number = header_number
        self.sort = sort        
        
        
    def load_data(self, separator = '\\'):
        self.files = os.listdir(self.datapath)
        
        IQ = []
        file_number = []
        for file in self.files:
            if file.endswith(self.file_extension):
                curve = np.genfromtxt(self.datapath+separator+file, skip_header=self.header_number)
                IQ.append(curve[:,1])
                if self.sort:
                    file_number.append(int(re.findall('\d{4}',file)[-1]))
                    
        if self.sort:
            sort_index = np.argsort(file_number)
            IQ = np.array(IQ)[sort_index,:]    
                
                
        self.Q = curve[:,0]
        self.IQ = np.array(IQ).T
        
    
    def sum_curves(self, clean=False, clean_thresshold=1e5, return_vectors=False):
        self.summed_intensity = np.sum(self.IQ,axis=0)
        if clean:
            index_remove = self.summed_intensity < clean_thresshold
            self.IQ = self.IQ[:,~index_remove]
            self.summed_intensity = np.sum(self.IQ,axis=0)
            
        if return_vectors:
            return self.summed_intensity
            
            
    def intensity_trends(self, tuples_list, produce_plot=True):       
        max_int = np.max(IQ)*1.05

        fig = plt.figure()
        plt.plot(self.Q, self.IQ)
        for tuples in tuples_list:
            plt.plot([tuples[0], tuples[0]], [0, max_int], '-b')
            plt.plot([tuples[1], tuples[1]], [0, max_int], '-b')

        intensities = []
        for tuples in tuples_list:
            Q_index = (Q >= tuples[0]) & (Q >= tuples[1])
            intensities.append(np.sum(IQ[Q_index,:], axis=0))
        
        self.intensities = intensities 


        if produce_plot:
            fig = plt.figure()
            num_tuples = len(tuples_list)
            for num in range(num_tuples):   
                for num2 in range(num_tuples-1):
                    if num > num2:
                        ax = plt.subplot2grid((num_tuples-1, num_tuples-1), (num2, num-1))
                        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                        plt.scatter(intensities[num], intensities[num2], s=20, facecolors='none', edgecolors='k')
                        plt.xlabel('Q = %s' % str(tuples_list[num]))
                        plt.ylabel('Q = %s' % str(tuples_list[num2]))
            #plt.tight_layout()
    

    def svd_analysis(self, norm_range=[8,9], num_components=4, curve_region=[0,1], 
                     produce_plot=True, return_vectors=False):
        data = self.IQ[:, curve_region[0]:curve_region[1]]
        Q_index = (self.Q >= norm_range[0]) & (self.Q <= norm_range[1])
        norm_curves = np.array([data[:,num]/np.trapz(data[Q_index,num], 
                                x=self.Q[Q_index]) for num in range(np.shape(data)[1])])
        
        norm_curves = np.array(norm_curves).T
        U, s, V = np.linalg.svd(norm_curves)
        
        if produce_plot:
            fig = plt.figure()
            plt.loglog(s, 'o')
        
            fig = plt.figure()
            for num in range(num_components):
                plt.plot(self.Q, U[:,num])
            
            fig = plt.figure()
            for num in range(num_components):
                plt.plot(V[num,:])
                print(np.linalg.norm(V[num,:]))
                
        if return_vectors:
            return norm_curves, U, s, V

        
        



    
    def return_files(self):
        return self.files
    
    def return_data(self):
        return self.Q, self.IQ    
       
                
        
        
    
    