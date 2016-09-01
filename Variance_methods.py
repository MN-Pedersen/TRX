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
import matplotlib.cm as cm


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
        
    
    def sum_curves(self, clean=False, clean_thresshold=1e7, keep_tail=True,
                   tail_intensity=5e4):
        self.summed_intensity = np.sum(self.IQ,axis=0)
        
        if keep_tail:
            index_keep = self.summed_intensity < tail_intensity
            self.tail = self.IQ[:,index_keep]
        
        
        if clean:
            index_remove = self.summed_intensity < clean_thresshold
            self.IQ = self.IQ[:,~index_remove]
            self.summed_intensity = np.sum(self.IQ,axis=0)
            
            
    def intensity_trends(self, tuples_list, produce_plot=True, use_subset=False):       
        max_int = np.max(self.IQ)*1.05

        fig = plt.figure()
        plt.plot(self.Q, self.IQ)
        for tuples in tuples_list:
            plt.plot([tuples[0], tuples[0]], [0, max_int], '-b')
            plt.plot([tuples[1], tuples[1]], [0, max_int], '-b')

        intensities = []
        for tuples in tuples_list:
            Q_index = (self.Q >= tuples[0]) & (self.Q >= tuples[1])
            if use_subset:
                intensities.append(np.sum(self.subset[Q_index,:], axis=0))
            else:    
                intensities.append(np.sum(self.IQ[Q_index,:], axis=0))
        
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
                    produce_plot=True, use_subset=False):
                        
        if use_subset:
            data = self.subset
        else:
            data = self.IQ[:, curve_region[0]:curve_region[1]]
            
        Q_index = (self.Q >= norm_range[0]) & (self.Q <= norm_range[1])
        norm_curves = np.array([data[:,num]/np.trapz(data[Q_index,num], 
                                x=self.Q[Q_index]) for num in range(np.shape(data)[1])])
        
        norm_curves = np.array(norm_curves).T
        self.U, self.s, self.V = np.linalg.svd(norm_curves)
        colors = [cm.jet(num/float(num_components),1) for num in range(num_components)]
        
        if produce_plot:
            fig = plt.figure()
            plt.loglog(self.s, 'o')
        
            fig = plt.figure()
            axes = plt.subplot(111)
            for num in range(num_components):
                plt.plot(self.Q, np.tile(-num/3, len(self.Q)), '-k', lw=2)
                plt.plot(self.Q, self.U[:,num]-num/3, color = colors[num], lw=2,
                         label = 'Comp %i' % int(num+1) )
            plotbox = axes.get_position()
            axes.set_position([plotbox.x0, plotbox.y0, plotbox.width * 0.8, plotbox.height])
            axes.legend(loc='upper left', bbox_to_anchor=(1, 1))
                
            
            num_frames = np.linspace(1,len(self.V),len(self.V))
            fig = plt.figure()
            axes = plt.subplot(111)
            for num in range(num_components):
                plt.plot(num_frames, np.tile(-num/3, len(self.V)), '-k', lw=2)
                plt.plot(num_frames, self.V[num,:]-num/3, color = colors[num], lw=2,
                         label = 'Comp %i' % int(num+1) )
            plotbox = axes.get_position()
            axes.set_position([plotbox.x0, plotbox.y0, plotbox.width * 0.8, plotbox.height])
            axes.legend(loc='upper left', bbox_to_anchor=(1, 1))
            
            
            
            
    def generate_subset(self, bandwidth=1, use_region=True, region=[0,1],
                            use_max_occ=False, intensity_region=[7,7.15]):
             
       frame_num = np.linspace(1,np.shape(self.IQ)[1],np.shape(self.IQ)[1])                   
       if use_region:
           self.subset = self.IQ[:, region[0]:region[1]]
           self.frames = frame_num[region[0]:region[1]]
                
       if use_max_occ:
           Q_index = (self.Q > intensity_region[0]) &  (self.Q < intensity_region[1]) 
           IQ_region = np.sum(self.IQ[Q_index,:], axis=0)
           hist, bin_edges = np.histogram(IQ_region, bins=100)
           maximum_index = np.argmax(hist)
           intensity = bin_edges[maximum_index] # not quite accurate
           min_intensity = intensity*(1-bandwidth/100)
           max_intensity = intensity*(1+bandwidth/100)
           curve_idx = (IQ_region >= min_intensity) & (IQ_region <= max_intensity)
           self.subset = self.IQ[:,curve_idx]
           self.frames = frame_num[curve_idx]
                
       fig = plt.figure()
       plt.plot(frame_num, np.sum(self.IQ, axis=0), 'ok')
       plt.plot(self.frames, np.sum(self.subset, axis=0), 'or')
            
       fig = plt.figure()
       plt.plot(self.Q, self.IQ, '-k')
       plt.plot(self.Q, self.subset, '-r')
       
       
       
    def cov_analysis(self, use_subset=True, scale_color=False, colorscale=0.8):
        if use_subset:
            data = self.subset
        else:
            data = self.IQ
        
        scaling =  np.linalg.lstsq(np.atleast_2d(data[:,-1]).T, data)
        scaled_IQ = data/scaling[0]

        self.differentials = np.array([scaled_IQ[:,num]-scaled_IQ[:,num+1] for num in range(0,np.shape(data)[1],2)])
        
        cov = np.cov(self.differentials.T)
        corr = np.corrcoef(self.differentials.T)
        
        fig = plt.figure()
        trans_data = np.log10(np.abs(cov))
        plt.imshow(trans_data)
        plt.title('log10 of absolute covariance')
        if scale_color:
            plt.clim(vmin=np.min(trans_data)*colorscale, vmax=np.max(trans_data)*colorscale)
        plt.colorbar()
        
        fig = plt.figure()
        plt.imshow(corr)
        plt.title('Correlation matrix')
        if scale_color:
            plt.clim(vmin=np.min(corr)*colorscale, vmax=np.max(corr)*colorscale)
        plt.colorbar()
        
        fig = plt.figure()
        plt.imshow(self.differentials)
        if scale_color:
            plt.clim(vmin=np.min(self.differentials)*colorscale, 
                     vmax=np.max(self.differentials)*colorscale)
        plt.colorbar()

                            
        
        
        
                
                
                
                


        
        
       
                
        
        
    
    