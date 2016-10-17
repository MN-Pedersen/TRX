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
import h5py
import fabio
import timeit

#%%

class Raw1D_dataset:
    """
    Class for analysis of azimuthally averaged scattering curves intended for 
    detector response characterisation. Please see the docstring for 
    the different methods for details
    
    Arguments
    -----    
    datapath : The path the data which is to be analysed
    
    file_extension : The file extension for the. Default is '.chi'
    
    separator : The separator used for directories. Default is '\\' (windows),
                use '/' for linux/mac.
                
    header_number : Number of lines to skip before reading the data. default is
                    24.
    
    QIQ : list of indices for Q-vector and IQ-vector, respectively. Default is
          [0, 1] which indicates first column is Q-vector and second column is 
          IQ
          
    Sort : T/F indicating whether the data should be used based on the file 
    number preceeding the file extension. Default is False
    
    
    Methods
    -----    
    load_data : loads all ascii files in the data path with the specified
                  file extension.
     
    sum_curves : Sums the intensity of the different curves. Can be used to
                   indentify outliers from e.g. storage ring injection
                   
    zinger_removal : Removes zingers by either replacing them with median or
                     excluding the entire curve
                   
    generate_subset : Method for producing a subset of the loaded data which
                        can be used instead of the entire data set in subsequent
                        methods
    
    intensity_trends : Tool for detecting non-linearities between different
                       Q-regions
                       
    svd_analysis : Implemtation of numpy's SVD
    
    cov_analysis : Method for analysing the covariance and correlation matrices
                   of the data.
    
    Notes
    -----
    
    No error handling yet
    
    """
    def __init__(self, datapath=None, file_extension='.chi', separator='\\',
                 header_number=24, QIQ=[0, 1], sort=False):
        self.datapath = datapath
        self.file_extension = file_extension
        self.header_number = header_number
        self.separator = separator
        self.QIQ = QIQ
        self.sort = sort        
        
        
    def load_data(self):
        self.files = os.listdir(self.datapath)
        
        IQ = []
        file_number = []
        for file in self.files:
            if file.endswith(self.file_extension):
                curve = np.genfromtxt(self.datapath+self.separator+file, skip_header=self.header_number)
                IQ.append(curve[:,self.QIQ [1]])
                if self.sort:
                    file_number.append(int(re.findall('\d{4}',file)[-1]))
                    
        if self.sort:
            sort_index = np.argsort(file_number)
            IQ = np.array(IQ)[sort_index,:]    
                
                
        self.Q = curve[:,self.QIQ[0]]
        self.IQ = np.array(IQ).T
        
    def preprocess_curves(self, clean=False, clean_thresshold=1e6, keep_tail=True,
                   tail_intensity=5e4, normalise=True, norm_range=[8,9]):
        self.summed_intensity = np.sum(self.IQ,axis=0)
        
        if keep_tail:
            index_keep = self.summed_intensity < tail_intensity
            self.tail = self.IQ[:,index_keep]
        
        
        if clean:
            index_remove = self.summed_intensity < clean_thresshold
            self.IQ = self.IQ[:,~index_remove]
            self.IQ = self.IQ[:770]
            self.Q = self.Q[:770]
            self.summed_intensity = np.sum(self.IQ,axis=0)
        
        if normalise:
            Q_index = (self.Q >= norm_range[0]) & (self.Q <= norm_range[1]) 
            self.IQ  = np.array([self.IQ[:,num]/np.trapz(self.IQ[Q_index,num], 
                                x=self.Q[Q_index]) for num in range(np.shape(self.IQ)[1])])
            self.IQ  = np.array(self.IQ).T
        
            
    def zinger_removal(self, std_cutoff=5, rejected='median'):
        if rejected is not 'median' and rejected is not 'exclude':
            raise ValueError('Rejected only accepts "median" or "exclude" as input')
        median = np.median(self.IQ, axis=1)
        std = np.std(self.IQ, axis=1)
        excluded = np.zeros(np.shape(self.IQ)[1], dtype=bool) # stays 0 if 'median is used'
        
        for num in range(np.shape(self.IQ)[1]):
            test = self.IQ[:,num] > median + std_cutoff*std
            
            if np.sum(test) > 0:
                excluded[num] = True
                if rejected is 'median':
                    self.IQ[test,num] = median[test]
                    
        if rejected is 'exclude':
             self.IQ = self.IQ[:,~excluded]
             
        
            
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
    

    def svd_analysis(self, num_components=4, produce_plot=True, 
                     use_subset=False, Raw=True):
         
        if use_subset and Raw:
            data = self.subset
        elif use_subset and not Raw:
            data = np.array([self.subset[:,num]-self.subset[:,num+1] for num in range(0,np.shape(self.subset)[1]-2,2)]).T
        elif not use_subset and Raw:
            data = self.IQ
        elif not use_subset and not Raw:
            data =  np.array([self.IQ[:,num]-self.IQ[:,num+1] for num in range(0,np.shape(self.IQ)[1]-2,2)]).T
            
        
        self.U, self.s, self.V = np.linalg.svd(data)
        colors = [cm.jet(num/float(num_components),1) for num in range(num_components)]
        
        if produce_plot:
            num_frames = np.linspace(1,len(self.V),len(self.V))
            fig = plt.figure()
            axes = plt.subplot(111)
            plt.title('Singular values')
            plt.loglog(num_frames[num_components-1:100], self.s[num_components-1:100],'o')
            for num in range(num_components):
                plt.loglog(num_frames[num], self.s[num], 'o', color=colors[num], label = 'Comp %i' % int(num+1))
            plotbox = axes.get_position()
            plt.ylabel('Singular values')
            plt.xlabel('Component')
            axes.set_position([plotbox.x0, plotbox.y0, plotbox.width * 0.8, plotbox.height])
            axes.legend(loc='upper left', bbox_to_anchor=(1, 1))
                
        
            fig = plt.figure()
            axes = plt.subplot(111)
            plt.title('Left singular vectors')            
            for num in range(num_components):
                plt.plot(self.Q, np.tile(-num/3, len(self.Q)), '-k', lw=2)
                plt.plot(self.Q, self.U[:,num]-num/3, color = colors[num], lw=3,
                         label = 'Comp %i' % int(num+1) )
            plotbox = axes.get_position()
            plt.ylabel('Normalized variance')
            plt.xlabel('Q (1/A)')
            axes.set_position([plotbox.x0, plotbox.y0, plotbox.width * 0.8, plotbox.height])
            axes.legend(loc='upper left', bbox_to_anchor=(1, 1))
                
            
            
            fig = plt.figure()
            axes = plt.subplot(111)
            plt.title('Right singular vectors')            
            for num in range(num_components):
                plt.plot(num_frames, np.tile(-num/3, len(self.V)), '-k', lw=2)
                plt.plot(num_frames, self.V[num,:]-num/3, color = colors[num], lw=3,
                         label = 'Comp %i' % int(num+1) )
            plotbox = axes.get_position()
            plt.ylabel('Weight')
            plt.xlabel('Frame')
            axes.set_position([plotbox.x0, plotbox.y0, plotbox.width * 0.8, plotbox.height])
            axes.legend(loc='upper left', bbox_to_anchor=(1, 1))
            
            
            
            

       
       
       
    def cov_analysis(self, use_subset=True, scale_color=False, colorscale=0.8):
        if use_subset:
            data = self.subset
        else:
            data = self.IQ
        
        scaling =  np.linalg.lstsq(np.atleast_2d(data[:,-1]).T, data)
        scaled_IQ = data/scaling[0]

        self.differentials = np.array([scaled_IQ[:,num]-scaled_IQ[:,num+1] for num in range(0,np.shape(data)[1]-2,2)])
        
        self.cov = np.cov(self.differentials.T)
        self.corr = np.corrcoef(self.differentials.T)
        
        fig = plt.figure()
        trans_data = np.log10(np.abs(self.cov))
        plt.imshow(trans_data)
        plt.title('log10 of absolute covariance')
        if scale_color:
            plt.clim(vmin=np.min(trans_data)*colorscale, vmax=np.max(trans_data)*colorscale)
        plt.colorbar()
        
        fig = plt.figure()
        plt.imshow(self.corr)
        plt.title('Correlation matrix')
        if scale_color:
            plt.clim(vmin=np.min(self.corr)*colorscale, vmax=np.max(self.corr)*colorscale)
        plt.colorbar()
        
        fig = plt.figure()
        plt.imshow(self.differentials)
        if scale_color:
            plt.clim(vmin=np.min(self.differentials)*colorscale, 
                     vmax=np.max(self.differentials)*colorscale)
        plt.colorbar()
        
        
    def point_spread(self, Q_list=None, num_points=20):
        indices = []
        for Q in Q_list:
            idx_all = self.Q == np.round(self.Q, 1)
            length_idx = len(self.Q[idx_all])
            index = np.linspace(0,self.Q-1,1)[idx_all]
            index = index[int(length_idx)/2]
            indices.append(index)
            
        data_points = np.linspace(1,num_points,1)
        
        
        
        
#%%
        
        
class Raw2D_dataset:
    """
    
    
    
    """
    def __init__(self, filename_path=None, datapath=None, mask_path=None, file_extension='.edf',
                 beam_center=[920,920], sort=True):
        self.file_path = filename_path
        self.datapath = datapath
        self.mask_path = mask_path
        self.file_extension = file_extension
        self.beam_center = beam_center        
        self.sort = sort
        
        
    
    def create_hdf5(self, separator='\\', image_size=[1918,1918]):
        self.files = os.listdir(self.datapath)
        self.mask = fabio.open(self.mask_path).data
        self.beam
        
        loadtxt = []
        file_number = []
        for file in self.files:
            if file.endswith(self.file_extension):
                loadtxt.append(self.datapath+separator+file)
                if self.sort:
                    file_number.append(int(re.findall('\d{4}',file)[-1]))
                    
        if self.sort:
            sort_index = np.argsort(file_number)
            loadtxt = np.array(loadtxt)[sort_index]

        with h5py.File(self.file_path, 'a') as file:
            dset = file.create_dataset('Raw_2D',shape=(image_size[0], 
                                    image_size[1], len(loadtxt)), dtype=np.int16) #, compression='gzip', fletcher32=True)
            start_time = timeit.default_timer()
            for num, full_path in enumerate(loadtxt):
                image = fabio.open(full_path).data # through header away
                dset[:,:,num] = image[1:image_size[0]+1, 1:image_size[0]+1]
                if np.mod(num,25) == 0:
                    elapsed = (timeit.default_timer() - start_time)/60
                    print('finished with %i frames in %0.2fmin' % (num,elapsed))
            
            
    def define_subset(self, coords=[25,75,25,75]):
        with h5py.File(self.file_path) as file:
            example = file['/Raw_2D'][:,:,0]
            
            plt.figure()
            plt.subplot(121)
            plt.imshow(example)
            plt.colorbar()
            
            plt.subplot(122)
            plt.imshow(example[coords[0]:coords[1],coords[2]:coords[3]])
            plt.colorbar()
            
            if '/Raw_subset' in file:
                del file['/Raw_subset']
            
            self.subset=file['/Raw_2D'][coords[0]:coords[1],coords[2]:coords[3],:]
            file.create_dataset('/Raw_subset', data=self.subset, dtype=np.int16)
            
    def pull_request_dummy(self):
        pass
    
        
        
        
#%% IMPLEMENT!!!

###########################
###########################
###########################
###########################
###########################



#%%
#
#data = C50k.subset
#Q = C50k.Q
#scaling =  np.linalg.lstsq(np.atleast_2d(data[:,-1]).T, data)
#scaled_IQ = data/scaling[0]

#cov = np.cov(scaled_IQ)
#corrcoef = np.corrcoef(scaled_IQ)

#fig = plt.figure(figsize=(14,10))
#trans_data = np.log10(np.abs(cov))
#plt.title('log10 of absolute Q-Covariance from total scattering curves\n30K counts data (representative)\nCurves have been scaled before comparison')
#plt.imshow(trans_data, extent=[min(Q),max(Q),max(Q),min(Q)])
#plt.xlabel('Q (1/A)')
#plt.ylabel('Q (1/A)')
#plt.colorbar()

#fit = plt.figure(figsize=(14,10))
#plt.title('Q-Correlation from total scattering curves\n30K counts data (representative)\nCurves have been scaled before comparison')
#plt.imshow(corrcoef, extent=[min(Q),max(Q),max(Q),min(Q)])
#plt.xlabel('Q (1/A)')
#plt.ylabel('Q (1/A)')
#plt.colorbar()


#%%

#indices = [50,100,138,300,500]
#buffer = 20

#C30 = np.corrcoef(C60k.differentials.T)
#C30 = corrcoef
#fig = plt.figure(figsize=(12,8))
#plt.title('Raw point spread function\nData from 30K counts')
#for index in indices:
    
    #plt.plot(corrcoef[index,:], label = 'Raw: Q=0.2f' % (Q[index]))

#    maximum = np.argmax(C30[index,:])
#    points = np.linspace(-buffer,buffer,2*buffer+1)
    

#    ax=plt.subplot(211)
#    plt.title('Differential point spread functions\nData from 60K counts')
#    plt.plot(Q[maximum-buffer:maximum+buffer+1], C30[index, maximum-buffer:maximum+buffer+1], '-', label = 'Diff: Q=%0.2f' % (Q[index]))
#    plt.xlabel('Q (1/A)')
#    plt.ylabel('Correlation')
#    plotbox = ax.get_position()
#    ax.set_position([plotbox.x0, plotbox.y0+plotbox.height*0.03, plotbox.width, plotbox.height*0.97])     
    
    
#    ax = plt.subplot(212)
#    plt.plot(points,C30[index, maximum-buffer:maximum+buffer+1],'-o', label = 'Diff: Q=%0.2f' % (Q[index]))#, markersize=10)
#    plt.xlabel('Data point')
#    plt.ylabel('Correlation')
#    plotbox = ax.get_position()
#    ax.set_position([plotbox.x0, plotbox.y0, plotbox.width*0.96, plotbox.height])
#    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))                
#                
                
                


        
        
       
                
        
        
    
    