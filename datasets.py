# -*- coding: utf-8 -*-
"""
Class objects for data sets

MNP
"""

#%% import libraries

import numpy as np





#%%

class total_scattering_dataset(Q=0, IQ=0, sigmaQ=0, detector_dist=35, energy=18, reference=None, list_delays=None):
    """
    Class for storing time resolved X-ray scattering data sets
    
    Parameters
    ----------
    Q : The momentum transfer vector (X-axis)
    
    IQ : The recorded scattering intensities (scale not considered). Must be 
         2D array where each column is a curve
    
    sigmaQ : The standard error of the measurents
    
    detector_dist : Sample to detector distances in units of mm
    
    reference : The time delay used as reference in the subtraction
    
    list_delays : the list of the time delays in the data set exactly as they
                  were recorded. The list should hence include repitions if one 
                  time point was measured several times in the series. Example:
                  '-3ns, 100ps, 200ps, -3ns, 100ps, 300ps'
    
    
    
    Attributes
    -------
    Same as above
    
    Notes
    -----
    None yet
    
    """
    
    def __init__(self, Q, IQ, sigmaQ, detector_dist, reference, list_delays):
        # the easy variables
        self.Q = Q
        self.detector_dist = detector_dist
        self.reference = reference
        
        # cleaning the delay list
        list_delays = [delay.replace(' ','') for delay in list_delays.split()]
        self.list_delays = list_delays
        self.unique_delays = np.unique(list_delays)
        
        # reshape the two matrices to prober shape
        self.IQ_raw = IQ #np.reshape(IQ,(len(Q),-1,len(list_delays)))
        self.sigmaQ_raw = sigmaQ #Snp.reshape(sigmaQ,(len(Q),-1,len(list_delays)))
    
    def sort_delay_(self):
        
        
        pass
        
    def calculate_mean_(self):
        """
        

        
        Attributes
        -----
        adds mean_dIQ, mean_dsQ which is the average differential scattering 
            intensity and differential standard error, respectivly. also 
            
        """        
        length_Q = len(self.Q)
        number_delays = len(self.unique_delays)
        mean_dIQ = np.zeros((length_Q, number_delays))
        mean_dsQ = np.zeros((length_Q, number_delays))
        
        for num, delay in range(self.unique_delays):
            pass
            
    
        
        

    def median_outlier_detection_(self, scale=2, reject_fraction=10):
        """
        
        Parameters
        -----
        Scale : Scale used to determine if data points (I(Q)) are outliers
        
        reject_fraction : The minimum percent of data points required to label
                        a data curve 'outlier'
                        
        Attributes
        -----
        adds 'inliers_indices' which is a 2D matrix (len(Q), len(list_delays))
        
        """
        
        pass
    
    def subtract_reference_(self):
        """
        
        The reference is calulated according to the following procedure:
        if image number i and j are two references i taken before the 
        non-reference image k, and j taken after, the reference to use for 
        image k is I_i + (I_j-I_i)/(j-i)*(k-i). This minimize the effect of 
        slow drifts.        
                
        
        Returns
        -----
        dIQ : The differences curves obtained by subtracting the weigthed
            average of reference from the total scattering curves
            
        Notes
        -----
        We use return instead of storing the result in the instance to limit
        memory consumption
        
        """
        reference_index = np.array(self.list_delays) == self.reference
        numbers_reference = np.linspace(0, len(self.list_delays), 
                                        dtype=int)[reference_index]
        number_blocks = np.floor(np.shape(self.IQ)[1]/len(self.list_delays))
        
        for num, delay in enumerate(self.list_delays):
            
            # produce array of the number of the frames
            file_numbers = np.linspace(0, number_blocks, dtype=int)+num
            
            # figure out what data curves to use
            try:
                number_i = numbers_reference[numbers_reference < num][-1]
            except IndexError: # there is no reference before in the list_delay
                number_i = numbers_reference[numbers_reference > num][0]
            
            numbers_reference_before = np.linspace(0, number_blocks, 
                                                       dtype=int)+number_i
                
            try:
                number_k = numbers_reference[numbers_reference > num][0]
            except IndexError: # there is no reference after in the list
                               # use the frame from the next round instead
                first_reference =  numbers_reference[0]
                number_k = number_i+(len(self.list_delays)-number_i)+\
                                                       first_reference
                                    
             
            numbers_reference_after = np.linspace(0, number_blocks, 
                                                   dtype=int)+number_k 
            # fetch the data
            #reference_before = self.IQ[number_i]
            
            
            
            
        
       
    def scale_data_(self, Qmin=0, Qmax=1):
        """
       
        Parameters
        -----
        Qmin : The lower bound of the range used for scaling
        
        Qmax : the upper bound of the range used for scaling
        
        Attributes
        -------
        Same as above, but scales IQ
       
        adds the 'Scaled_flag'
       
        """
       
        pass
   