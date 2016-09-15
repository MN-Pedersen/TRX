# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 20:29:18 2016

@author: kurt
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import re
import seaborn as sns
sns.set(font_scale=1.2)

#%%

class difference_curves:
        def __init__(self, datapath=None, file_extension='.dat', separator='\\',
                 header_number=2):
                     
            self.datapath = datapath
            self.file_extension = file_extension
            self.header_number = header_number
            self.separator = separator

        
        def load_data(self):
            self.files = os.listdir(self.datapath)
            
            Q = []
            IQ = []
            errors = []
            labels = []
            for file in self.files:
                if file.startswith('diff_av'):
                    delay = file.split(sep='_')[-1]
                    delay = delay[:-len(self.file_extension)]
                    labels.append(delay)
                    data = np.genfromtxt(self.datapath+self.separator+file)
                    Q.append(data[:,0])
                    IQ.append(data[:,1])
                    errors.append(data[:,2])
            
            # *** sort the data in order of time delay ***
            time_units = ['ps','ns','us','ms']
            time_conversions = [1e-12,1e-9,1e-6,1e-3]
            time_seconds = np.zeros(len(labels))
            
            # convert into np.array for fancy pancy indexing
            Q = np.array(Q).T
            IQ = np.array(IQ).T
            errors = np.array(errors).T
            
            for idx1, delay in enumerate(labels):
                if delay == 'off':
                    time_seconds[idx1] = -1e1
                    break
                else:
                    for idx2, unit in enumerate(time_units):
                        if delay.endswith(unit):
                            delay_number_float = float(re.findall("-?\d+\.?\d*",delay)[0])
                            time_seconds[idx1] = delay_number_float*time_conversions[idx2]
                            break # time in seconds have been calculated - go to next entry
            
            # store the data in the instance                
            indices = np.argsort(time_seconds)
            self.time_seconds = time_seconds[indices]
            self.Q = Q[:,indices]
            self.IQ = IQ[:,indices]
            self.errors = errors[:,indices]
            self.labels = np.array(labels)[indices]
            
        def SVD(self, num_components=5, UQ_multiplier=0):
            # check that all Q-vectors are similar
            Qsum = np.sum(self.Q, axis=0)
            checks = [np.allclose(Qsum[num],Qsum[0]) for num in range(1,len(Qsum))]
            if np.sum(checks) < len(checks):
                raise ValueError('Differing Q-vectors\nPlease sanitise your input data')
            
            self.U, self.s, self.V = np.linalg.svd(self.IQ) # excepts data as ROWS
            
            colours = [cm.jet(num/float(num_components),1) for num in range(num_components)]            
            
            fig = plt.figure(figsize=(12,8))
            plt.title('Singular values')
            plt.xlabel('Component')
            plt.ylabel('Singular value')
            plt.loglog(self.s,'o')
            
            fig = plt.figure(figsize=(12,8))
            axes = plt.subplot(111)
            plt.title('Left singular vectors')
            plt.xlabel('Q (1/A)')
            plt.ylabel('Variance')
            for num in range(num_components):
                plt.plot(self.Q[:,0], np.tile(-num, len(self.Q[:,0])), '-k', lw=2)
                if UQ_multiplier == 0: # integer comparisons are ok
                    plt.plot(self.Q[:,0], self.U[:,num]-num, color = colours[num], lw=2,
                         label = 'Comp %i' % int(num+1) )
                    plt.ylabel('Variance') # overwritten in loop which is ok 
                else:
                    weight = self.Q[:,0]*UQ_multiplier
                    plt.plot(self.Q[:,0], self.U[:,num]*weight-num, color = colours[num], lw=2,
                         label = 'Comp %i' % int(num+1) )
                    plt.ylabel('Variance scaled by Q*%i' % UQ_multiplier) # overwritten in loop which is ok 
            plotbox = axes.get_position()
            axes.set_position([plotbox.x0, plotbox.y0, plotbox.width * 0.8, plotbox.height])
            axes.legend(loc='upper left', bbox_to_anchor=(1, 1))
            
            
            fig = plt.figure(figsize=(12,8))
            axes = plt.subplot(111)
            plt.title('Right singular vectors')
            plt.xlabel('Delay (seconds)')
            plt.ylabel('Weigth')
            for num in range(num_components):
                plt.semilogx(self.time_seconds, np.tile(-num, len(self.V)), '-k', lw=2)
                plt.semilogx(self.time_seconds, self.V[num,:]-num, color = colours[num], lw=2,
                         label = 'Comp %i' % int(num+1) )
            plotbox = axes.get_position()
            axes.set_position([plotbox.x0, plotbox.y0, plotbox.width * 0.8, plotbox.height])
            axes.legend(loc='upper left', bbox_to_anchor=(1, 1))
                
    
                
                

                    

