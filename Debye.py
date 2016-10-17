# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 18:02:17 2016

@author: denis
"""

from math import pi
import numpy as np
import matplotlib.pyplot as plt


class Debye:
    
    def __init__(self):
        pass
        #self.ZXYZ = ZXYZ
        # dummy comment2
        #self.q = q        
        
    def DebyeScat_fromZXYZ(self, ZXYZ, q):
        
        Elements = self.getElements(ZXYZ)
        atomForm = self.getAtomicFormFactor(Elements,q)
                
        S = np.zeros(q.shape)
        for i,item in enumerate(ZXYZ):
            xyz_i = np.array(item[1:])
            f_i = atomForm[item[0]]
            
            S += f_i**2
            
            for j,jtem in enumerate(ZXYZ[:i]):
                xyz_j = np.array(jtem[1:])
                r_ij = np.sqrt(np.sum((xyz_i - xyz_j)**2))                
                f_j = atomForm[jtem[0]]
                
                S[q!=0] += 2*f_i[q!=0]*f_j[q!=0]*np.sin(q[q!=0]*r_ij)/(q[q!=0]*r_ij)
                S[q==0] += 2*f_i[q==0]*f_j[q==0]
        
        return S
        
    
    def ZXYZtoGR(self, ZXYZ, Rmin = 0, Rmax = 1e2, dR = 1e-2):
        
        Elements = self.getElements(ZXYZ)
        Rpts = (Rmax-Rmin)/(dR)
        
        r = np.linspace(Rmin,Rmax,Rpts+1)
        r_bins = np.linspace(Rmin-dR/2, Rmax+dR/2, Rpts+2)
        
        gr = list()
        for i,item in enumerate(Elements):
            for j,jtem in enumerate(Elements[:i+1]):
                xyz_loc = np.array(list(x[1:] for x in ZXYZ if x[0]==item or x[0]==jtem))
                dist = np.sqrt(np.subtract(xyz_loc[:,[0]],xyz_loc[:,[0]].T)**2 + \
                               np.subtract(xyz_loc[:,[1]],xyz_loc[:,[1]].T)**2 + \
                               np.subtract(xyz_loc[:,[2]],xyz_loc[:,[2]].T)**2).ravel()
                
                gr_ij = np.histogram(dist,r_bins)[0]
                gr.append([[item,jtem],gr_ij])
        return r, gr
                        
    


    def DebyeScat_fromGR(self, r, gr, q):
        Elements = list(set(x[0][0] for x in gr))
        atomForm = self.getAtomicFormFactor(Elements,q)        
        
        QR = q[np.newaxis].T*r[np.newaxis]
        Asin = np.sin(QR)/QR
        Asin[QR==0] = 1;
        
        S = np.zeros(q.shape)
        for item in gr:
            f_i = atomForm[item[0][0]][np.newaxis]
            f_j = atomForm[item[0][1]][np.newaxis]
            S += np.squeeze(f_i.T*f_j.T*np.dot(Asin,item[1][np.newaxis].T))
        
        return S
            
        
    
    
    def getAtomicFormFactor(self,Elements,q):
        # Returns a dictionary of atomic formfactors using the element name as a key        
        
        ###        
        # define the scattering vector used for calculations of atomic formfactor:
        s=q/(4*pi)
        
        # The most modern and advanced parameterization of f0 form factors is
        # given by WaasKirf file.
        fname = 'f0_WaasKirf.dat'
        
        # This is a functional form for this parameterization
        formFunc = lambda s,a: np.sum(np.reshape(a[:5],[5,1])*np.exp(-a[6:,np.newaxis]*s**2),axis=0)+a[5]
            
        with open(fname) as f:
            content = f.readlines()    
        # Read the coefficients for the calculation:
        atomData = list()
        for i,x in enumerate(content):
            if x[0:2]=='#S':
                atomName = x.rstrip().split()[-1]
                if any([atomName==x for x in Elements]):
                    atomCoef = content[i+3].rstrip()
                    atomCoef = np.fromstring(atomCoef, sep=' ')
                    atomData.append([atomName, atomCoef])
    
    
        atomData.sort(key=lambda x: Elements.index(x[0]))
        atomForm = {}
        for x in atomData:
            atomForm[x[0]] = formFunc(s,x[1])
    
        return atomForm




    def getElements(self,ZXYZ):    
        Elements = list(set(x[0] for x in ZXYZ))
        return Elements