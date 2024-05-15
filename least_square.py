# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 19:47:40 2022

@author: LuizF
"""

import numpy as np

class least_square(object):
    
    """
    Function which adjust a set of data with least squares method. 
    Like a pure cosine in data: amplitude*cos(omega*x + phase)
        
    Parameters
    ----------
        x : array
            Values of variable independent 
            
        y : array 
            Values of variable dependendent
            
        period: Number (integer or float)
            Period for compute angular frequency od oscillation. Default: 1
    Atributes:
    ---------
        get_values: Array
        Return values adjust the input data
        
        mean: Number (float)
        Return mean, phase, amplitude
        
    Examples:
    ----------
        y =  np.array([31,35,37,33,28,20,16,15,18,23,31])
        x =  np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])

        fitting = least_square(x, y)
        fitting.results()
        
        >>> (25.638625278740392, -1.1034612682536515, 11.043342567715348)
        fitting.get_values()
        
        >>> array([30.61374721, 35.45866844, 36.55266696, 33.4778725 , 27.40875206,
       20.66350334, 15.81858211, 14.7245836 , 17.79937805, 23.8684985 ,
       30.61374721])    
    """
    
    
    def __init__(self, x, y, period = 1):
        
        #declare variable dependent and independent
        self.x = x
        self.y = y
        
        #compute de frequency angular for period 
        omega = 2*np.pi / period
        
        #Configure the array max(rows_x) X (3): vertical stack and transpose
        matrix_sol = np.vstack([np.ones(len(self.x)), np.cos(omega*self.x), 
                                np.sin(omega*self.x)]).T
        
        #solve the linear system Ax = B
        self.mean_, a, b = np.linalg.lstsq(matrix_sol, self.y, rcond=None)[0]
    
        #Compute de amplitude 
        self.amplitude_ = np.sqrt(a**2 + b**2)
        
        # Conditions for phase 
        if a > 0:
            phase = np.arctan(- b / a) 
        elif a < 0 and b > 0:
            phase = - (np.arctan(- b / a) - np.pi)
            
        elif a < 0 and b < 0:
            phase =  (np.arctan(- b / a) + np.pi)
            
        elif a == 0 and b > 0:
            phase = np.pi / 2
            
        elif a == 0 and b < 0:
            phase = - np.pi / 2
            
        else:
            "When the amplitude is very smaller, the phase not is well defined"
            #Debug the value error for a == 0 and b == 0
            raise ValueError("Phase it is arbitrary")
            
        #change variable
        self.phase_ = phase
        self.omega = omega
        self.phase2_ = np.arctan2(b, a)
       
    @property
    def get_values(self):
        "The values for adjust"
        return self.mean_ + self.amplitude*np.cos(self.omega*self.x + 
                                                  self.phase_)
    @property
    def mean(self):
        return self.mean_  
    @property
    def phase(self):
        return self.phase_
    @property
    def amplitude(self):
        return self.amplitude_
    @property
    def phase2(self):
        "This return values negative (inverted)"
        return self.phase2_
    
