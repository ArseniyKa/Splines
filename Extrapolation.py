#!/usr/bin/env python
# coding: utf-8

import numpy as np
from enum import Enum
import Approximation as ap


class ExtrapolatedType(Enum):
    EXPONENT = 1
    PARABOLA = 2
    LINEAR = 3
    NEITHER = 4
    
    
class Extrapolation:
    def __init__(self, xo, discrete_size, extr_type):
        self.xo= xo
        self.discrete_size= discrete_size
        self.temp =  np.linspace(0, xo, discrete_size)
        self.spline = np.zeros(discrete_size)
        self.extr_type = extr_type
    
    def setData(self, coeff_vector):
        self.yo = ap.calculateCurveValue(self.xo, coeff_vector)
        self.dyo = ap.calculateDerivative(self.xo, coeff_vector)
        
     
    def exponent(self):   
        B = self.dyo / self.yo
        A = self.yo * np.exp(-1*B * self.xo)
#         print(self.extr_temp)
        self.spline = A * np.exp(B * self.extr_temp)
        self.coeff_vec = ([[A, B]])
      
    
    def parabola(self):
        a=  self.dyo/(2*self.xo) 
        c =  self.yo - self.dyo * self.xo / 2 
        self.spline = a * self.extr_temp**2 + c
        self.coeff_vec = ([[c, 0, a]])
       
    
    def linear(self):
        a_lin = self.yo - self.dyo*self.xo
        b_lin = self.dyo
        self.spline = a_lin + b_lin * self.temp
        self.coeff_vec = ([[a_lin, b_lin]])
        
        
    def extrapolation(self, input_vec):        
            self.setData(input_vec)
               
            if self.extr_type == ExtrapolatedType.LINEAR:
                self.linear()
            elif self.extr_type == ExtrapolatedType.PARABOLA:
                self.parabola()
            elif self.extr_type == ExtrapolatedType.EXPONENT:
                self.exponent()  
       
    

        


