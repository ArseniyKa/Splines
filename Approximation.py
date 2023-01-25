#!/usr/bin/env python
# coding: utf-8

# In[1]:


# !pip install pandas
# !pip install odfpy
import pandas as pd
import numpy as np
import matplotlib
from decimal import Decimal
from numpy import linalg as LA

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

class ValueSettings:
      def __init__(self, xo, yo, dyo, x1, y1, constr_number):
        self.xo = xo
        self.yo = yo
        self.dyo = dyo    
        self.x1 = x1
        self.y1 = y1   
        self.constr_number = constr_number
        
    

class Approximation:
    def __init__(self, polinom_size):
        self.polinom_size = polinom_size
        
    def simpleFitting(self, x, y):
            polinom_size = self.polinom_size
            
            vandermond = np.array([x **i for i in range(polinom_size)]).T
            XX = vandermond.T.dot(vandermond)
            XY = vandermond.T.dot(y)
            coeff_vector = LA.solve(XX, XY)

            return coeff_vector
    
    
    def conditionalFitting(self, x, y, value_settings):
        polinom_size = self.polinom_size
        constr_number = value_settings.constr_number
        
        xo=value_settings.xo
        yo = value_settings.yo
        dyo = value_settings.dyo
        
        x1 = value_settings.x1
        y1  = value_settings.y1
       
        xo_polinom = np.array([ xo**i for i in range(polinom_size) ])
        dxo_polinom = np.array([ i * xo**(i-1) if xo!=0  else 0 for i in range(polinom_size) ])
        x1_polinom = np.array([ x1**i for i in range(polinom_size) ])
           
    ############################################### Fill initial conditions and Lagrange matrixes
        conditions = xo_polinom
        conditions = np.append(conditions, dxo_polinom) if constr_number >=2 else conditions
        conditions = np.append(conditions, x1_polinom) if constr_number ==3 else conditions
        conditions = conditions.reshape(constr_number, polinom_size)
        
        lagrange_multypliers =  -0.5 * conditions.T
        conditions = np.hstack(( conditions, np.zeros((constr_number,constr_number)) ))

     ####################### Simple matrix fitting   
        vandermond = np.array([x **i for i in range(polinom_size)]).T
        XX = vandermond.T.dot(vandermond)
        XY = vandermond.T.dot(y)
        
        XY = np.append(XY, yo)
        XY = np.append(XY, dyo) if constr_number >= 2 else XY
        XY = np.append(XY, y1) if constr_number == 3 else XY

    ##################### Create matrix for Algebraic system solving
        
        A = np.hstack((XX,lagrange_multypliers))
        A = np.vstack((A, conditions))
        coeff_vector = LA.solve(A, XY)
        
 ### to remove Lagrange variables
        return coeff_vector[:polinom_size]


    ###################################### FUNCTION CURVE
def calculateCurveValue(xo, coeff_vector):
        value = sum([ coeff_vector[n] * xo**n for n in range(coeff_vector.size) ])
        return value
    ########################################### DERIVATIVE CURVE

def calculateDerivative(xo, coeff_vector):
        derivative = sum([ n * coeff_vector[n] * xo**(n-1) for n in range(1, coeff_vector.size) ])
        return derivative

    ########################################### SECOND DERIVATIVE CURVE

def calculateSecondDerivative(x, coeff_vector):
        second_derivative = sum([ n* (n - 1) * coeff_vector[n] * xo**(n-2) for n in range(2, coeff_vector.size) ])
        return second_derivative
    ###############################################################

