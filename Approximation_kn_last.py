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
#         self.constr_number = constr_number
        
    def simpleFitting(self, x, y):
            polinom_size = self.polinom_size
            vandermond_tr = np.zeros((polinom_size,x.size))

            for i in range(polinom_size):
                  vandermond_tr[i] = x**i

            vandermond = vandermond_tr.transpose()

            XX = vandermond_tr.dot(vandermond)
            XX_inv = np.linalg.inv(XX)
            XY = vandermond_tr.dot(y)
            coeff_vector = XX_inv.dot(XY)

            return coeff_vector
    
    
    def conditionalFitting(self, x, y, value_settings):
        polinom_size = self.polinom_size
        
        xo=value_settings.xo
        yo = value_settings.yo
        dyo = value_settings.dyo
        
        x1 = value_settings.x1
        y1  = value_settings.y1
        
        constr_number = value_settings.constr_number

        vandermond_tr = np.zeros((polinom_size,x.size))
        lagrange_multypliers_matrix = np.zeros((constr_number, polinom_size))
        init_conditial_matrix = np.zeros((constr_number, polinom_size))
        function_vec = np.zeros(polinom_size)
        derivative_vec = np.zeros(polinom_size)
        
        function_vec2 = np.zeros(polinom_size)

        for i in range(polinom_size):
            function_vec[i] = xo**i
            derivative_vec[i] = i * xo**(i-1) if constr_number >1 and xo!=0 else 0
            vandermond_tr[i] = x**i
            
            function_vec2[i] = x1**i if constr_number==3 else 0

    ############################################### Fill initial conditional and Lagrange matrixes
        if constr_number<3:
            for j in range(constr_number):
                init_conditial_matrix[j] = function_vec if j==0 else derivative_vec
                lagrange_multypliers_matrix[j] = -0.5 * ( function_vec if j==0 else derivative_vec )
        else:
            for j in range(constr_number):
                if j==0:
                    init_conditial_matrix[j] = function_vec 
                    lagrange_multypliers_matrix[j] = -0.5 * function_vec 
                elif j==1:
                    init_conditial_matrix[j] = derivative_vec 
                    lagrange_multypliers_matrix[j] = -0.5 * derivative_vec 
                else:
                    init_conditial_matrix[j] = function_vec2 
                    lagrange_multypliers_matrix[j] = -0.5 * function_vec2 
                    


     ####################### Simple matrix fitting   

        vandermond = vandermond_tr.transpose()        
        XX = vandermond_tr.dot(vandermond)
        XY = vandermond_tr.dot(y)
        
        XY = np.append(XY, yo)
        XY = np.append(XY, dyo) if constr_number >= 2 else XY
        XY = np.append(XY, y1) if constr_number == 3 else XY


    ##################### Create matrix for Algebraic system solving

        lagrange_multypliers_matrix = lagrange_multypliers_matrix.transpose()
        init_conditial_matrix = np.hstack(( init_conditial_matrix, np.zeros((constr_number,constr_number)) ))

        A = np.hstack((XX,lagrange_multypliers_matrix))
        A = np.vstack((A, init_conditial_matrix))

    #     print(lagrange_multypliers_matrix)
    #     print("------")
    #     print(A)


        coeff_vector = LA.solve(A, XY)
        ### to remove Lagrange variables
        coeff_vector = coeff_vector[:polinom_size]
    #     print ("coeff_vector", coeff_vector)
        return coeff_vector


    ###################################### FUNCTION CURVE
def calculateCurveValue(xo, coeff_vector):
        # print("coeff_vector", coeff_vector)
        value = 0
        for n in range(coeff_vector.size):
            value += coeff_vector[n] * xo**n

        return value
    ########################################### DERIVATIVE CURVE

def calculateDerivative(xo, coeff_vector):
        derivative = 0

        for n in range(1, coeff_vector.size):
            derivative += n* coeff_vector[n] * xo**(n-1)

        return derivative

    ########################################### SECOND DERIVATIVE CURVE

def calculateSecondDerivative(x, coeff_vector):
        second_derivative = 0

        for n in range(2, coeff_number):
            second_derivative += n* (n - 1) * coeff_vector[n] * xo**(n-2)

        return second_derivative
    ###############################################################

