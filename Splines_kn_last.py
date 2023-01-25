#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
from decimal import Decimal
import Approximation as ap
import Extrapolation as ex
from enum import Enum

class ExtrapolatedType(Enum):
    EXPONENT = 1
    PARABOLA = 2
    LINEAR = 3
    NEITHER = 4


class CurveType(Enum):
    SATURATION1 = 1
    SATURATION2 = 2
    PERMEABILITY1 = 3   ### which increases
    PERMEABILITY2 = 4
            

class Splines:
    
    def splitData(self, x, y,  spline_number, x_board):
        if (spline_number ==1):
            X = [x[np.logical_and(x_board[i-1 ] > x, x >= x_board[i])] for i in range(1, x_board.size)]
            Y = [y[np.logical_and(x_board[i-1] > x, x >= x_board[i])] for i in range(1, x_board.size)]
        else:
######################################################
            index_first = np.logical_and(x_board[0] > x, x >= (x_board[1] + x_board[2])/2 )
            index_last = np.logical_and((x_board[-2] + x_board[-3])/2 > x, x >= x_board[-1] )

            indexes = [np.logical_and((x_board[i - 1] + x_board[i])/2 > x, x >= (x_board[i+1] + x_board[i+2])/2 ) for i in range(1, x_board.size - 2)]
            indexes.insert(0, index_first)
            indexes.append(index_last)

            X = [x[ind] for ind in indexes]
            Y = [y[ind] for ind in indexes]
            
            return [X, Y]
            
            
    def __init__(self, x, y, spline_number, coeff_number, discrete_size, begin, end, last_spline_number, last_coeff_number, extr_type, curve_type):
        self.spline_number = spline_number
        self.coeff_number = coeff_number
        self.discrete_size = discrete_size      
        self.begin = begin        
        self.end = end   
        last_spline_number = last_spline_number if curve_type != CurveType.PERMEABILITY1 else 1 
        self.last_spline_number = last_spline_number   
        self.last_coeff_number = last_coeff_number
        
        self.x_board = np.linspace(begin, end, spline_number+1)        
        self.coeff_size = spline_number if (spline_number ==1 or last_spline_number == 1) else spline_number - 1
        self.coeff_matrix =np.zeros((self.coeff_size ,coeff_number), np.float64)
        self.curve_type = curve_type
        
        self.X, self.Y =  self.splitData(x, y,  spline_number, self.x_board)
        if last_spline_number > 1:
            self.last_x_board = np.linspace(self.x_board[-2], end, last_spline_number + 1)
            self.last_coeff_matrix =np.zeros((last_spline_number,last_coeff_number))
            
            self.last_X, self.last_Y = self.splitData(x, y,  last_spline_number, self.last_x_board)
            
        
        self.extr_type = extr_type if curve_type == CurveType.SATURATION1 else ExtrapolatedType.NEITHER
        if extr_type != ExtrapolatedType.NEITHER:
            self.extr = ex.Extrapolation(end, discrete_size)
            
                 
               
        
  #############################################CALCULATE COEFFICIENTS MATRIX          
    def calculateCoefficients(self):
        yo =0
        dyo=0
        curve_type = self.curve_type
        right_value = 1 if curve_type ==CurveType.PERMEABILITY1  else 0
        permeability_flag = True if curve_type ==CurveType.PERMEABILITY1   else False
        left_value = 0
        
        apr = ap.Approximation(self.coeff_number)
        
        for i in range(self.coeff_size):
            if i == 0:
                func_settings = ap.ValueSettings(1, right_value, 0,  0,0, 1)
                
            elif i== self.coeff_size - 1 and permeability_flag == True:
                func_settings = ap.ValueSettings(self.x_board[i], yo, dyo,  0, left_value, 3)
                                    
            else:
                func_settings = ap.ValueSettings(self.x_board[i], yo, dyo, 0,0, 2)
                
                
            self.coeff_matrix[i, :] = apr.conditionalFitting(self.X[i], self.Y[i], func_settings)
            yo = ap.calculateCurveValue(self.x_board[i + 1], self.coeff_matrix[i, :]);
            dyo = ap.calculateDerivative(self.x_board[i + 1], self.coeff_matrix[i, :]);
            
        return self.coeff_matrix
                
    
 ######## CALCULATE LAST COEFF MATRIX    
    def calculateLastCoefficients(self):
        curve_type = self.curve_type
        last_x_board = self.last_x_board
        yo = ap.calculateCurveValue(self.x_board[-2], self.coeff_matrix[-1, :]);
        dyo = ap.calculateDerivative(self.x_board[-2], self.coeff_matrix[-1, :]);
        
        apr= ap.Approximation(self.last_coeff_number)
        
        
        for i in range(self.last_spline_number):
            if i == self.last_spline_number - 1:
                x1  = np.min(self.last_X[i]) if curve_type == CurveType.SATURATION1 else 0
                y1 = np.max(self.last_Y[i]) if curve_type == CurveType.SATURATION1 else 1
                func_settings = ap.ValueSettings(last_x_board[i], yo, dyo,  x1, y1, 3)
            else:
                func_settings = ap.ValueSettings(last_x_board[i], yo, dyo,  0,0, 2)
                
            self.last_coeff_matrix[i, :] = apr.conditionalFitting(self.last_X[i], self.last_Y[i], func_settings)
            yo = ap.calculateCurveValue(last_x_board[i + 1], self.last_coeff_matrix[i, :]);
            dyo = ap.calculateDerivative(last_x_board[i + 1], self.last_coeff_matrix[i, :]);
            
        return self.last_coeff_matrix
    
    def calculateAllCoefficients(self):
        self.calculateCoefficients()
        coeff_list=[elem for elem in self.coeff_matrix]
            
        
        if  self.last_spline_number > 1:
            self.calculateLastCoefficients()
            size = self.last_coeff_matrix.size
            coeff_list.extend( self.last_coeff_matrix)
  
        if  self.extr_type != ExtrapolatedType.NEITHER:
              self.extrapolation()  
              coeff_list.extend( np.array(self.extr.coeff_vec) )
        

        return coeff_list
        
    ################################################  EXTRAPOLATION
    def extrapolation(self):
            if self.last_spline_number ==1:
                self.extr.setData(self.coeff_matrix[-1, :])
            else:
                self.extr.setData(self.last_coeff_matrix[-1, :])
            
            if self.extr_type == ExtrapolatedType.LINEAR:
                self.extr.linear()
            elif self.extr_type == ExtrapolatedType.PARABOLA:
                self.extr.parabola()
            elif self.extr_type == ExtrapolatedType.EXPONENT:
                self.extr.exponent()    
    
    def createSplines(self):
        curve_type = self.curve_type
        x_board = self.x_board
        
        right_side = self.end if self.spline_number==1 or self.last_spline_number==1 else self.x_board[-2]
        
        splines = np.zeros((self.coeff_size, self.discrete_size ))
        temp = np.linspace(self.begin, right_side, splines.size)
        
        for i in range(self.coeff_size):
            X_temp = temp[np.logical_and(temp <= x_board[i], temp >= x_board[i+1])]
            splines[i] = ap.calculateCurveValue(X_temp, self.coeff_matrix[i, :])
           
        splines = splines.reshape(splines.size)
        
        for i in range(splines.size):
            if (splines[i]<0):
                splines[i] = 0
                
        all_splines = [] 
        all_splines.extend([[temp, splines]])
   ######################## LAST SPLINES         
        if (self.last_spline_number>1):
            last_x_board = self.last_x_board
            last_splines = np.zeros((self.last_spline_number, self.discrete_size), np.float64)
            last_temp = np.linspace(last_x_board[0], last_x_board[-1], last_splines.size) 
                      
            for i in range(self.last_spline_number):
                X_temp = last_temp[np.logical_and(last_temp <= last_x_board[i],last_temp >= last_x_board[i+1])]
                last_splines[i] = ap.calculateCurveValue(X_temp, self.last_coeff_matrix[i, :])
                
            last_splines = last_splines.reshape(last_splines.size)
            all_splines.extend([[last_temp, last_splines]])
            
  ######################## EXTRAPOLATED SPLINE  
        if self.extr_type != ExtrapolatedType.NEITHER:
              all_splines.extend([[self.extr.extr_temp, self.extr.extr_spline]])
   
          
        return all_splines
        
        
       
       





