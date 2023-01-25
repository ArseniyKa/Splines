#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
from decimal import Decimal
import Approximation as ap
import Extrapolation as ex
from enum import Enum

class CurveType(Enum):
    SATURATION1 = 1
    PERMEABILITY1 = 2   ### which increases
    PERMEABILITY2 = 3
            

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
        
        
    def reverseSplitData(self, x, y,  spline_number, x_board):
        if (spline_number ==1):
            X = [x[np.logical_and(x_board[i-1 ] < x, x <= x_board[i])] for i in range(1, x_board.size)]
            Y = [y[np.logical_and(x_board[i-1] < x, x <= x_board[i])] for i in range(1, x_board.size)]
        else:
######################################################
            index_first = np.logical_and(x_board[0] < x, x <= (x_board[1] + x_board[2])/2 )
            index_last = np.logical_and((x_board[-2] + x_board[-3])/2 < x, x <= x_board[-1] )

            indexes = [np.logical_and((x_board[i - 1] + x_board[i])/2 < x, x <= (x_board[i+1] + x_board[i+2])/2 ) for i in range(1, x_board.size - 2)]
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
        self.curve_type = curve_type
        self.extr_type = extr_type
        
        self.last_spline_number = last_spline_number
        self.last_coeff_number = last_coeff_number  
        
        self.x_board = np.linspace(begin, end, spline_number+1) 
        
        self.X, self.Y =  self.splitData(x, y,  spline_number, self.x_board) if curve_type !=CurveType.PERMEABILITY1 else self.reverseSplitData(x, y,  spline_number, self.x_board)
        
        if curve_type == CurveType.SATURATION1:         
            self.last_spline_flag= True if   last_spline_number >1 and spline_number >1 else False
            self.extrapolation_flag = True if extr_type !=ex.ExtrapolatedType.NEITHER else False
        else:
            self.last_spline_flag =False
            self.extrapolation_flag =False
        
        if curve_type ==CurveType.PERMEABILITY1:
            self.reverse_flag=False
        else:
            self.reverse_flag = False
        
                                                                        
                         
   #############################################CALCULATE COEFFICIENTS MATRIX 
    def getCoefficients(self, last_flag):
  
        curve_type = self.curve_type
 ############# CALCULATE LAST COEFF MATRIX 
        if last_flag:               
                yo = ap.calculateCurveValue(self.x_board[-2], self.coeff_matrix[-1]);
                dyo = ap.calculateDerivative(self.x_board[-2], self.coeff_matrix[-1]);
                
                x_board = self.last_x_board
                coeff_number =  self.last_coeff_number
                size = self.last_spline_number
                
                x1 = self.last_X[-1][-1]
                y1 = self.last_Y[-1][-1]
                
                X = self.last_X
                Y = self.last_Y
                flag =True
        else: 
            xo= 1 if  curve_type !=CurveType.PERMEABILITY1 else  0
            yo = 0
            dyo=0
                                                                      
            x_board = self.x_board
            coeff_number =  self.coeff_number
            size = self.size
            
            x1 =0  if  curve_type !=CurveType.PERMEABILITY1 else  1
            y1 = 1
                                    
            X = self.X
            Y = self.Y
            flag = True if curve_type !=CurveType.SATURATION1 else False
            
        
        coeff_matrix =np.zeros(( size ,coeff_number), np.float64)
        apr = ap.Approximation(coeff_number)
        
        for i in range(size):
            if i == 0 and not last_flag:
                func_settings = ap.ValueSettings(xo, yo, 0,  0,0, 1)               
            elif i== size - 1 and flag:
                func_settings = ap.ValueSettings(x_board[i], yo, dyo,  x1, y1, 3)                                 
            else:
                func_settings = ap.ValueSettings(x_board[i], yo, dyo, 0,0, 2)
                
                
            coeff_matrix[i] = apr.conditionalFitting(X[i], Y[i], func_settings)
            yo = ap.calculateCurveValue(x_board[i + 1],  coeff_matrix[i]);
            dyo = ap.calculateDerivative(x_board[i + 1], coeff_matrix[i]);
        
        return coeff_matrix
            
            
                
 
    
    def calculateCoefficients(self):
        self.size = self.spline_number -1 if self.last_spline_flag else self.spline_number
        self.coeff_matrix = self.getCoefficients(False)
        coeff_list = [elem for elem in self.coeff_matrix]    
        
        if  self.last_spline_flag:
            self.last_x_board = np.linspace(self.x_board[-2], self.end, self.last_spline_number + 1)  
            self.last_X, self.last_Y = self.splitData(self.X[-1], self.Y[-1],  self.last_spline_number, self.last_x_board)
            
            self.last_coeff_matrix = self.getCoefficients(True)
            coeff_list.extend(np.array(self.last_coeff_matrix))
           
            
      ################################################  EXTRAPOLATION
        if self.extrapolation_flag:
                self.extr = ex.Extrapolation(self.end, self.discrete_size, self.extr_type)
                coeff_vec = self.last_coeff_matrix[-1] if self.last_spline_flag else self.coeff_matrix[-1]
                self.extr.extrapolation(coeff_vec)  
                coeff_list.extend( np.array(self.extr.coeff_vec) )
        

        return coeff_list
             
    def makeSplines(self, coeff_matrix, x_board, size, left, right):
        discrete_size = self.discrete_size
        
        temp = np.linspace(left, right, size * discrete_size)
        
        if not self.reverse_flag:
            indexes = [np.logical_and( x_board[i] >= temp, temp >= x_board[i+1] ) for i in range(size)]
        else:
            indexes = [np.logical_and( x_board[i] <= temp, temp <= x_board[i+1] ) for i in range(size)]
        
        X_temp = [temp[ind] for ind in indexes]
        splines = np.array([ap.calculateCurveValue(X_temp[i], coeff_matrix[i]) for i in range(size)])         
        splines = np.array([elem if elem>0  else 0  for spline in splines for elem in spline])
        
        return [temp, splines.reshape(-1)]
        
        
    def createSplines(self):
      
        x_board = self.x_board
        size = self.size
        discrete_size = self.discrete_size
        
        left = self.begin
        right = self.end if self.spline_number==1 or not self.last_spline_flag else self.x_board[-2]
        temp, splines = self.makeSplines(self.coeff_matrix, x_board, size, left, right)           
        all_splines = [[temp, splines]] 

   ######################## LAST SPLINES         
        if self.last_spline_flag:
            last_x_board = self.last_x_board
            size = self.last_spline_number
            left = last_x_board[0]
            right = last_x_board[-1]
            
            last_temp, last_splines = self.makeSplines(self.last_coeff_matrix, last_x_board, size, left, right)
            all_splines.extend([[ last_temp, last_splines ]])
            
  ######################## EXTRAPOLATED SPLINE  
        if  self.extrapolation_flag:
              all_splines.extend([[self.extr.temp, self.extr.spline]])
   
          
        return all_splines
        
        
       
       





