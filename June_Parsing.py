#!/usr/bin/env python
# coding: utf-8

import Splines as sp
import Extrapolation as ex
import numpy as np
import matplotlib.pyplot as plt


class Settings:
       def __init__(self, spline_number, coeff_number, discrete_size, last_spline_number, last_coeff_number):
            self.spline_number = spline_number
            self.coeff_number = coeff_number
            
            self.discrete_size = discrete_size   
  
            self.last_spline_number = last_spline_number   
            self.last_coeff_number = last_coeff_number
        
def splineStep(x, y, curve_type, full_name, name, settings ):
    spline_number = settings.spline_number
    coeff_number = settings.coeff_number

    last_spline_number = settings.last_spline_number
    last_coeff_number = settings.last_coeff_number

    discrete_size= settings.discrete_size
    if curve_type != sp.CurveType.PERMEABILITY1: 
        begin = 1
        end = np.min(x) * 0.99 if curve_type == sp.CurveType.SATURATION1 else 0
    else:
        begin = 0
        end = 1
    
    extr_type = ex.ExtrapolatedType.LINEAR
 #######################################################


    spl = sp.Splines(x,y, spline_number, coeff_number,
discrete_size, begin, end, last_spline_number, last_coeff_number, extr_type, curve_type)

    all_coeff_matrix = spl.calculateCoefficients()
    writeCoefficients(spl, all_coeff_matrix, name)
    all_splines = spl.createSplines()
   
    
################################################ MAKE PLOTS    
    temp = all_splines[0][0]
    splines =all_splines[0][1]
    plt.figure()
    plt.plot(temp, splines, "b", linewidth=2)


    if last_spline_number >1 and curve_type == sp.CurveType.SATURATION1:
        last_temp = all_splines[1][0]
        last_splines = all_splines[1][1]
        plt.plot(last_temp, last_splines, "c", linewidth=2)


    if extr_type != ex.ExtrapolatedType.NEITHER and curve_type == sp.CurveType.SATURATION1:
        extr_temp =all_splines[-1][0]
        extr_spline = all_splines[-1][1]
        plt.plot(extr_temp, extr_spline, "r", linewidth=2)

   
    plt.plot(x, y, 'bo')
    plt.grid(True)
    plt.savefig(full_name + '.png')
    plt.show()
    plt.close()
    
def interStep(x, y, curve_type, full_name, name, settings ):
    x = x[1::]
    y = y[1::]
    # # print(df)
    print(x)
    a_coeff = [ (y[i+1] - y[i])/(x[i+1] - x[i]) for i in range(0, x.size-1)]
    b_coeff = [y[i] - a_coeff[i]*x[i] for i in range(0, x.size-1)]
    print("a coeff\n", len(a_coeff))
    print("b coeff\n", len(b_coeff))
    coeff_vec = [[ b_coeff[i], a_coeff[i] ] for i in range (len(a_coeff))]
    print("coeff_vec\n", coeff_vec)

    print(len(coeff_vec))
    ##### THIS IS WRITE FILE 
#     with open("spline_coeff"+ name + ".txt", "w+") as f:
#             np.savetxt(f, np.array(([x.size - 1]) ) )
#             for i in range ( len(coeff_vec) - 1, -1, -1 ):
#                     np.savetxt(f, np.array(( [ len(coeff_vec[i]) ] ) ) )
#                     np.savetxt(f, np.array(([x[i+1], x[i] ]) ) )
#                     np.savetxt( f, coeff_vec[i] )
                    
#     temp, splines = self.makeSplines(self.coeff_matrix, x_board, size, left, right)           
#         all_splines = [[temp, splines]] 
    splines =[]
    X = []
    
    for i in range(x.size-1):
        temp = np.linspace(x[i], x[i+1], 1000)
        one = coeff_vec[i][0] + coeff_vec[i][1] * temp
        splines.extend(one)
        X.extend(temp)
        
    plt.figure()    
#     if curve_type == sp.CurveType.SATURATION1:
    plt.plot(x, y, 'bo')
    plt.plot(X, splines, 'c')
    plt.grid(True)
        #     plt.savefig(full_name + '.png')
    plt.show()
    plt.close()


#     spl = sp.Splines(x,y, spline_number, coeff_number,
# discrete_size, begin, end, last_spline_number, last_coeff_number, extr_type, curve_type)

#     all_coeff_matrix = spl.calculateCoefficients()
#     writeCoefficients(spl, all_coeff_matrix, name)
#     all_splines = spl.createSplines()
   
    


# In[ ]:


def makeCurves(excel_data, prefix, type_name,  ind, settings):
    for i in range(3):
        print("iteration ", i)
        if i==0:
            x=excel_data[:,2]
            y=excel_data[:,0]
            curve_type = sp.CurveType.SATURATION1
            name = type_name + 'Pc_'+ ind

        elif i==1:
            x=excel_data[:,2]
            y=excel_data[:,4] 
            curve_type = sp.CurveType.PERMEABILITY1
            name = type_name +  'Kw_'+ ind

        else:
            x=excel_data[:,2]
            y=excel_data[:,3]
            curve_type = sp.CurveType.PERMEABILITY2
            name = type_name + 'Kn_' + ind
           
        full_name =prefix + name
        splineStep(x,y, curve_type, full_name, name, settings)
#         interStep(x,y, curve_type, full_name, name, settings)


# In[ ]:


# ########## CREATE FILE

def writeCoefficients(spl, coeff_list, name):
#     print("x_board\n", spl.x_board)
#     print("last_x_board\n", spl.last_x_board)
    if (spl.curve_type == sp.CurveType.SATURATION1 and spl.last_spline_number != 1):
        boards = np.concatenate( (spl.x_board[:-2], spl.last_x_board ), axis=0)
        if spl.extr_type !=ex.ExtrapolatedType.NEITHER:
            boards = np.append(boards,[0])
    else:
        boards = spl.x_board
    
#     print("coefficient matrix\n" , len(all_coeff_matrix) )
#     print("coefficient matrix\n" , (coeff_list) )
#     print("boards\n", boards)
    with open("Data/coefficients/" + name + ".txt", "w+") as f:
        np.savetxt(f, np.array(([boards.size - 1]) ) )
        if (spl.curve_type != sp.CurveType.PERMEABILITY1):
            
            for i in range ( len(coeff_list) - 1, -1, -1 ):
                np.savetxt(f, np.array(( [coeff_list[i].size] ) ) )
                np.savetxt(f, np.array(([boards[i+1], boards[i] ]) ) )
                np.savetxt( f, coeff_list[i] )
        else:
            for i in range ( len(coeff_list) ):
                np.savetxt(f, np.array(( [coeff_list[i].size] ) ) )
                np.savetxt(f, np.array(([boards[i+1], boards[i] ]) ) )
                np.savetxt( f, coeff_list[i] )
            
#             f.write("------------------\n")
#         f.write("------------------\n")


def checkSeparatedTable(elem, dirname):
    value = 0
    with open(dirname + "/" + elem) as file:
        line = file.readlines()
        for i in range(len(line)-1):
            if line[i]=="\n" and line[i+1] !="\n":
                value = i
#                 print (value)
    return value

