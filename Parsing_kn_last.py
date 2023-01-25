#!/usr/bin/env python
# coding: utf-8

import Splines as sp
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
    begin = 1
    end = np.min(x) * 0.99 if curve_type == sp.CurveType.SATURATION1 else 0
    
    extr_type = sp.ExtrapolatedType.LINEAR
 #######################################################


    spl = sp.Splines(x,y, spline_number, coeff_number,
discrete_size, begin, end, last_spline_number, last_coeff_number, extr_type, curve_type)

    all_coeff_matrix = spl.calculateAllCoefficients()
    writeCoefficients(spl, all_coeff_matrix, name)
    all_splines = spl.createSplines()
   
    
################################################ MAKE PLOTS    
    temp = all_splines[0][0]
    splines =all_splines[0][1]
    plt.figure()
    plt.plot(temp, splines, "b", linewidth=2)


    if spl.last_spline_number >1:
        last_temp = all_splines[1][0]
        last_splines = all_splines[1][1]
        plt.plot(last_temp, last_splines, "c", linewidth=2)


    if extr_type != sp.ExtrapolatedType.NEITHER and curve_type == sp.CurveType.SATURATION1:
        extr_temp =all_splines[-1][0]
        extr_spline = all_splines[-1][1]
        plt.plot(extr_temp, extr_spline, "r", linewidth=2)

   
    plt.plot(x, y, 'bo')
    plt.grid(True)
    plt.savefig(full_name + '.png')
#     plt.show()
    plt.close()


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


# In[ ]:


# ########## CREATE FILE

def writeCoefficients(spl, coeff_list, name):
#     print("x_board\n", spl.x_board)
#     print("last_x_board\n", spl.last_x_board)
    if (spl.curve_type == sp.CurveType.SATURATION1 and spl.last_spline_number != 1):
        boards = np.concatenate( (spl.x_board[:-2], spl.last_x_board ), axis=0)
        if spl.extr_type !=sp.ExtrapolatedType.NEITHER:
            boards = np.append(boards,[0])
    else:
        boards = spl.x_board
    
#     print("coefficient matrix\n" , len(all_coeff_matrix) )
#     print("coefficient matrix\n" , (coeff_list) )
#     print("boards\n", boards)

#     with open("Data/coefficients/" + name + ".txt", "w+") as f:
#         np.savetxt(f, np.array(([boards.size - 1]) ) )
#         for i in range ( len(coeff_list) - 1, -1, -1 ):
#             np.savetxt(f, np.array(( [coeff_list[i].size] ) ) )
#             np.savetxt(f, np.array(([boards[i+1], boards[i] ]) ) )
#             np.savetxt( f, coeff_list[i] )
# #             f.write("------------------\n")
# #         f.write("------------------\n")


def checkSeparatedTable(elem, dirname):
    value = 0
    with open(dirname + "/" + elem) as file:
        line = file.readlines()
        for i in range(len(line)-1):
            if line[i]=="\n" and line[i+1] !="\n":
                value = i
#                 print (value)
    return value

