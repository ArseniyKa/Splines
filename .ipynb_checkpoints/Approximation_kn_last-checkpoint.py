{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "e7540fa7",
   "metadata": {},
   "outputs": [],
   "source": [
    "# !pip install pandas\n",
    "# !pip install odfpy\n",
    "import pandas as pd\n",
    "import numpy as np\n",
    "import matplotlib\n",
    "from decimal import Decimal\n",
    "from numpy import linalg as LA\n",
    "\n",
    "matplotlib.use('TkAgg')\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "def simpleFitting(x, y, polinom_size):\n",
    "        vandermond_tr = np.zeros((polinom_size,x.size))\n",
    "        \n",
    "        for i in range(polinom_size):\n",
    "              vandermond_tr[i] = x**i\n",
    "        \n",
    "        vandermond = vandermond_tr.transpose()\n",
    "        \n",
    "        XX = vandermond_tr.dot(vandermond)\n",
    "        XX_inv = np.linalg.inv(XX)\n",
    "        XY = vandermond_tr.dot(y)\n",
    "        coeff_vector = XX_inv.dot(XY)\n",
    "\n",
    "        return coeff_vector\n",
    "    \n",
    "    \n",
    "def conditionalFitting(x, y, xo, yo, dyo, polinom_size, constr_number):\n",
    "    \n",
    "    vandermond_tr = np.zeros((polinom_size,x.size))\n",
    "    lagrange_multypliers_matrix = np.zeros((constr_number, polinom_size))\n",
    "    init_conditial_matrix = np.zeros((constr_number, polinom_size))\n",
    "    function_vec = np.zeros(polinom_size)\n",
    "    derivative_vec = np.zeros(polinom_size)\n",
    "    \n",
    "    for i in range(polinom_size):\n",
    "        function_vec[i] = xo**i\n",
    "        derivative_vec[i] = i * xo**(i-1)\n",
    "        vandermond_tr[i] = x**i\n",
    "\n",
    "############################################### Fill initial conditional and Lagrange matrixes\n",
    "    for j in range(constr_number):\n",
    "        init_conditial_matrix[j] = function_vec if j==0 else derivative_vec\n",
    "        lagrange_multypliers_matrix[j] = -0.5 * ( function_vec if j==0 else derivative_vec )\n",
    "\n",
    "\n",
    " ####################### Simple matrix fitting   \n",
    "\n",
    "    vandermond = vandermond_tr.transpose()        \n",
    "    XX = vandermond_tr.dot(vandermond)\n",
    "    XY = vandermond_tr.dot(y)\n",
    "    XY = np.append(XY, yo)\n",
    "    XY = np.append(XY, dyo) if constr_number == 2 else XY\n",
    "\n",
    "\n",
    "##################### Create matrix for Algebraic system solving\n",
    "\n",
    "    lagrange_multypliers_matrix = lagrange_multypliers_matrix.transpose()\n",
    "    init_conditial_matrix = np.hstack(( init_conditial_matrix, np.zeros((constr_number,constr_number)) ))\n",
    "    \n",
    "    A = np.hstack((XX,lagrange_multypliers_matrix))\n",
    "    A = np.vstack((A, init_conditial_matrix))\n",
    "\n",
    "#     print(lagrange_multypliers_matrix)\n",
    "#     print(\"------\")\n",
    "#     print(A)\n",
    "\n",
    "\n",
    "    coeff_vector = LA.solve(A, XY)\n",
    "    ### to remove Lagrange variables\n",
    "    coeff_vector = coeff_vector[:polinom_size]\n",
    "#     print (\"coeff_vector\", coeff_vector)\n",
    "    return coeff_vector\n",
    "\n",
    "\n",
    "###################################### FUNCTION CURVE\n",
    "def calculateCurveValue(xo, coeff_vector):\n",
    "    # print(\"coeff_vector\", coeff_vector)\n",
    "    value = 0\n",
    "    for n in range(coeff_vector.size):\n",
    "        value += coeff_vector[n] * xo**n\n",
    "\n",
    "    return value\n",
    "########################################### DERIVATIVE CURVE\n",
    "\n",
    "def calculateDerivative(xo, coeff_vector):\n",
    "    derivative = 0\n",
    "    \n",
    "    for n in range(1, coeff_vector.size):\n",
    "        derivative += n* coeff_vector[n] * xo**(n-1)\n",
    "        \n",
    "    return derivative\n",
    "\n",
    "########################################### SECOND DERIVATIVE CURVE\n",
    "\n",
    "def calculateSecondDerivative(x, coeff_vector):\n",
    "    second_derivative = 0\n",
    "\n",
    "    for n in range(2, coeff_number):\n",
    "        second_derivative += n* (n - 1) * coeff_vector[n] * xo**(n-2)\n",
    "        \n",
    "    return second_derivative\n",
    "###############################################################"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "0e64adbd",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.10"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}