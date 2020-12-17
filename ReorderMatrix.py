###############################################
#                Dingqi Zhang                 #
#             School of Pharmacy              #
#  Macau University of Science and Technology #
#               December, 2020                #
###############################################
import os
import pandas as pd
import numpy as np
import math
from scipy.stats import pearsonr
import operator

curdir = os.path.dirname(os.path.realpath(__file__))
os.chdir(curdir)

matrix = pd.read_excel("matrix.xlsx", index_col=0)

def calculate_potential(value, x, y):
    distance = abs(x-y)/math.sqrt(2)
    mass = 100-value
    return mass*distance

def matrix_energy(matrix):
    size = matrix.shape
    if size[0] != size[1]:
        raise Exception("matrix is not a square!")
    size = size[0]
    energy = 0
    for x in range(size):
        for y in range(size):
            value = float(matrix.iat[x,y])
            energy += calculate_potential(value,x,y)
    return energy

def reordercols(matrix):

    if matrix.shape[0] != matrix.shape[1]:
        raise Exception("Number of rows does not equal number of columns!")

    compound_list = matrix.columns.to_list()

    col_order = compound_list
    
    for colnum in range(len(compound_list)-1):
        selected = np.array(matrix.iloc[:,colnum])
        pearsons = {}
        for other in range(colnum+1,len(compound_list)):
            the_other = np.array(matrix.iloc[:,other])
            pearson = pearsonr(selected, the_other)
            pearsons[other] = pearson[0]
        most_cor_col = max(pearsons.items(), key=operator.itemgetter(1))[0]
        if pearsons[most_cor_col]>0.5:
            col_order_remaining = col_order[colnum+1:]
            col_order_remaining.pop(most_cor_col-colnum-1)
            col_order = col_order[:colnum+1]+[col_order[most_cor_col]]+col_order_remaining
            matrix = matrix[col_order]
    
    return matrix

def main(matrix):
    
    if matrix.shape[0] != matrix.shape[1]:
        raise Exception("Number of rows does not equal number of columns!")

    compound_list = matrix.columns.to_list()

    all_matrices = []
    for compound in range (0,len(compound_list)):
        comp_order = compound_list
        comp_order[0], comp_order[compound] = comp_order[compound], comp_order[0]

        col_reordered_mat = reordercols(matrix[comp_order])
        row_reordered_mat = reordercols(matrix.transpose()[comp_order]).transpose()
        all_matrices.append(col_reordered_mat.columns.to_list())
        all_matrices.append(row_reordered_mat.columns.to_list())
    
    record = open("matrix_permutations.txt","w",encoding="utf-8")
    min_energy_mat = (matrix.columns.to_list(),matrix_energy(matrix))
    record.write("\t".join(str(x) for x in range(0,len(compound_list)))+"\t"+"matrix energy\n")
    for reordered_mat in all_matrices:
        energy = matrix_energy(matrix[reordered_mat].reindex(reordered_mat))
        record.write("\t".join(reordered_mat)+"\t"+str(energy)+"\n")
        if energy < min_energy_mat[1]:
            min_energy_mat = (reordered_mat, energy)
    
    print(min_energy_mat)

    matrix = matrix[min_energy_mat[0]].reindex(min_energy_mat[0])

    print(matrix)
    matrix.to_excel("min_energy_matrix.xlsx")
    record.close()

main(matrix)