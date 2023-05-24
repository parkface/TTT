import numpy as np
import math
from pdb import set_trace

import matplotlib.pyplot as plt



def make_adjacency_matrix(N=1024, scheme="singleton", save=True):
    assert (N&(N-1) == 0) and (N>=4), f"{N} should be a powers of two and larger than 2"
    assert scheme in ["singleton", "normal"], f"{scheme} is not supported" 
    # make adjacency matrix of BUs
    # N: poly length
    adj_mat = np.zeros((N//2,N//2),dtype=np.int32)
    # adj_mat[i,j] = 1 if i-th BU takes input from j-th BU or gives output to j-th BU 
    for i in range(N//2):
        # i-th BU gives outputs to j-th BUs
        first_output_BU_idx = i//2
        second_output_BU_idx = i//2 + N//4
        # i-th BU takes inputs from j-th BUs
        if 2*i < N//2:
            first_input_BU_idx = 2*i
            second_input_BU_idx = 2*i+1
        else:
            first_input_BU_idx = 2*i-N//2
            second_input_BU_idx = 2*i+1-N//2
        adj_mat[i,first_output_BU_idx] = 1
        adj_mat[i,second_output_BU_idx] = 1
        adj_mat[i,first_input_BU_idx] = 1
        adj_mat[i,second_input_BU_idx] = 1
    
    if save == True:
        np.savetxt(f"BU_adjacency_matrix_using_{scheme}_for_{N}length_poly.csv", adj_mat, delimiter=",", fmt="%d")
    return adj_mat

def convert_adjacency_matrix(adj_mat, N=1024, scheme = "half_fold",save=True):
    # i-th BU and i+n//2 -th BU is located at the same PE, where adj_mat is a nxn matrix
    n, _= adj_mat.shape

    new_mat = np.zeros((n//2, n//2), dtype=np.int32)
    
    if scheme=="half_fold":
        for i in range(n):
            for j in range(n):
                if (i >= n//2) and (j >= n//2):
                    new_mat[i-n//2, j-n//2] = adj_mat[i,j] or adj_mat[i-n//2,j-n//2]
                elif (i >= n//2) and (j < n//2):
                    new_mat[i-n//2, j] = adj_mat[i,j] or adj_mat[i-n//2,j]
                elif (i < n//2) and (j >= n//2):
                    new_mat[i, j-n//2] = adj_mat[i,j] or adj_mat[i,j-n//2]
                else:
                    new_mat[i,j] = adj_mat[i,j]
    elif scheme=="neighbor_merge":
        for i in range(n//2):
            for j in range(n//2):
                new_mat[i,j] = adj_mat[2*i,2*j] or adj_mat[2*i+1,2*j] or adj_mat[2*i,2*j+1] or adj_mat[2*i+1,2*j+1]

    if save == True:
        np.savetxt(f"BU_adjacency_matrix_using_{scheme}_for_{N}length_poly_half_folded.csv", adj_mat, delimiter=",", fmt="%d")
    return new_mat

def main():
    N = 64
    #scheme = "half_fold"
    FFT = "singleton"
    scheme = "neighbor_merge"
    adj_mat = make_adjacency_matrix(N)
    new_mat = convert_adjacency_matrix(adj_mat, N, scheme) 
    plt.matshow(adj_mat)
    plt.savefig(f"BU_adjacency_matrix_using_{FFT}_for_{N}length_poly.png", bbox_inches="tight")
    plt.matshow(new_mat)
    plt.savefig(f"BU_adjacency_matrix_using_{FFT}_for_{N}length_poly_{scheme}.png", bbox_inches="tight")
    plt.show()


if __name__=="__main__":
    main()
    
