
## 1 Score matrix
#To find the local alignment of b with a the Smith-Waterman calculates a scoring matrix first. 
#The following code calculates this matrix for two strings a and b with linear gap costs.

import numpy as np
import itertools

def matrix(a,b, match_score = 2, gap_cost = 2): #a = column length, b = row length 
    H = np.zeros((len(a)+1,len(b)+1), np.int)

    for i,j in itertools.product(range(1,H.shape[0]), range(1, H.shape[1])): #product is equivalent to nested for loop
        #we skip i=0 nd j=0 because they are not important for thi algorithm

        match = H[i-1,j-1] + (match_score if a[i-1] == b[i-1] else -match_score) 
        delete = H[i-1,j] - gap_cost
        insert = H[i,j-1] -gap_cost
        H[i,j] = max(match, delete, insert)

    return H

##2 backtracing 
# From the calculated scoring matrix to calculate the optimal alignment of b with a. Since b will not simply
# be a substring of a in most cases, some version of b that includes gaps and deletions needs to be calculated
def traceback(H, b, b_="", old_i = 0):
    H_flip = np.flip(H, 0) #up is down and vice versa
    H_flip = np.flip(H_flip,1) #left becomes right and vice versa
    #same as H_flip = np.flip(np.flip(H,0),1) for the complete reversion of the matrix

    i_,j_ = np.unravel_index(H_flip.argmax(), H_flip.shape) #i_ = down, j = right
    # Converts a flat index or array of flat indices into a tuple of coordinate arrays.

    i,j = np.subtract(H.shape, (i_ + 1, j_ + 1))#subtra

    if H[i, j] == 0:
        return b_, j
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_

    return traceback(H[0:i, 0:j], b, b_, i)

def smith_waterman(a, b, match_score=3, gap_cost=2):
    a, b = a.upper(), b.upper()
    H = matrix(a, b, match_score, gap_cost)
    b_, pos = traceback(H, b)
    return pos, pos + len(b_)




seq1 = "CANE"
seq2 = "CERVO"
m = matrix(seq1, seq2)
print(m)
print("")
c,j_= traceback(m, seq2)
print(c)
print(j_)