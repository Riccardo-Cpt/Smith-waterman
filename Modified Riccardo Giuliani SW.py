#Riccardo Giuliani 213785

# When a bifurcation is found, that means that there are two possible ways to follow, this algorithm chooses the diagonal one. 
# It discards the up and/or left direction

#The algorithm allows to to print alignments in decreasing or increasing oreder of gaps depending on the sort parameter we choose. For this project sort_type is set = 2 to order by increase number of gaps.
#Length and max_gaps filters are set to 0 and inf respectively beacause they are not requested

import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "Smith Waterman")
parser.add_argument("sequence1", type = str, help = "First sequence for alignment")
parser.add_argument("sequence2", type = str, help = "Second sequence for alignment")
parser.add_argument("match", type = int, help = "Match score")
parser.add_argument("mismatch", type = int, help = "Mismatch score")
parser.add_argument("gap", type = int, help = "Gap penalty")
#parser.add_argument("minimumscore", type = int, help = "Minimum score")
#parser.add_argument("minimumlength", type = int, help = "Minimum length")
#parser.add_argument("maximumgaps", type = int, help = "Maximum number of gaps desired")  It can be changed at the end of the script


#substitution matrix
def sub_matrix(match, mismatch, gap):
    sM = {"match_score" : match, "mismatch_score": mismatch, "gap_penalty": gap}
    return sM
    
# scoring matrix
def scor_mat(sequence_1, sequence_2, sM):
    n = len(sequence_1)
    m = len(sequence_2)
    M = np.zeros((m+1,n+1)) #rows = seq2, columns = seq1
        
    for i in range(1, M.shape[0]):
            for j in range(1, M.shape[1]):
                
                if sequence_2[i-1] == sequence_1[j-1]:
                    align = sM["match_score"]
                    
                else:
                    align = sM["mismatch_score"]
                
                M[i,j] = max(M[i-1,j-1] + align , M[i-1, j] + sM["gap_penalty"],
                              M[i, j-1] + sM["gap_penalty"], 0)
                
    return M    

def print_matrix(M, seq1, seq2):
    firstrow = list(map(str, M[0].astype(int)))
    info = "\t-\t" + "\t".join(list(seq1)) + "\n" + "-\t" + "\t".join(firstrow) + "\n"
    for i in range(len(seq2)):
        row = list(map(str, M[i + 1].astype(int)))
        info = info + seq2[i] + "\t" + "\t".join(row) + "\n"
    print(info)


#Max finder
def findMax(M):
    val = np.amax(M)
    i,j = np.where(M == val)
    return i,j

#Backtracking
def traceBack(M,sM, seq1, seq2, i_,j_):
    al1, al2, n_gap = [],[],[]
    
    for k in range(len(i_)):
        
        i = i_[k]; j = j_[k]
        max_val = M[i,j]
        align1, align2 = "",""
        gaps = 0
        
        while M[i,j] > 0:
            if max_val == 0:
                print("No complementarity between sequences")

            if M[i-1,j-1] + sM["match_score"] == M[i,j] and seq1[j-1] == seq2[i-1] \
            or M[i-1, j-1] + sM["mismatch_score"] == M[i,j]:
                align1 = align1 + seq1[j-1]; align2 = align2 + seq2[i-1]
                i-=1; j-=1

            elif M[i-1,j]+sM["gap_penalty"] == M[i,j]:
                align1 = align1+"-"
                align2 = align2+seq2[i-1]
                i-=1
                gaps +=1    

            elif M[i,j-1] + sM["gap_penalty"] == M[i,j]:
                align1 = align1 + seq1[j-1]
                align2 = align2 + "-"
                j-=1
                gaps +=1

            else:
                print("Error: not valid alignament")
        
        al1.append(align1[::-1])
        al2.append(align2[::-1])
        n_gap.append(gaps)
        
    return al1, al2, n_gap
                        
#Filtering
def FilterAligns(M,sM, seq1, seq2, min_length, min_score, max_gaps):
    L = M.flatten()
    L = np.unique(L)
    L = np.sort(L)

    alignments = []
   
    for score in L[::-1]:
        if score >= min_score:
            i_,j_ = np.where(M == score)
            a1, a2, gaps = traceBack(M,sM,seq1, seq2, i_,j_)
            
            for k in range(len(a1)):
                if len(a1[k]) >= min_length and gaps[k] < max_gaps:                   
                    alignments.append({"score" : score, "align1" : a1[k], "align2" : a2[k],
                                            "length" : len(a1[k]), "gaps": gaps[k]})
            
    return alignments
        

def sort_alignments(alignments, sort_type):
    if sort_type == 1: # sort by the number of gaps in decreasing order
        sort = sorted(alignments, key = lambda x: x["gaps"], reverse= True)
        return sort
    elif sort_type == 2: # sort by the number of gaps in increasing order
        sort = sorted(alignments, key = lambda x: x["gaps"], reverse= False)
        return sort


def printscores(filter_sorted_align):
    for align in filter_sorted_align:
        matches, mismatches, gaps = 0,0,0
        align1 = align["align1"]
        align2 = align["align2"]
        align_graphic = []

        for i, j in zip(align1, align2):
            if i == j:
                align_graphic.append('|') 
                matches += 1
            elif '-' in (i, j):
                align_graphic.append(' ')
                gaps += 1
            else:
                align_graphic.append('x') 
                mismatches += 1
        print("{}\n{}\n{}\n{} matches, {} mismatches and {} gaps. Length: {} Score: {} \n".format(align2,
         "".join(align_graphic), align1, matches, mismatches, gaps, align["length"], align["score"]))



def SmithWaterman(seq1, seq2, match = 3, mismatch = -3, gap = -2 ):

    sM = sub_matrix(match,mismatch,gap)
    M = scor_mat(seq1,seq2, sM)
    print_matrix(M, seq1, seq2)
    #For score and length the filter is set >=n, if I want more than 6 I put 7
    High_scores = FilterAligns(M,sM,seq1,seq2,min_length = 0,min_score = 7,max_gaps = float('inf') )
    sort = sort_alignments(High_scores, 2)
    printscores(sort)

# Argparse definitions
args = parser.parse_args()
print(args)
seq2 = args.sequence1
seq1 = args.sequence2
match = args.match
mismatch = args.mismatch
gap = args.gap

if __name__ == "__main__":
    SmithWaterman(seq1, seq2, match, mismatch, gap)
