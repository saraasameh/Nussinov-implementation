# -*- coding: utf-8 -*-
"""
Created on Sun May 30 18:36:33 2021

@author: pc
"""
import numpy as np

def score (N1 ,N2):
    N1 = str(N1.upper())
    N2 = str(N2.upper())
    score1 =0
    
    if( N1== 'U' and N2 == 'A') or (N1== 'A' and N2 == 'U')or(N1== 'C' and N2 == 'G')or(N1== 'G' and N2 == 'C'):
        score1=1
        
    return score1


def fill (rna):
    N = len(rna)
    #myMatrix = np.zeros((N, N), dtype=int)
    
    myMatrix = np.empty([N, N])
    myMatrix[:] = np.NAN

    myMatrix[range(N), range(N)] = 0 #fill main diagonal
    myMatrix[range(1, N), range( N-1 )] = 0
    
    
    for k in range(N):
        for i, j in zip(range(0, N - k), range(k+1, N)): #fill diagonal diagonal
            
            maximum=0
            if j-i > 1:
                for m in range(i + 1, j):
                    total = myMatrix[i, m] + myMatrix[m + 1, j]
                    if total > maximum:
                        maximum= total
                            
            myMatrix[i][j] = max(
                myMatrix[i + 1, j],
                myMatrix[i, j - 1],
                myMatrix[i + 1, j - 1] + score( rna[i], rna[j]),
                maximum
                )
            
    print(myMatrix)
    return myMatrix


def traceback(i, j, matrix, sequence,pairs):

    if j <= i:
        return

    elif matrix[i,j] == matrix[i,j-1]:
        traceback(i, j-1, matrix, sequence,pairs)
    elif matrix[i,j] == matrix[i+1,j]:
        traceback(i+1, j, matrix, sequence,pairs)
    elif matrix[i,j] == matrix[i+1,j-1] + score(sequence[i], sequence[j]):
        pairs.append((i, j))
        traceback(i + 1, j - 1, matrix, sequence,pairs)
        
    else:
        for k in range(i + 1, j-1):
            if matrix[i, j] == matrix[i,k] + matrix[k + 1,j]:
                traceback(i, k, matrix, sequence,pairs)
                traceback(k+1, j, matrix, sequence,pairs)
                break
                
            
    return pairs



def drawbrackets(seq,pairs):
        
        output = ["."] * len(seq)
        
        for a, b in pairs:
            if a <= b:
                output[a] = "("
                output[b] = ")"
            else:
                output[b] = "(" 
                output[a] = ")"
                
        return "".join(output)


seq='AUGAGGUCAUGCAAU'
filled_matrix = fill(seq)
pairs =[]
t= traceback(0, len(seq)-1, filled_matrix ,seq, pairs)
print(t)
print (drawbrackets(seq, t ))

#CGGACCCAGACUUUC
#GGGAAAUCC

#UAACGUACUGGAGUA
#GGAAUUAGUUAACC