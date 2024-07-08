import numpy as np
import pandas as pd

'''
This function is same as nw() function in Question1(c)_(d).py, populates the scoring matrix using dynamic programming (DP). 
But in this after matrix population, It is also doin the trace back to print the alignment.
'''

def nw(x, y, match=2, mismatch=-1, gap=-3):
    nx = len(x)
    ny = len(y)
    
    '''
    In local alignment, every cell has a minimum value of 0. Therefore, we avoid 
    penalizing the first row and first column with a gap penalty.
    '''
    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))
    
    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3
    P[0, :] = 4

    # Matrix filling.
    '''
    The logic is same here as in the global alignment F matrix filling.
    F matrix filling:
    1) Matrix F is filled starting from cell (1,1), indicating matches, mismatches, and gaps in both sequences.
    2) Matches fetch their score from the diagonal cell, while mismatches add 
       their respective penalty from diagnal cell.
    3) Gaps in the second sequence derive their score from the cell above, with an added penalty.
    4) Gaps in the first sequence fetch their score from the left cell, along with a penalty.
    5) Each cell in F receives the maximum score among the three options, considering a minimum value of 0 for local alignment.

    explained in detail global alignment.

    The logic is same here as in the global alignment F matrix filling
    P matrix filling:

    1) Initially, Matrix P filled the first row with 4 and the first column with 3, indicating values for gaps 
    in the first and second sequences, respectively. This setup suggests 2 for matches/mismatches, 
    as confirmed by later traceback.

    2) After computing F[i,j], we populate P[i,j] based on the filled value. 
    Since F[i,j] can come from multiple cells, we fill P to track multiple paths 
    by adding 2, 3, or 4 depending on the source.

    3) If F[i,j] originates from the diagonal, we assign 2 to P[i,j].
       If F[i,j] originates from the left cell, we assign 3 to P[i,j].
       If F[i,j] originates from the upper cell, we assign 4 to P[i,j].

    '''
    max_score = -1
    max_i = 0
    max_j = 0
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            a = 0
            match_score = F[i - 1, j - 1] + (match if x[i - 1] == y[j - 1] else mismatch)
            up_score = F[i - 1, j] + gap
            left_score = F[i, j - 1] + gap
            F[i, j] = max(0, match_score, up_score, left_score, a)
            
            # Update max_score and the indices for local alignment
            if F[i, j] >= max_score:
                max_score = F[i, j]
                max_i = i
                max_j = j
        
            if F[i, j] == match_score:
                P[i, j] = 2
            elif F[i, j] == up_score:
                P[i, j] = 3
            elif F[i, j] == left_score:
                P[i, j] = 4

    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)
    print()

    # Trace through the optimal local alignment.
    '''
    The traceback process starts from the cell with the maximum score (max_i, max_j).
    We initialize two empty lists, rx and ry, to store the aligned sequences.
    While traversing through the traceback path, we iterate until we reach either the top or left boundary of the matrix
    or encounter a cell with a score of 0, indicating the end of the alignment.
    Within the inner loop, we check the value of the traceback matrix P[i, j].
    Depending on its value, we decide whether to align characters, insert a gap in one of the sequences, or both.
    If P[i, j] indicates a match/mismatch (2, 5, 6, 9), we append characters from both sequences to rx and ry,
    decrement i and j, and move diagonally in the matrix.
    If P[i, j] indicates a gap in sequence x (3, 5, 7, 9), we append a character from sequence x and a gap to rx and ry,
    decrementing only i to move upwards in the matrix.
    If P[i, j] indicates a gap in sequence y (4, 6, 7, 9), we append a gap and a character from sequence y to rx and ry,
    decrementing only j to move leftwards in the matrix.
    After completing the traceback for one alignment path, we reverse the rx and ry strings and append them as a tuple
    to the alignments list.
    Finally, we return the list of alignments containing all the optimal local alignments found.
    '''
    i = max_i
    j = max_j
    alignments = []
    while i > 0 and j > 0 and F[i, j] != 0:
        rx = []
        ry = []
        while i > 0 and j > 0 and P[i, j] != 0:
            if P[i, j] in [2, 5, 6, 9]:
                rx.append(x[i - 1])
                ry.append(y[j - 1])
                i -= 1
                j -= 1
            elif P[i, j] in [3, 5, 7, 9]:
                rx.append(x[i - 1])
                ry.append('-')
                i -= 1
            elif P[i, j] in [4, 6, 7, 9]:
                rx.append('-')
                ry.append(y[j - 1])
                j -= 1
        
        # Reverse the strings and store alignment
        rx = ''.join(rx)[::-1]
        ry = ''.join(ry)[::-1]
        alignments.append((rx, ry))

    return alignments

'''
This function calculates the alignment score by iterating over the sequences. 
It identifies matches, mismatches, and gaps, adding corresponding scores passed as arguments. 
Finally, it returns the total score after iterating over both sequences.
'''
def score(rx, ry, match, mismatch, gap):
    score = 0
    for i in range(len(rx)):
        if rx[i] == ry[i]:
            score += match
        elif rx[i] == '-' or ry[i] == '-':
            score += gap
        else:
            score += mismatch
    return score

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA" 

alignments = nw(seq1, seq2)
for rx, ry in alignments:
    print(rx)
    print(ry)
    a = score(rx, ry, 2, -1, -3) # Score calculation
    print("Score:", a)
