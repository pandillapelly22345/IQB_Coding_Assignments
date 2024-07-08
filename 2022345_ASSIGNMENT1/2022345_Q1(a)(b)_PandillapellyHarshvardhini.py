import numpy as np
import pandas as pd

def nw(x, y, match=2, mismatch=-3, gap=-1):
    nx = len(x)
    ny = len(y)
    
    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))
    """
    In the below two lines there was an error i.e '-' in front of gap. So I modified np.linspace(0, -gap * nx, nx + 1) to np.linspace(0, gap * nx, nx + 1)
    Previously, the expression had negative values due to the -gap term, resulting in positive values from the multiplication of two negative numbers. 
    However, since the gap is given as negative, there was no need for the negation. Therefore, I adjusted the expression to reflect this, 
    ensuring that the gaps are correctly represented.
    """
    '''
    Filling the first row and first column with the gap penality.
    '''
    F[:, 0] = np.linspace(0, gap * nx, nx + 1)
    F[0, :] = np.linspace(0, gap * ny, ny + 1)

    # Pointers to trace through an optimal alignment.
    """
    Previously, the P matrix was initialized after filling the F matrix. However, since the P matrix is used for tracing back the alignments, 
    it needs to be filled simultaneously with the F matrix. Therefore, it should be declared and initialized at the same time as the F matrix.
    """
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3
    P[0, :] = 4

    # Matrix filling.
    '''
    F matrix filling: 

    1) The matrix filling begins from cell (1,1) after the first row and first column are initialized with gap penalties. 
       Each cell in matrix F represents a match, mismatch, or gap in the sequences.

    2) For a match, the score is obtained from the diagonal cell (i-1, j-1) and the match score is added. 
       For a mismatch, the mismatch score is added. Match or mismatch can be checked using i-1 and j-1,
       if they have same charecterthen match, otherwise mismatch.

    3) If a cell represents a gap in the second sequence, the score is obtained from the 
       upper cell (i-1, j) and the gap score is added.

    4) If a cell represents a gap in the first sequence, the score is obtained from the 
       left cell (i, j-1) and the gap score is added.

    5) Finally, each cell in matrix F(i,j) is filled with the maximum score among these three possibilities
       from point 2, 3 and 4.

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
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            match_score = F[i - 1, j - 1] + (match if x[i - 1] == y[j - 1] else mismatch)
            up_score = F[i - 1, j] + gap
            left_score = F[i, j - 1] + gap
            F[i, j] = max(match_score, up_score, left_score)
            
            if F[i, j] == match_score:
                P[i, j] += 2
            if F[i, j] == up_score:
                P[i, j] += 3
            if F[i, j] == left_score:
                P[i, j] += 4

    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)
    print()

    # Trace through an optimal alignment.
    '''
    Explanation of traceback logic:
    - We start the traceback process from the bottom-right corner of the F matrix, which corresponds to the end of the sequences.
    - We initialize two empty lists, rx and ry, to store the aligned sequences.
    - The while loop continues until we reach either the top or left boundary of the matrix, ensuring that we trace back the entire alignment.
    - Within the loop, we examine the value of the traceback matrix P[i, j] to determine the direction of the traceback.
    - If P[i, j] indicates a match or mismatch (2, 5, 6, 9), we append characters from both sequences to rx and ry,
      decrementing i and j to move diagonally towards the top-left corner of the matrix.
    - If P[i, j] indicates a gap in sequence x (3, 5, 7, 9), we append a character from sequence x and a gap to rx and ry,
      decrementing only i to move upwards in the matrix.
    - If P[i, j] indicates a gap in sequence y (4, 6, 7, 9), we append a gap and a character from sequence y to rx and ry,
      decrementing only j to move leftwards in the matrix.
    - After completing the traceback, we reverse the rx and ry strings to obtain the correct order of characters.
    - Finally, we return the aligned sequences as a single string joined by a newline character.
    '''

    i = nx
    j = ny
    rx = []
    ry = []

    while i > 0 or j > 0:
        if P[i, j] in [2, 5, 6, 9]:
            rx.append(x[i - 1])
            ry.append(y[j-1])
            i -= 1
            j-=1
        elif P[i, j] in [3, 5, 7, 9]:
            rx.append(x[i - 1])
            ry.append('-')
            i -= 1
        elif P[i, j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j - 1])
            j -= 1

    # Reverse the strings.
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]

    # Calculate score
    score_val = score(rx, ry, match, mismatch, gap)

    return '\n'.join([rx, ry]) + f"\nScore: {score_val}"

'''
This function calculates the alignment score by iterating over the sequences. 
It identifies matches, mismatches, and gaps, adding corresponding scores passed as arguments. 
Finally, it returns the total score after iterating over both sequences.
'''
def score(rx, ry, match, mismatch, gap):
    score = 0
    for i in range(0, len(rx)):
        if (rx[i] == ry[i]):
            score += match
        elif rx[i] == '-' or ry[i] == '-':
            score += gap
        else:
            score += mismatch
    return score

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA" 

print(nw(seq1, seq2))