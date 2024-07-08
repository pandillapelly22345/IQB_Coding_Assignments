import numpy as np
import pandas as pd


'''
This function is same as nw() function in Question1(a)_(b).py, populates the scoring matrix using dynamic programming (DP). 
But in this after matrix population, it invokes the traceback() function to trace all optimal alignments.
'''
def nw(x, y, match=2, mismatch=-3, gap=-1):
    nx = len(x)
    ny = len(y)
    
    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, gap * nx, nx + 1)
    F[0, :] = np.linspace(0, gap * ny, ny + 1)

    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3
    P[0, :] = 4

    # Matrix filling.
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            match_score = F[i - 1, j - 1] + (match if x[i - 1] == y[j - 1] else mismatch)
            left_score = F[i - 1, j] + gap
            up_score = F[i, j - 1] + gap
            F[i, j] = max(match_score, left_score, up_score)
            
            if F[i, j] == match_score:
                P[i, j] += 2
            if F[i, j] == left_score:
                P[i, j] += 3
            if F[i, j] == up_score:
                P[i, j] += 4

    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)

    # Trace through optimal alignments.
    alignments = []
    traceback([], [], nx, ny, x, y, P, alignments)
    
    return alignments
'''
This function recursively finds all the optimal alignments, following the 
same tracing-back approach as in the previous implementation (Q1(a).py).

In this we handling cells with multiple incoming paths. So for that we are using recursion.

The base case occurs when both i and j become less than or equal to zero, indicating the end of matrix traversal. 
At this point, we store the current optimal alignment into string rx_str for i(rows) and ry_str j(columns) then we 
append this to alignment list to show the alignments further.
'''
def traceback(rx, ry, i, j, x, y, P, alignments):
    if i <= 0 and j <= 0:
        rx_str = ''.join(rx)[::-1]
        ry_str = ''.join(ry)[::-1]
        alignments.append((rx_str, ry_str))
        return
    
    if P[i, j] in [2, 5, 6, 9]:
        p = rx + [x[i - 1]]
        q = ry + [y[j - 1]]
        traceback(p, q, i - 1, j - 1, x, y, P, alignments)
    if P[i, j] in [3, 5, 7, 9]:
        p = rx + [x[i - 1]]
        q = ry + ['-']
        traceback(p, q, i - 1, j, x, y, P, alignments)
    if P[i, j] in [4, 6, 7, 9]:
        p = rx + ['-']
        q = ry + [y[j - 1]]
        traceback(p, q, i, j - 1, x, y, P, alignments)


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

alignments = nw(seq1, seq2)
max_score = []
print("\nOptimal Alignments:")
for alignment in alignments:
    rx, ry = alignment
    print(rx)
    print(ry)
    a = score(rx, ry, 2, -3, -1) # Score calculation
    max_score.append(a)
    print("score:", a)  
    print()

# Find the maximum score
max_value = max(max_score)

# Find indices of all alignments with the maximum score
max_indices = [i for i, score in enumerate(max_score) if score == max_value]

# Print the best alignments
'''
It will print the best alignment among the optimal alignments. For example if there the 6 optimal alignment from
which 4 alignments score is 8 and remaining 2 alignment score is 3 then it will print the best alignment(s) i.e.
4 alignment with their score in this case.
'''
print("\nBest Alignment(s):")
for index in max_indices:
    print("Alignment:", alignments[index])
    print("Score:", max_value)
    print()
