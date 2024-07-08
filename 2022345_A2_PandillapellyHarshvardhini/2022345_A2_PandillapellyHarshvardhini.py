# Pandillapelly Harshvardhini(2022345)

protein_sequence = "MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAPSTVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPAKSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKTHMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFKCQICPYASRNSSQLTVHLRSHTASELDDDVPKANCLSTESTDTPKAPVITLPSEAREQMATLGERTFNCCYPGCHFKTVHGMKDLDRHLRIHTGDKPHKCEFCDKCFSRKDNLTMHMRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDERPYKCQLCPYASRNSSQLTVHLRSHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCDVRCTMKANLKSHIRIKHTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRVHSRVHCKDRPFKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVTQRAFRCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPLEPSQDL"

alpha_helix = {
    'E': 1.53, 
    'A': 1.45, 
    'L': 1.34, 
    'H': 1.24, 
    'M': 1.20, 
    'Q': 1.17, 
    'W': 1.14, 
    'V': 1.14, 
    'F': 1.12, 
    'K': 1.07, 
    'I': 1.00, 
    'D': 0.98, 
    'T': 0.82, 
    'S': 0.79, 
    'R': 0.79,
    'C': 0.77, 
    'N': 0.73, 
    'Y': 0.61, 
    'P': 0.59,
    'G': 0.53
}

beta_strand = {
    'M': 1.67,
    'V': 1.65,
    'I': 1.60,
    'C': 1.30,
    'Y': 1.29,
    'F': 1.28,
    'Q': 1.23,
    'L': 1.22,
    'T': 1.20,
    'W': 1.19,
    'A': 0.97,
    'R': 0.90,
    'G': 0.81,
    'D': 0.80,
    'K': 0.74,
    'S': 0.72,
    'H': 0.71,
    'N': 0.65,
    'P': 0.62,
    'E': 0.26
}

"""
Name of amino acid with its symbol:
'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 
'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser',
'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly'
"""


""" 
generateValidWindows() function generates valid windows of a size of 6 based on certain criteria.
It takes three parameters: window_size (size of the window), helixORstrand (whether it's for helix or beta strand),
and window_score (minimum score required for a window to be considered valid).It iterates over the protein sequence,
considers each window of the specified size, calculates the score of each window, and if the score meets the criteria i.e 
should be greater than 4 for helix and 3 for beta (window_score), it adds the window to the list of valid windows.
"""
def generateValidWindows(window_size, helixORstrand, window_score):
    valid_window_range = []
    current_index = 0
    for si in range(0, len(protein_sequence)):
        si = current_index
        ei = si + window_size
        
        if ei > len(protein_sequence):
            break
        
        window_sequence = protein_sequence[si:ei]
        
        counter = sum(1 for char in window_sequence if (alpha_helix[char] if helixORstrand == 1 else beta_strand[char]) >= 1)
        
        if counter >= window_score:
            valid_window_range.append([si, ei - 1])
        
        current_index += 1
        
    return valid_window_range

"""
calculateScore() This function calculates the score of a given subsequence of the protein sequence. 
It takes three parameters: rs (start index of the subsequence), re (end index of the subsequence), and 
helixORstrand (whether it's for helix or beta strand). 
It iterates over the particular given subsequence, retrieves the score of each character from the respective dictionary
(alpha_helix or beta_strand).
"""

def calculateScore(rs,re,helixORstrand):
    score=0
    subsequence=protein_sequence[rs:re+1]
    for char in subsequence:
        if(helixORstrand==1):
            scr=alpha_helix[char]
        else:
            scr=beta_strand[char]
        score=score+scr
    return score


"""
isExtendable() This function checks if a given subsequence (defined by its start and end indices) is extendable based on a minimum score
threshold. It calls the calculateScore function to determine the score of the subsequence and compares it with the threshold. 
If the score is greater than or equal to the threshold, it returns True, indicating that the subsequence is extendable
"""
def isExtendable(rs,re,helixORstrand):
    score = calculateScore(rs,re,helixORstrand)
    if(round(score,2)<4):
        return False
    else:
        return True
    
'''
left_start = ls
left_end = le
right_start = rs
right_end = re
old_start = os
old_end = oe
end_index = ei
start_index = si
'''

"""
extendWindowLeft() This function extends the window to the left (towards smaller indices) to find the maximum extendable window. 
It takes four parameters: si (start index of the window), ei (end index of the window), extended_window_length 
(length of the extended window), and helixORstrand (whether it's for helix or beta strand). It iterates over the protein 
sequence to find the maximum extendable window towards the left by checking if each subsequence is extendable using the 
isExtendable function.
Once it finds the maximum extendable window, it returns the start and end indices of the extended window.
"""

def extendWindowLeft(si, ei, extended_window_length, helixORstrand):
    os = si
    ls = si
    le = ls + (extended_window_length - 1)
    
    while ls >= 0:
        ls -= 1
        le -= 1
        
        if ls < 0 or not isExtendable(ls, le, helixORstrand):
            os = ls + 1
            # oe=re-1
            ewr = extendWindowRight(si, ei, extended_window_length, helixORstrand)
            ewr[1] = ewr[2] - 1
            break
            
    return [os, ei, ewr[1]]


"""
extendWindowRight() This function extends the window to the right (towards larger indices) to find the maximum extendable window. 
It has the same parameters as extendWindowLeft and similar functionality but extends the window towards the right.
"""
def extendWindowRight(si, ei, extended_window_length, helixORstrand):
    oe = ei
    rs = ei - (extended_window_length - 1)
    re = oe
    
    while re < len(protein_sequence):
        rs += 1
        re += 1
        
        if re >= len(protein_sequence) or not isExtendable(rs, re, helixORstrand):
            oe = re - 1
            break
            
    return [si, oe, re]


"""
extendWindow() This function orchestrates the extension process for all valid windows. It takes three parameters: 
valid_window_range (list of valid windows), extended_window_length (length of the extended window), and helixORstrand 
(whether it's for helix or beta strand). It iterates over each valid window, extends it to the left and right using 
extendWindowLeft and extendWindowRight functions, and adds the extended window to the final list of extended windows.
it checks for each amino acid addig to window that it is extendible or not.
"""
def extendWindow(valid_window_range, extended_window_length, helixORstrand):
    final_extended_window = []
        
    for p in valid_window_range:
        si = p[0]
        ei = p[1]
        
        left_extended_window = extendWindowLeft(si, ei, extended_window_length, helixORstrand)
        right_extended_window = extendWindowLeft(si, ei, extended_window_length, helixORstrand)
        
        final_extended_window.append([left_extended_window[0], right_extended_window[2]])
        
    return final_extended_window


"""
initializeArrays() This function initializes arrays for helix, beta strand, and the final merged sequence. 
It creates empty arrays of the specified length i.e length of protien sequence and returns them.
"""
def initializeArrays(sequence_length):
    h_array = ['X'] * sequence_length
    b_array = ['X'] * sequence_length
    final_array = ['X'] * sequence_length
    return h_array, b_array, final_array


"""
final_Helix_or_beta_sequence() This function creates arrays indicating helix or beta strand sequences based on the extended windows 
it takes extended window valus from extended_helix_windows and extended_beta_strand_windows. 
It iterates over the extended windows and sets the corresponding positions in the array to indicate the presence of helix or 
beta strand in the sequence.
It takes following parameters
extended_windows i.e.extended_helix_windows and extended_beta_strand_windows, window_type i.e. H or S and sequence length i.e. 
length of protien sequence.
"""
def final_Helix_or_beta_sequence(extended_windows, window_type, sequence_length):
    window_array = ['X'] * sequence_length
    for start, end in extended_windows:
        for ind in range(start, end + 1):
            window_array[ind] = window_type
    return window_array


"""
Merged_Helix_Beta() function merges helix and beta strand sequences into a final sequence, resolving conflicts. 
It iterates over the arrays indicating helix and beta strand sequences, resolving conflicts where both helix and 
beta strand are indicated. It returns a merged array indicating the final sequence with resolved conflicts.
It checks whenever there is merge between H and S then this fuction decides chich should come.
"""
def Merged_Helix_Beta(h_array, b_array, sequence_length):
    final_array = ['X'] * sequence_length
    ind = 0
    while ind < sequence_length:
        if h_array[ind] == "H" and b_array[ind] == "X":
            final_array[ind] = "H"
            ind += 1
        elif h_array[ind] == "X" and b_array[ind] == "S":
            final_array[ind] = "S"
            ind += 1
        elif h_array[ind] == "H" and b_array[ind] == "S":
            conflict = []
            while ind < sequence_length and h_array[ind] == "H" and b_array[ind] == "S":
                conflict.append(ind)
                ind += 1
            hel_score = sum(alpha_helix[protein_sequence[index]] for index in conflict)
            beta_score = sum(beta_strand[protein_sequence[index]] for index in conflict)
            if hel_score >= beta_score:
                for index in conflict:
                    final_array[index] = "H"
            else:
                for index in conflict:
                    final_array[index] = "S"
        else:
            ind += 1
    return final_array


"""
printSequenceWithWindows() This function basically helps to visualize the sequence with helix region or beta strand region it basically align extended 
region in protien sequence to its corresponding region in protein sequence.
In this question this function aligns protein sequence with helix array and protein sequence with beta strand array which came
from this function - final_Helix_or_beta_sequence.
it also aligns the final merged beta an helix array to protein sequence.
"""
def printSequenceWithWindows(sequence, window_array):
    for i in range(0, len(sequence), 80):
        sequence_line = list(sequence[i:i + 80])
        window_line = [' '] * 80
        for ind in range(i, min(i + 80, len(sequence))):
            if window_array[ind] != 'X':
                window_line[ind - i] = window_array[ind]
        # print("\n")
        print(''.join(sequence_line))
        print(''.join(window_line))
        print("\n")


"""
printThreeSeq() This function helps to visualize three sequences/array at the same time. Like in QUSTION1(c) part.
"""
def printThreeSeq(sequence, helix_array, strand_array):
    for i in range(0, len(sequence), 80):
        sequence_line = list(sequence[i:i + 80])
        helix_line = [' '] * 80
        strand_line = [' '] * 80
        start_index = i
        for ind in range(start_index, start_index + len(sequence_line)):
            if helix_array[ind] == "H":
                helix_line[ind - start_index] = "H"
            if strand_array[ind] == "S":
                strand_line[ind - start_index] = "S"
        # print("\n")
        print(''.join(sequence_line))
        print(''.join(helix_line))
        print(''.join(strand_line))
        print("\n")

# MAIN FUNCTION

# Generate valid windows for helix of size 6
helix_windows = generateValidWindows(6, 1, 4)

# Generate valid windows for beta strand of size 6
beta_strand_windows = generateValidWindows(5, 2, 3)

# for each helix window(calculated above) it gives extended helix windows
extended_helix_windows = extendWindow(helix_windows, 4, 1)

# for each beta strand window(calculated above) it gives extended beta strand windows
extended_beta_strand_windows = extendWindow(beta_strand_windows, 4, 2)

# Initialize arrays
final_h_array, final_b_array, final_array = initializeArrays(len(protein_sequence))

# Process helix and beta strand windows
final_h_array = final_Helix_or_beta_sequence(extended_helix_windows, 'H', len(protein_sequence))
final_b_array = final_Helix_or_beta_sequence(extended_beta_strand_windows, 'S', len(protein_sequence))

# Merged_Helix_Beta seq
final_array = Merged_Helix_Beta(final_h_array, final_b_array, len(protein_sequence))

# Print sequence with final windows
print("QUESTION 1(a)")
print()
print("Helix Windows:")
print()
printSequenceWithWindows(protein_sequence, final_h_array)

print("QUESTION 1(b)")
print()
print("Beta Strand Windows:")
print()
printSequenceWithWindows(protein_sequence, final_b_array)

print("QUESTION 1(c)")
print()
print("Helix and Beta strands:")
print()
printThreeSeq(protein_sequence, final_h_array, final_b_array)

print("QUESTION 1(c)")
print()
print("Final Windows (Merged Helix and beta strands):")
print()
printSequenceWithWindows(protein_sequence, final_array)
