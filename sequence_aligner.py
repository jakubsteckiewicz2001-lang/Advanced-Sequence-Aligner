import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- 1. FUNKCJA OBLICZAJĄCA PUNKTY (MACIERZ) ---
def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((n + 1, m + 1))
    
    for i in range(n + 1): score_matrix[i][0] = i * gap
    for j in range(m + 1): score_matrix[0][j] = j * gap
        
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j-1]:
                diag = score_matrix[i-1][j-1] + match
            else:
                diag = score_matrix[i-1][j-1] + mismatch
            score_matrix[i][j] = max(diag, score_matrix[i-1][j] + gap, score_matrix[i][j-1] + gap)
    return score_matrix

# --- 2. FUNKCJA ODTWARZAJĄCA DOPASOWANIE (TRACEBACK) ---
def get_alignment(seq1, seq2, score_matrix, gap=-2):
    align1, align2 = "", ""
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        current_score = score_matrix[i][j]
        if i > 0 and j > 0 and (score_matrix[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -1) == current_score):
            align1 += seq1[i-1]; align2 += seq2[j-1]
            i -= 1; j -= 1
        elif i > 0 and score_matrix[i-1][j] + gap == current_score:
            align1 += seq1[i-1]; align2 += "-"
            i -= 1
        else:
            align1 += "-"; align2 += seq2[j-1]
            j -= 1
    return align1[::-1], align2[::-1]

# --- 3. FUNKCJA RYSUJĄCA WYKRES ---
def plot_matrix(matrix, seq1, seq2):
    plt.figure(figsize=(10, 8))
    x_labels = ["-"] + list(seq2)
    y_labels = ["-"] + list(seq1)
    sns.heatmap(matrix, annot=True, fmt=".0f", cmap="YlGnBu", xticklabels=x_labels, yticklabels=y_labels)
    plt.title("Needleman-Wunsch Alignment Matrix")
    plt.show()

# --- 4. URUCHOMIENIE I WYNIKI ---
s1 = "GATTACA"
s2 = "GCATGCUA"

# Najpierw liczymy macierz!
matrix = needleman_wunsch(s1, s2)

# Potem robimy dopasowanie na podstawie tej macierzy
a1, a2 = get_alignment(s1, s2, matrix)

print("Macierz dopasowania zapisana pomyślnie.")
print(f"Final Alignment Score: {matrix[len(s1)][len(s2)]}")
print("\nNajlepsze dopasowanie (Alignment):")
print(a1)
print(a2)

# Na końcu pokazujemy wykres
plot_matrix(matrix, s1, s2)