
#scores a peptide in terms of MHC affinity.
def score_peptide(peptide)
    return 1

#Goes over a string starting with peptides of length n1 to peptides of length n2  (i.e. 8 to 12 in length) and returns
#a list of lists where each list gives us the scores for a sliding window over one peptide length (i.e. all peptides length 8).
# input:
#   n1: starting peptide length
#   n2: ending peptide length
#   peptide string: i.e. "AMRTSYKLMNPRSTN"
#   mutation pos: character position (counting from 0) of the mutation.
def sliding_window(n1, n2, peptide_string, mutation pos):
    return None



peptide_scores = sliding_window(8, 12, "AMRTSYKLMNPRSTNANDFKRPN", 7)












