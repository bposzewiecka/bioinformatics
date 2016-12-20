from itertools import product
from operator import add, sub
from Bio import SeqIO
import sys

def CarilloLippmanMSA(seqs):

# Number of sequences
    dims = len(seqs)

# Sequences lengths
    lengths = [ len(seq) for seq in seqs]

# Cartesian product of "dims" dimension. 
# Each dimension has two elements: 0 and 1. Element containing only zeroes is excluded.
    deltas = list(product([0, 1], repeat = dims))[1:] 

# Cartesian product equivalent to "dims" nested for-loops. 
# The number of iterations of each loop is given by the length of sequence plus 1
    def getProduct(seqs):
        ranges = map(range, map(lambda s: len(s) + 1, seqs))
        return product(*ranges)  

# Score of two letters s1 and s2 from the alphabet with a gap symbol '-'
    def pairScore(s1, s2):
        if s1 == '-' and s2 == '-': 
            return 0
        elif s1 == s2: 
            return 1
        else:
            return -1

# Letter for each sequence given it's positions (vector point) 
# and vector delta (0 denotes a gap symbol and 1 denotes symbol from the alphabet)
    def getSymbols(point, delta): 
        symbols = [ seqs[i][d - 1 ]  for i, d in  enumerate(point) ]

        for i, d in enumerate(delta): 
            if d == 0: symbols[i] = '-'

        return symbols

# Sum of pairs score of letters and gap symbols '-'. Symbols are given by getSymbol function on point and delta
    def getScore(point, delta): 

        symbols = getSymbols(point, delta)

        pairs = [ (symbols[i], symbols[j])  for i in range(dims) for j in range(dims)  if i < j ]

        return sum([ pairScore(s1, s2) for s1, s2 in pairs])

# Smith–Waterman algorithm applied to i-th and j-th sequence
# If reversed is True algorithm align reversed sequences. The default reverse value is False
# Function returns a dictionary of alignment scores where keys are pairs of lengths of all prefixes of those two sequences
    def getPairAlignment(i, j, reverse = False):
        S = {}

        s1 = seqs[i]
        s2 = seqs[j]        

        if reverse:
            s1 = s1[::-1]
            s2 = s2[::-1]

        l1 = len(s1)
        l2 = len(s2)

        for i in range(l1 + 1):
            S[i, 0] = 0

        for i in range(l2 + 1):
            S[0, i] = 0

        for i in range(l1):
            for j in range(l2):
                S[i + 1, j + 1] = max(S[i, j + 1] + pairScore(s1[i], '-'), 
                                      S[i + 1, j] + pairScore('-', s2[j]), 
                                      S[i , j] + pairScore(s1[i], s2[j]), 0)
        return S

# Function generates dictionary of "getPairAlignment" function results for every pair of sequences
    def getPairAlignments():
        pA = {}

        for i in range(dims):
            for j in range(dims):
                if i > j:
                    pA[i,j] = getPairAlignment(i, j)
        return pA

# Simple heuristic local MSA based on ungapped alignment
# Function returns the maximum score of local alignment
    def getMSAHeuristicAlignment():
        H = {}

        for point in getProduct(seqs):

            prevPoint = tuple([ p - 1 for p in point ])

            if -1 in prevPoint:
                continue

            score = getScore(point, (1,) * dims)
            H[point] = max(score + H.get(prevPoint, 0), 0)

        return max(H.values())

# Compute dictionary of list of constraints for each pairs of sequences. 
    def getMSALipmanConstraints(H, pA):

        # Compute maximum local alignment score for all sequence pairs
        pAValDict = { key :max(pA[key].values())  for key in pA }  

        # Compute sum of all maximum local alignment scores for all sequence pairs
        pASum = sum(pAValDict.values())

        # Constraints 
        constraints = {}

        for k in range(dims):
            for l in range(dims):
                if k  >  l:

                    # Carillo and Lippman lower bound
                    bound =  pAValDict[k, l] + H - pASum

                    # alignment of k-th and l-th sequences
                    forward = pA[k, l]

                    #alignment of k-th and l-th sequences reversed
                    backward = getPairAlignment(k,l, reverse = True)

                    #compute best local alignment score of sequence k and l that goes through i and j

                    pathsVals = {}

                    for i, j in forward:
                        pathsVals[i,j] = forward[i,j] + backward[ lengths[k] - i , lengths[l] - j ]

                    # add to constraints all points that have best local alignment score greater or equal to Carillo and Lippman lower bound
                    constraints[k,l] = [ (i, j) for i, j in pathsVals if pathsVals[i, j] >= bound]

        return constraints

# Checks if point meet constraints computed using Carillo and Lippman lower bounds
    def meetConstraints(point, constraints):
     #   return True
        for i in range(dims):
            for j in range(dims):
                if i > j:
                    if (point[i],point[j]) not in constraints[i,j]:
                        return False
        return True

# Computes local alignment using Carillo and Lippman lower bounds constraints
    def getMSALipmanAlignment(H, pA):

        # Local alignment scores
        S = { (0,) * dims : 0 }

        # Backtrack information
        B = {}

        constraints =  getMSALipmanConstraints(H, pA)

        for point in getProduct(seqs):
       
        # If point does not meet constraints given by Carillo and Lippman lower bounds
        # it is not taken into consideration
            if not  meetConstraints(point, constraints):
                continue

            prevPointsAndDeltas = [ (tuple(map(sub, point, delta)) , delta) for delta in deltas ]

            for prevPoint, delta in   prevPointsAndDeltas: 

                if -1 in prevPoint:
                    continue

                score = getScore(point, delta)

                if score + S.get(prevPoint,0) > S.get(point,0):
                    S[point] = score + S.get(prevPoint,0)
                    B[point] = delta

        return S, B

# Finds best local MSA using backtrack information
# Returns a list of local alignments
    def backtrack(S, B):
        l = []

        point = max(S, key=S.get)

        while point in B:
            l.append(getSymbols(point, B[point]))
            point  = tuple(map(sub, point, B[point]))
  
        return [ ''.join(reversed(s))  for s  in list(zip(*l)) ] 

    pA = getPairAlignments()
    H = getMSAHeuristicAlignment()

    S, B = getMSALipmanAlignment(H, pA)

    return backtrack(S, B)
    

seqs = [record.seq for record in list(SeqIO.parse("test.fasta", "fasta")) ] 

for i in CarilloLippmanMSA(seqs):
    print(i)
