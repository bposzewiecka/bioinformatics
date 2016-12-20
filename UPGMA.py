from Bio import SeqIO

GAP_PENALTY = -1
MATCH_SCORE= 0
MISMATCH_SCORE  = -1

# returns list of aligned profiles using Needleman-Wunsch algorithm
def alignment(seqs1, seqs2): 

	# score of two triplets
	def score(r1, r2):
		if r1 ==  r2:  return MATCH_SCORE
		if r1 == '---'  or  r2 == '---': return GAP_PENALTY 		
		return MISMATCH_SCORE

	#score for i-th comumn of seq1 and j-th column of seq2
	def score_seqs(i, j):
		return sum( [ score(seq1[i], seq2[j]) for seq1  in seqs1 for seq2  in seqs2 ])

	seqsT1 = list(zip(*seqs1))
	seqsT2 = list(zip(*seqs2))

	l1 = len(seqs1[0])
	l2 = len(seqs2[0])

	sl1 = len(seqs1)
	sl2 = len(seqs2)

	A = [x[:] for x in [[0]*(l2 + 1)]*(l1 + 1)]
	B = [x[:] for x in [[()]*(l2 + 1)]*(l1 + 1)]

	gp =  GAP_PENALTY * sl1 * sl2 

	for i in range(l1):
		A[i + 1][0] =  gp* (i + 1)
		B[i + 1][0] = (-1, 0)

	for i in range(l2):
		A[0][i + 1] = gp * (i + 1)
		B[0][i + 1] = (0, -1)

	for i in range(l1):
		for j in range(l2):

			A[i + 1][j + 1] = A[i][j] + score_seqs(i,j)
			B[i + 1][j + 1] = (-1, -1)

			if A[i + 1][j] + gp > A[i + 1][j + 1] :
				A[i + 1][j + 1] = A[i + 1][j]  + gp 
				B[i + 1][j + 1] = (0, -1)

			if A[i][j + 1] + gp > A[i + 1][j + 1] :
				A[i + 1][j + 1] = A[i][j + 1]  + gp 
				B[i + 1][j + 1] = (-1, 0)

# Finds aligned sequences using backtrack information
	def backtrack():		
		i = l1
		j = l2

		l = []

		while i != 0 or j != 0:
			if B[i][j] == (-1, -1):
				i = i - 1
				j = j - 1
				l.append( seqsT1[i] + seqsT2[j])
			elif B[i][j] == (-1, 0):
				i = i - 1
				l.append( seqsT1[i] + ('---',) * sl2)
			else:
				j = j -1
				l.append(  ('---',) * sl1 + seqsT2[j] )

		return [ list(reversed(a)) for a in zip(*l)]

	return backtrack()

# Distance of two sequences (seq1, seq2)
# We align two sequences and return the percent of mismatches between this alignments
def tripleDist(seq1, seq2):
	[alignSeq1, alignSeq2] =  alignment([seq1], [seq2])
	return sum([  l1 != l2   for  l1, l2  in zip(alignSeq1, alignSeq2) ] ) / len(alignSeq1)

# Flatmap on tuple of int
def flatmap(t):
	if  type( t ) == int: 
		return [t]
	else:
		t1, t2 =  t
		a1 = flatmap(t1)
		a2 = flatmap(t2)

		a1.extend(a2)
		return a1

# returnslist of three DNA triplet sequences created form DNA sequence (taking into account 3 ORFs)
def getTripletSeq(seq):

	l = len(seq)

# alignment can contain gaps that are not multiple of 3 at the begining and at the end of alignment
	tripletSeq1 = seq + '-' * ((3 - l) % 3)
	tripletSeq2 = '-'  + seq + '-' * ((3 - l- 1) % 3)
	tripletSeq3 = '--' + seq + '-' * ((3 - l - 2) % 3)

# alignment contains only gaps that are multuple of 3
#	tripletSeq1 = seq[0: l // 3 * 3]
#	tripletSeq2 = seq[1: (l - 1) // 3 * 3 + 1]
#	tripletSeq3 = seq[2: (l - 2) // 3 * 3 + 2]

	def pack(tripleSeq):
		return [ tripleSeq[3 * i: 3 * i + 3 ] for i in range(len(tripleSeq) // 3)]
			
	return [pack(tripletSeq1), pack(tripletSeq2), pack(tripletSeq3)]

# MSA heuristic using UPGMA
def UPGMA(seqs):

	tripletSeqs = []

	#create list of triplets. 3 DNA triplet are created from 1 input DNA
	for seq in seqs:
		for triplet in getTripletSeq(seq):
			tripletSeqs.append(triplet)

	tlen = len(tripletSeqs)

	# every three triple DNA sequences comes from the same DNA sequence
	# this function check if two diffrent DNA triplet sequences comes from the same DNA sequence
	def isSibling(e1, e2):
		if type( e1 ) != int or type( e2 ) != int:
			return False

		return e1 != e2 and e1 // 3 == e2 // 3

	# compute distances for all DNA triplet sequences that are not siblings
	dist = { (i,j) : tripleDist(tripletSeqs[i], tripletSeqs[j])   
					for i in range(tlen) for j in range(tlen)  
					if (i // 3 != j // 3) and len(tripletSeqs[i]) != 0  and len(tripletSeqs[j]) != 0 }


	# dictionary where keys are trees and values contains the numeber of nodes in those trees
	elems = { i : 1 for i in range(tlen) if len(tripletSeqs[i]) != 0 }

	# dictionary of aligned profiles of all triplets
	upgma = { i: [triplet]  for i, triplet  in enumerate(tripletSeqs) }

	while len(elems) != 1:

		# find elements with minimal distance (m1 i m2)
		m = min(dist, key=dist.get)
		m1, m2 = m

		# delete distances from this element from dictionaty
		del dist[m2, m1]
		del dist[m1, m2]

		# delete distances between elements that are siblings of m1 or m2
		for e1, e2 in list(dist):
			if isSibling(m1, e1) or isSibling(m1, e2) or isSibling(m2, e1) or isSibling(m2, e2):
				del dist[e1, e2]

		# delete siblings of m1 and m2 from element list
		for e in list(elems):
			if isSibling(m1, e) or isSibling(m2, e):
				del elems[e]

		# align two profiles
		upgma[m] =  alignment(upgma[m1], upgma[m2])

		# delete profiles of m1 and m2 or its siblings
		for e in list(upgma): 
			if isSibling(m1, e) or isSibling(m2, e) or e == m1 or e == m2:
				del upgma[e]

		# update the distance dictionary using UPGMA method with Arithmetic Mean
		for elem in elems:
			if elem not in (m1, m2):
				d = (dist[m1, elem] * elems[m1] + dist[m2, elem] * elems[m2]) / ( elems[m1] +   elems[m2])
				dist[m, elem] = d
				dist[elem, m] = d

		# update number of elements in new tree
		elems[m] = elems[m1] + elems[m2]

		#delete elements m1 and m2 from dictionaty elems
		del elems[m1]
		del elems[m2]

		# delete all distances form element m1 or m2
		for e1, e2 in list(dist):
			if e1 == m1 or e1 == m2 or e2 == m1 or e2 == m2: del dist[e1, e2]
		
	aligns = [''.join(i) for i in list(upgma.values())[0] ]
	align_tree = list(upgma)[0]

	result = [ (idxs // 3,  idxs % 3, align) for idxs, align in zip(flatmap(align_tree), aligns)]
	
	#result is list of tuples each containing sequence number, number of ORF (0, 1 or 2) and aligned sequences sorted by sequence number

	return sorted(result,  key=lambda x: x[0]) 


seqs = [ str(record.seq) for record in list(SeqIO.parse("test.fasta", "fasta")) ] 

for i in UPGMA(seqs):
    print(i)
