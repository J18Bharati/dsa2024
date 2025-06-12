import sys

def backward_pairs(rna_sequence, index):
	#Input the sequence and index of a base to find all possible pairings with lower index bases.
	n = len(rna_sequence)
	if index-4 <= 0:
		yield (None, None, None)
		return
	for i in range(0, index-4):
		pair = "".join(sorted(rna_sequence[i]+rna_sequence[index]))
		if pair == "AU" or pair == "CG":
			yield (i, index, rna_sequence[i]+rna_sequence[index])
	yield (None, None, None)

def solve(rna_sequence):
	n = len(rna_sequence) #Length of the string.
	mem_mat = [[(0, []) for _ in range(n)] for _ in range(n)] #Matric to store the soltions as a tuple containing the max value and the pairs.
	mem_backpairs = [list(backward_pairs(rna_sequence, i)) for i in range(n)] #For each base, possible pairings with previous bases.

	def opt(i, j):
		if i >= j-4: #Return the default tuple since the conditions are not met.
			return (0, [])
		if mem_mat[i][j] != (0, []): #If we already have a solution in memory, return that.
			return mem_mat[i][j]

		#These are the solutions for different possible values of t. 
		solutions = [(1 + opt(i, t-1)[0] + opt(t+1, j-1)[0], [(t, j)] + opt(i, t-1)[1] + opt(t+1, j-1)[1]) for t, _, _ in mem_backpairs[j] if t!=None and t>=i]
		
		#Append to the list of possible solutions the other case.
		solutions.append(opt(i, j-1))

		#Find the maximum, store it and return it.
		mem_mat[i][j] = max(solutions)
		return mem_mat[i][j]
	
	#We start from smaller intervals and build up, referencing solutions to smaller intervals in the memory matrix.
	for interval_length in range(5, n): 
		for start_point in range(n - interval_length):
			opt(start_point, start_point + interval_length)

	return mem_mat[0][n-1]#max([max(mem_mat[i]) for i in range(n)])


def main():
	#Optionally provide the sring as a command line argument. Defaults to the problem string.
	try:
		rna_sequence = sys.argv[1].upper()
	except IndexError:
		rna_sequence = "AUGGCUACCGGUCGAUUGAGCGCCAAUGUAAUCAUU"

	value, pairs = solve(rna_sequence) #Get the answers.

	print("Max number of pairs:", value)

	for pair in pairs:
		print(pair, rna_sequence[pair[0]] + "-" + rna_sequence[pair[1]])

if __name__ == "__main__":
	main()