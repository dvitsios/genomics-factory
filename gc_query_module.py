'''
-------------------------------------------------------------------------------
SYNOPSIS

	Provides functionality for queyring the GC field and retrieving
	values of all possible ref-alt pairs for an alt allele at a given index.
	
	It also allows inference of the total number of alleles from the length 
	of a GC comma-separated list

AUTHORS         
	Dimitrios Vitsios
'''
import sys
	

def get_num_of_alleles_from_gc(gc_len):
	"""
		Infer number of alleles (both ref and alternative) from length of GC comma-separated list
		Input: 
		      - Length of GC field
			Possible values are: Σ(i), i=1,...,N, where N = number_of_alleles
			so gc_len = 3, 6, 10, 15, 21, etc.

		Output:
		      - Number of alleles (both ref and alternative): 
			the lowest 'i' that gives [Sum{i}, i=1,...N] == gc_len
	"""

	num_alleles = None
	sum = 0

	for i in range(1, gc_len):
		sum += i
		if sum == gc_len:
			num_alleles = i
			return num_alleles


def get_gc_indices_for_alt_allele(num_alleles, alt_idx=0, verbose=True):
	"""
		Retrieve all alternative alleles at a specific index, 
		e.g. for alt-allele: D at index: 2, it returns the indices of 
		AD, BD, CD and DD values from GC
	"""

	# --- Sanity check for input arguments ---
	if (not isinstance(num_alleles, int)) or num_alleles < 2:
		raise ValueError('num_alleles must be integer >= 2')
	
	if (not isinstance(alt_idx, int)) or alt_idx < 0:
		raise ValueError('alt_idx must be non negative integer')

	if alt_idx >  num_alleles - 2:
		raise ValueError('alt_idx must be <= (num_alleles - 2)')
	# -----------------------------------------


	# list of uppercase letters: ['A', 'B', 'C', ...]
	generic_allele_lookup_table = list(map(chr, range(65, 91)))

	allele_idx = alt_idx + 1	
	cur_pairs_list = []

	for i in range(num_alleles):
		for j in range(i+1):
			ref_alt_pair = generic_allele_lookup_table[i] + generic_allele_lookup_table[j]
			ref_alt_pair = ref_alt_pair[::-1]
			cur_pairs_list.append(ref_alt_pair)


	alt_allele = generic_allele_lookup_table[ allele_idx ]
	indices = [i for i,s in enumerate(cur_pairs_list) if alt_allele in s]

	if verbose:
		print("   Alt-allele:",alt_allele)
		print("   Indices:", indices)
		print("   All pairs:", cur_pairs_list)
		print("   Selected pairs:", [cur_pairs_list[i] for i in indices])

	return indices
	

def get_limited_gc_indices_for_alt_allele(num_alleles, alt_idx=0, verbose=True):
	"""
		For the alt allele at a specific index (e.g. D at index: 2)
		it returns the indices of AD, and DD values from GC
	"""

	# --- Sanity check for input arguments ---
	if (not isinstance(num_alleles, int)) or num_alleles < 2:
		raise ValueError('num_alleles must be integer >= 2')
	
	if (not isinstance(alt_idx, int)) or alt_idx < 0:
		raise ValueError('alt_idx must be non negative integer')

	if alt_idx >  num_alleles - 2:
		raise ValueError('alt_idx must be <= (num_alleles - 2)')
	# -----------------------------------------


	# list of uppercase letters: ['A', 'B', 'C', ...]
	generic_allele_lookup_table = list(map(chr, range(65, 91)))

	allele_idx = alt_idx + 1	
	cur_pairs_list = []

	for i in range(num_alleles):
		for j in range(i+1):
			ref_alt_pair = generic_allele_lookup_table[i] + generic_allele_lookup_table[j]
			ref_alt_pair = ref_alt_pair[::-1]
			cur_pairs_list.append(ref_alt_pair)


	alt_allele = generic_allele_lookup_table[ allele_idx ]
	limited_pairs = []
	limited_pairs.append('A' + alt_allele) # AB
	limited_pairs.append(alt_allele * 2)  # BB

	indices = [i for i,pair in enumerate(cur_pairs_list) if pair in limited_pairs]

	if verbose:
		print("   Alt-allele:",alt_allele)
		print("   Indices:", indices)
		print("   All pairs:", cur_pairs_list)
		print("   Limited pairs:", limited_pairs)
		print("   Selected pairs:", [cur_pairs_list[i] for i in indices])

	return indices


def run_examples():
	# Examples
	GC_list = [[5601,3018,5795], [685, 13609, 272, 26, 17, 0, 6, 3, 0, 0]]

	for i in range(len(GC_list)):

		GC = GC_list[i]
		print("\n>> Example " + str(i+1) + ":", GC)

		num_alleles = get_num_of_alleles_from_gc(len(GC))
		print(" - Num. of alleles:", num_alleles)

		for alt_idx in range(num_alleles-1):
			print("\n - Look up alt-allele at index:", alt_idx)
			print("   ... get_gc_indices_for_alt_allele(", num_alleles, ",", alt_idx, "):\n")
			ret_idxs = get_gc_indices_for_alt_allele(num_alleles, alt_idx)
			ret_idxs = get_limited_gc_indices_for_alt_allele(num_alleles, alt_idx)
	

if __name__ == '__main__':

	print(">> get_gc_indices_for_alt_allele(num_alleles = 4, alt_idx = 1):\n")
	ret_idxs = get_gc_indices_for_alt_allele(4, 2)
	print(ret_idxs)

	print('+++++++++++++++++++++')
	ret_limited = get_limited_gc_indices_for_alt_allele(4, 1)
	print(ret_limited)

	run_examples()
