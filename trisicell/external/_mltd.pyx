cdef extern from "mltd/mltd.h":
	ctypedef struct MLTDResult:
		int distance
		int similarity
		double normalized_similarity
	cdef MLTDResult calc_mltd(const char* tree1, const char* tree2)

def run_mltd(file1, file2):
	cdef bytes file1_bytes = file1.encode()
	cdef char* cfile1 = file1_bytes
	cdef bytes file2_bytes = file2.encode()
	cdef char* cfile2 = file2_bytes
	cdef MLTDResult result = calc_mltd(cfile1, cfile2)
	result = {
		'distance': result.distance,
		'similarity': result.similarity,
		'normalized_similarity': result.normalized_similarity,
	}
	return result
