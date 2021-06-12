from libc.stdlib cimport free, malloc


cdef extern from "scite/findBestTrees.h":
	cdef int main_in_c(int argc, char **argv)

def run_scite(cmd):
	cdef char **c_argv
	args = [bytes(x, encoding="utf8") for x in cmd]
	c_argv = <char**>malloc(sizeof(char*) * len(args))
	for idx, s in enumerate(args):
		c_argv[idx] = s
	main_in_c(len(args), c_argv)
	return None
