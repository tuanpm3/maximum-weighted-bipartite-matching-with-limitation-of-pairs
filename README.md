maximum-weighted-bipartite-matching-with-limitation-of-pairs
============================================================

C++ program to compute the maximum weighted bipartite matching of a graph but with a limitation of matching pairs

#### Input
* matrix : a matrix (n x m) which represent the value weighted bipartite graph.
* k : limitation of matching pairs, k <= min(n, m)

#### Output
A string that represent the matching in which total of the weight is maximun.
If there are more than one result, the lowest in lexicographicaly order will be selected.

#### Dependencies
This C++ program use permutation algorithm combine together with hungarian algorithm.
I use implementation of hungarian algorithm from https://github.com/mcximing/hungarian-algorithm-cpp.git but with some modification.

#### Hungarian algorithm
https://en.wikipedia.org/wiki/Hungarian_algorithm

#### Building
If you are building directly from this repository you need to build on Linux with 
* g++ 
* make

Then, to build the program you just need to do one of these:

* make
* make debug
* make release

#### Example
matrix = {
	{1, 1, 4},
	{3, 1, 8},
	{9, 1, 8},
	{5, 6, 10}
}

k = 3.

==> output string = "(1-2;2-0;3-1)"
==> max sum of weight = matrix[1][2] + matrix[2][0] + matrix[3][1] = 8 + 9 + 6 = 23.
