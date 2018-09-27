///////////////////////////////////////////////////////////////////////////////
// by TuanPM, 2018
// 
#include <iostream>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include "hungarian.hpp"

using namespace std;

/**
 * @brief maximum_weighted_bipartite_matching_with_limitation_of_pairs
 * @param matrix graph of value
 * @param k limitation of number of pair(row, col)
 * @return string represent the matching. ex: "(1-2;2-0;3-1)"
 */
string maximum_weighted_bipartite_matching_with_limitation_of_pairs(vector< vector<int> > matrix, int k)
{
	int maxSum = 0;
	int n = matrix.size();
	vector<string> matchingCollection;
	HungarianAlgorithm hungarian;
	char buffer [1000] = {0};
	string matching;
	vector<int> assignment;	
	
	// just use hungarian algorithm 1 time
	if (n == k) {
		maxSum = hungarian.SolveMax(matrix, assignment);
		std::cout << "max = " << maxSum << "\n";			
		matching = "(";
		for(unsigned int i = 0;i < assignment.size(); i++) {
			if (i) matching.append(";");
			sprintf(buffer, "%d-%d", i, assignment[i]);
			memset(buffer, 0x00, sizeof(buffer));
			matching.append(buffer);
		}
		matching.append(")");		
		return matching;
	}
	
	// Create a vector for permuation
	vector<int> pv(k, 1);
	pv.resize(n, 0);
	
	// Do permuation
	do {
		// extract a a (k)x(n) matrix from matrix		
		vector< vector<int> > submatrix;
		vector<int> rowMapping;
		for (unsigned int i = 0; i < pv.size(); i++) {
			if (pv[i]) {
				submatrix.push_back(matrix[i]);
				rowMapping.push_back(i); // this row of submatrix mapping thi row i of matrix
			}
		}
		
		// at this point we have submatrix - a (k)x(n) matrix
		// n is number of row of the matrix
		// solve hungarian match on submatrix		
		int sum = hungarian.SolveMax(submatrix, assignment);
		
		// 
		if (maxSum > sum) continue;
		
		// create matching string
		matching = "(";
		for(unsigned int i = 0;i < assignment.size(); i++) {
			if (i) matching.append(";");
			sprintf(buffer, "%d-%d", rowMapping[i], assignment[i]);
			matching.append(buffer);
			memset(buffer, 0x00, sizeof(buffer));
			matching.append(buffer);
		}
		matching.append(")");
		
		// check sum
		if (maxSum < sum){
			matchingCollection.clear();
			matchingCollection.push_back(matching);
			maxSum = sum;
		} else { // maxSum == sum
			// lexicographical_compare 
			if (matchingCollection.size() == 1 
				&& std::lexicographical_compare(matching.begin(), matching.end(), matchingCollection[0].begin(), matchingCollection[0].end())) {
				matchingCollection.clear();
				matchingCollection.push_back(matching);
			}
		}
				
		matching.clear(); 
		assignment.clear();
	} while ( std::prev_permutation(pv.begin(), pv.end()) );

	std::cout << "max = " << maxSum << "\n";

	return matchingCollection[0];
}

// Test
int main(void)
{
	// -----------------		  
    // please use "-std=c++11" for this initialization of vector.
	vector< vector<int> > matrix = 
						{ 	
							{1, 1, 4},
							{3, 1, 8},
							{9, 1, 8},
							{5, 6, 10},
						};
	std::cout << "assignment = " << maximum_weighted_bipartite_matching_with_limitation_of_pairs(matrix, 3) << "\n";
	return 0;
}