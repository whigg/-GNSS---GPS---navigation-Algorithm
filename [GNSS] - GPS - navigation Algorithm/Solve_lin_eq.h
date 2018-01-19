#pragma once
#include "Header.h"
#include <vector>
using namespace std;
namespace ublas = boost::numeric::ublas;
#include "Inverse_Matrix.h"

vector<double> operator=(const vector<double>& v,const SatCoord& Sc) {

}

pair <vector<vector<double> >, double> SolveLinearEq(const vector<SatCoord>& SatCoords, const vector<vector<double> > & TETA,
	const double& delay, const vector<int>&  VecDeltTimeGps) {
	int N = SatCoords.size();
	vector<double>  R(N);
	vector<vector<double> > H(N,vector<double>(4));
	for (int i = 0; i < N; i++)
		H[i][4] = 1;
	vector<double> tempR(3);
	vector<double> tempi(3);
	for (int sat = 0; sat < N; sat++) {
		tempR = tempi+tempR;
	}

	pair <vector<vector<double> >, double> t;
	return t;
};