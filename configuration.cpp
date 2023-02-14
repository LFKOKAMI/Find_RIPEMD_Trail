#include "configuration.h"
#include "gurobi_c++.h"
//#include "/opt/gurobi910/linux64/include/gurobi_c++.h"
#include <string>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

/*
*model z=x+y
*(a0,a1): x[i]
*(a2,a3): y[i]
*(a4,a5): carry[i]
*(a6,a7): z[i]
*(a8,a9): carry[i+1]
*/
void Primitive::addition(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
	GRBVar& a8, GRBVar& a9) {

	int ma[34][11] = {
		0,1,0,1,0,1,0,-1,0,0,1,
		0,-1,0,-1,0,-1,0,1,0,0,3,
		0,0,0,0,-1,1,0,0,0,0,1,
		0,0,0,0,0,0,0,0,-1,1,1,
		-1,1,0,0,0,0,0,0,0,0,1,
		0,0,-1,1,0,0,0,0,0,0,1,
		-1,0,0,0,0,0,0,0,1,-1,2,
		1,0,-1,0,1,0,0,0,0,-1,2,
		0,0,1,0,1,0,0,0,-1,0,1,
		0,0,0,0,0,0,-1,1,0,0,1,
		1,0,0,0,1,-1,-1,0,0,0,2,
		0,0,1,-1,1,0,-1,0,0,0,2,
		1,0,1,0,-1,0,0,0,0,-1,2,
		0,-1,0,1,0,-1,0,-1,0,0,3,
		0,-1,-1,0,0,1,0,-1,0,0,3,
		0,1,0,-1,0,-1,0,-1,0,0,3,
		0,-1,0,1,0,1,0,1,0,0,1,
		1,-1,1,0,0,0,-1,0,0,0,2,
		0,0,-1,0,-1,0,1,0,1,0,2,
		0,1,0,1,0,0,0,0,0,-1,1,
		1,0,0,0,0,0,-1,0,0,-1,2,
		1,-1,1,-1,1,0,0,0,0,1,2,
		0,0,0,0,0,0,1,-1,-1,0,2,
		1,0,1,0,1,-1,0,1,0,1,1,
		0,0,0,0,0,1,0,-1,0,-1,2,
		-1,0,0,0,-1,0,1,0,0,1,2,
		-1,0,-1,0,0,0,1,0,0,1,2,
		0,1,-1,0,0,1,1,0,0,0,1,
		0,0,1,0,0,0,-1,0,0,-1,2,
		0,0,0,0,1,0,-1,0,0,-1,2,
		-1,0,0,0,0,1,1,-1,0,0,2,
		0,1,0,-1,0,1,0,1,0,0,1,
		0,1,0,1,-1,0,1,0,0,0,1,
		-1,0,-1,0,-1,0,0,0,0,1,3,
	};

	for (int i = 0; i < 34; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			a8 * ma[i][8] + a9 * ma[i][9] + ma[i][10]-1 >= 0
		);
	}
}


/*
*model z=x+0
*(a0,a1): x[i]
*(a2,a3)=(0,0): y[i] (delete)
*(a4,a5): carry[i]
*(a6,a7): z[i]
*(a8,a9): carry[i+1]
*/
void Primitive::addZero(GRBModel& model,
	GRBVar& a0, GRBVar& a1, int a2, int a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
	GRBVar& a8, GRBVar& a9) {

	int ma[34][11] = {
		0,1,0,1,0,1,0,-1,0,0,1,
		0,-1,0,-1,0,-1,0,1,0,0,3,
		0,0,0,0,-1,1,0,0,0,0,1,
		0,0,0,0,0,0,0,0,-1,1,1,
		-1,1,0,0,0,0,0,0,0,0,1,
		0,0,-1,1,0,0,0,0,0,0,1,
		-1,0,0,0,0,0,0,0,1,-1,2,
		1,0,-1,0,1,0,0,0,0,-1,2,
		0,0,1,0,1,0,0,0,-1,0,1,
		0,0,0,0,0,0,-1,1,0,0,1,
		1,0,0,0,1,-1,-1,0,0,0,2,
		0,0,1,-1,1,0,-1,0,0,0,2,
		1,0,1,0,-1,0,0,0,0,-1,2,
		0,-1,0,1,0,-1,0,-1,0,0,3,
		0,-1,-1,0,0,1,0,-1,0,0,3,
		0,1,0,-1,0,-1,0,-1,0,0,3,
		0,-1,0,1,0,1,0,1,0,0,1,
		1,-1,1,0,0,0,-1,0,0,0,2,
		0,0,-1,0,-1,0,1,0,1,0,2,
		0,1,0,1,0,0,0,0,0,-1,1,
		1,0,0,0,0,0,-1,0,0,-1,2,
		1,-1,1,-1,1,0,0,0,0,1,2,
		0,0,0,0,0,0,1,-1,-1,0,2,
		1,0,1,0,1,-1,0,1,0,1,1,
		0,0,0,0,0,1,0,-1,0,-1,2,
		-1,0,0,0,-1,0,1,0,0,1,2,
		-1,0,-1,0,0,0,1,0,0,1,2,
		0,1,-1,0,0,1,1,0,0,0,1,
		0,0,1,0,0,0,-1,0,0,-1,2,
		0,0,0,0,1,0,-1,0,0,-1,2,
		-1,0,0,0,0,1,1,-1,0,0,2,
		0,1,0,-1,0,1,0,1,0,0,1,
		0,1,0,1,-1,0,1,0,0,0,1,
		-1,0,-1,0,-1,0,0,0,0,1,3,
	};

	for (int i = 0; i < 34; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			a8 * ma[i][8] + a9 * ma[i][9] + ma[i][10] - 1 >= 0
		);
	}
}

/*
*model the expansion of the signed difference
*(a0,a1): z[i]
*(a2,a3): c[i]
*(a4,a5): z'[i]
*(a6,a7): c[i+1]
*/
void Primitive::expansion(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7) {
	int ma[19][9] = {
		0,1,0,1,0,-1,0,0,1,
		0,0,0,0,0,0,-1,1,1,
		0,0,-1,1,0,0,0,0,1,
		0,0,0,0,-1,1,0,0,1,
		0,-1,0,-1,0,-1,0,0,3,
		0,1,0,-1,0,1,0,0,1,
		0,0,-1,0,0,0,1,-1,2,
		-1,0,1,0,0,1,0,-1,2,
		-1,0,0,0,-1,0,0,-1,3,
		1,-1,0,0,0,0,-1,0,2,
		-1,0,0,-1,0,-1,0,0,3,
		0,1,1,0,1,0,0,-1,1,
		0,0,0,0,-1,0,-1,0,2,
		0,0,0,1,1,0,1,-1,1,
		-1,0,0,1,1,0,0,1,1,
		1,0,1,0,-1,0,0,1,1,
		1,-1,1,0,0,1,0,1,1,
		-1,0,-1,0,0,0,0,1,2,
		0,0,-1,0,1,-1,0,1,2,
	};

	for (int i = 0; i < 19; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			ma[i][8] - 1 >= 0
		);
	}
}

/*
*model 0=x-y (another expansion)
*(a0,a1): x[i]
*(a2,a3): y[i]
*(a4,a5): c[i]
*(a6,a7): c[i+1]
*/
void Primitive::zeroAddition(
	GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7
) {
	int ma[19][9] = {
		0,1,0,-1,0,1,0,0,1,
		0,0,0,0,-1,1,0,0,1,
		0,0,-1,0,0,0,-1,0,2,
		0,0,-1,1,0,0,0,0,1,
		0,-1,0,-1,0,-1,0,0,3,
		0,0,0,0,0,0,-1,1,1,
		0,1,0,1,0,-1,0,0,1,
		1,0,0,1,0,0,-1,0,1,
		0,0,0,0,-1,0,1,-1,2,
		-1,1,0,0,0,0,0,0,1,
		1,0,1,-1,1,0,0,-1,2,
		-1,0,0,1,1,0,0,-1,2,
		0,0,1,0,0,1,1,-1,1,
		-1,0,1,0,0,1,0,1,1,
		1,0,-1,0,1,0,0,1,1,
		1,-1,0,1,1,0,0,1,1,
		-1,0,0,0,-1,0,0,1,2,
		0,0,1,-1,-1,0,0,1,2,
		-1,0,0,0,0,0,1,-1,2,
	};

	for (int i = 0; i < 19; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			ma[i][8] - 1 >= 0
		);
	}
}

/*
*Simple rotation (the first strategy)
*(a0,a1): u[i]
*(a2,a3): t[i]
*(a4,a5): mu[i]
*(a6,a7): tau[i]
*(a8,a9): iota[i]
*/
void Primitive::rotation(
	GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
	GRBVar& a8, GRBVar& a9
) {
	int ma[20][11] = {
		1,0,0,0,-1,0,0,0,0,1,1,
		0,1,-1,0,0,0,1,0,0,0,1,
		0,1,0,0,0,-1,0,0,0,0,1,
		0,0,0,-1,0,0,0,1,0,0,1,
		0,0,0,1,0,0,0,-1,0,0,1,
		0,0,0,0,0,0,0,0,-1,1,1,
		0,0,0,0,-1,1,0,0,0,0,1,
		-1,0,0,0,-1,0,0,0,0,-1,3,
		1,0,0,0,0,0,0,0,-1,0,1,
		0,0,-1,0,0,0,0,1,0,0,1,
		0,1,0,0,0,0,0,0,0,-1,1,
		0,0,0,0,0,-1,0,-1,0,0,2,
		0,0,0,0,1,0,1,0,1,-1,1,
		0,-1,1,0,0,1,1,0,0,0,1,
		-1,0,0,0,1,0,1,0,0,1,1,
		0,0,0,0,0,0,-1,1,0,0,1,
		-1,0,0,0,0,0,-1,0,0,-1,3,
		0,-1,-1,0,0,0,-1,0,0,0,3,
		-1,1,0,0,0,0,0,0,0,0,1,
		1,0,1,0,0,0,-1,0,0,1,1,
	};

	for (int i = 0; i < 20; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			a8 * ma[i][8] + a9 * ma[i][9] + ma[i][10] - 1 >= 0
		);
	}
}

/*
*the second strategy for the rotation part
*model addition and expansion at the same time (costly)
*try to avoid it as possible
*mainly to handle some cases where we need the low-probability transitions
*model z=x+y (direct sum)
*(a0,a1): x[i]
*(a2,a3): y[i]
*(a4,a5): c[i]
*(a6,a7): z[i]
*(a8,a9): c[i+1]
*/
void Primitive::additionExpansion(
	GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
	GRBVar& a8, GRBVar& a9
) {
	int ma[40][11] = {
		-1,0,-1,0,-1,0,1,0,0,0,3,
		0,1,0,1,0,1,0,-1,0,0,1,
		0,-1,0,-1,0,-1,0,1,0,0,3,
		-1,1,0,0,0,0,0,0,0,0,1,
		0,0,0,0,-1,1,0,0,0,0,1,
		0,0,-1,1,0,0,0,0,0,0,1,
		0,0,0,0,0,0,0,0,-1,1,1,
		0,0,0,0,1,0,-1,0,-1,0,2,
		0,0,0,0,0,0,-1,1,0,0,1,
		0,1,0,-1,0,-1,0,-1,0,0,3,
		0,-1,0,1,0,-1,0,-1,0,0,3,
		0,-1,0,-1,0,1,0,-1,0,0,3,
		-1,0,0,0,-1,0,0,0,1,0,2,
		0,0,1,0,0,0,-1,0,-1,0,2,
		1,-1,0,0,1,0,0,0,-1,0,2,
		0,0,1,0,1,-1,0,0,-1,0,2,
		0,0,0,-1,0,-1,0,-1,0,1,3,
		0,0,-1,0,0,0,1,0,1,-1,2,
		-1,0,-1,0,0,0,0,0,1,0,2,
		0,1,1,0,0,1,1,0,0,-1,1,
		-1,0,0,0,0,0,1,0,1,-1,2,
		1,-1,1,0,-1,0,1,0,0,-1,3,
		1,0,-1,0,-1,0,-1,0,0,0,3,
		0,1,0,0,-1,0,0,0,1,-1,2,
		1,0,1,0,1,0,-1,0,0,1,1,
		0,1,0,1,0,0,1,0,1,-1,1,
		0,0,0,1,0,1,1,0,1,-1,1,
		0,1,1,-1,1,0,0,1,0,1,1,
		1,-1,1,0,0,1,0,1,0,1,1,
		1,0,0,1,1,-1,0,1,0,1,1,
		0,1,0,1,-1,0,1,0,0,1,1,
		0,1,-1,0,0,1,1,0,0,1,1,
		-1,0,0,1,0,1,1,0,0,1,1,
		0,1,1,0,0,0,0,1,-1,0,1,
		0,1,0,0,1,0,0,1,-1,0,1,
		0,0,1,0,0,1,0,1,-1,0,1,
		1,0,1,-1,1,-1,-1,0,0,0,3,
		0,1,-1,0,0,0,0,0,1,-1,2,
		0,0,-1,0,-1,0,0,0,0,1,2,
		-1,0,0,0,0,1,0,0,1,-1,2,
	};

	for (int i = 0; i < 40; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			a8 * ma[i][8] + a9 * ma[i][9] + ma[i][10] - 1 >= 0
		);
	}
}

/*
*model w=x^y^z (fast filtering model)
*(a0,a1): x[i]
*(a2,a3): y[i]
*(a4,a5): z[i]
*(a6,a7): w[i]
*/
void Primitive::XOR(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7) {
	int ma[20][9] = {
		1,0,-1,0,-1,0,-1,0,3,
		-1,0,1,0,-1,0,-1,0,3,
		-1,0,-1,0,1,0,-1,0,3,
		-1,0,-1,0,-1,0,1,0,3,
		0,1,0,1,0,1,0,-1,1,
		0,1,0,1,0,-1,0,1,1,
		0,1,0,-1,0,1,0,1,1,
		0,-1,0,1,0,1,0,1,1,
		0,1,0,-1,0,-1,0,-1,3,
		0,-1,0,1,0,-1,0,-1,3,
		0,-1,0,-1,0,1,0,-1,3,
		0,-1,0,-1,0,-1,0,1,3,
		-1,1,0,0,0,0,0,0,1,
		0,0,-1,1,0,0,0,0,1,
		0,0,0,0,-1,1,0,0,1,
		0,0,0,0,0,0,-1,1,1,
		1,0,1,-1,1,-1,-1,0,3,
		1,0,1,-1,-1,0,1,-1,3,
		1,0,-1,0,1,-1,1,-1,3,
		-1,0,1,0,1,-1,1,-1,3,
	};

	for (int i = 0; i < 20; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			ma[i][8] - 1 >= 0
		);
	}
}

/*
model w=x^y^z 
these inequalties only correspond to how to capture the conditions,
so, it should be used together with XOR()...
(a0,a1): x[i]
(a2,a3): y[i]
(a4,a5): z[i]
(a6,a7): w[i]
x,y,z  : value (monitoring variables)
*/
void Primitive::XORCut(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
	GRBVar& x, GRBVar& y, GRBVar& z) {

	int ma[24][12] = {
		0,1,0,1,1,0,1,-1,1,-1,0,2,
		0,1,0,1,1,0,1,-1,-1,1,0,2,
		0,1,1,0,0,1,1,-1,1,0,-1,2,
		1,-1,0,1,0,0,1,-1,0,1,-1,3,
		0,1,1,0,0,1,1,-1,-1,0,1,2,
		1,-1,0,1,0,0,1,-1,0,-1,1,3,
		0,1,0,0,1,-1,-1,0,-1,-1,0,4,
		0,1,0,0,1,-1,-1,0,1,1,0,2,
		0,1,1,0,0,1,-1,0,-1,0,-1,3,
		1,0,0,1,0,1,-1,0,0,-1,-1,3,
		0,1,1,0,0,1,-1,0,1,0,1,1,
		1,0,0,1,0,1,-1,0,0,1,1,1,
		0,1,0,0,-1,0,1,-1,-1,-1,0,4,
		0,1,0,0,-1,0,1,-1,1,1,0,2,
		0,1,-1,0,0,0,1,-1,-1,0,-1,4,
		-1,0,0,1,0,0,1,-1,0,-1,-1,4,
		0,1,-1,0,0,0,1,-1,1,0,1,2,
		-1,0,0,1,0,0,1,-1,0,1,1,2,
		0,0,0,1,-1,0,-1,0,1,-1,0,3,
		0,0,0,1,-1,0,-1,0,-1,1,0,3,
		0,0,-1,0,0,1,-1,0,1,0,-1,3,
		-1,0,0,0,0,1,-1,0,0,1,-1,3,
		0,0,-1,0,0,1,-1,0,-1,0,1,3,
		-1,0,0,0,0,1,-1,0,0,-1,1,3,
	};

	for (int i = 0; i < 24; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			x * ma[i][8] + y * ma[i][9] + z * ma[i][10] +
			ma[i][11] - 1 >= 0
		).set(GRB_IntAttr_Lazy, 1);
	}
}

/*
*model w=x ^ ( y | (~z) ) (fast filtering model)
*(a0,a1): x[i]
*(a2,a3): y[i]
*(a4,a5): z[i]
*(a6,a7): w[i]
*/
void Primitive::ONX(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7) {
	
	int ma[18][9] = {
		-1,0,0,0,0,-1,-1,0,3,
		-1,0,0,-1,0,0,-1,0,3,
		0,0,0,0,0,0,-1,1,1,
		0,1,1,-1,-1,0,0,1,2,
		0,1,-1,0,1,-1,0,1,2,
		-1,1,0,0,0,0,0,0,1,
		0,0,-1,1,0,0,0,0,1,
		0,0,0,0,-1,1,0,0,1,
		0,1,0,1,0,1,0,-1,1,
		0,-1,0,1,0,1,0,1,1,
		1,-1,0,0,0,-1,1,-1,3,
		1,-1,0,-1,0,0,1,-1,3,
		0,-1,1,-1,-1,0,0,-1,4,
		0,-1,-1,0,1,-1,0,-1,4,
		0,1,-1,0,-1,0,0,-1,3,
		0,1,1,-1,1,-1,0,-1,3,
		0,-1,-1,0,-1,0,0,1,3,
		0,-1,1,-1,1,-1,0,1,3,
	};

	for (int i = 0; i < 18; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			ma[i][8] - 1 >= 0
		);
	}
}

/*
*model w=x ^ ( y | (~z) )
these inequalties only correspond to how to capture the conditions,
so, it should be used together with ONX()...
*(a0,a1): x[i]
*(a2,a3): y[i]
*(a4,a5): z[i]
*(a6,a7): w[i]
*x,y,z  : value (monitoring variables)
*/
void Primitive::ONXCut(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
	GRBVar& x, GRBVar& y, GRBVar &z) {

	int ma[22][12] = {
		-1,0,0,0,0,0,-1,0,0,0,1,2,
		-1,0,0,0,0,0,-1,0,0,-1,0,3,
		0,1,0,-1,0,1,0,1,0,0,-1,2,
		0,1,0,1,0,-1,0,1,0,1,0,1,
		0,1,0,1,0,0,0,-1,0,-1,0,2,
		0,-1,0,1,0,-1,0,-1,0,1,0,3,
		0,-1,0,-1,0,1,0,-1,0,0,-1,4,
		1,0,0,0,0,1,1,-1,0,0,1,1,
		1,-1,0,1,0,0,1,0,0,-1,0,2,
		0,1,0,0,-1,0,1,-1,-1,0,0,3,
		0,1,-1,0,0,0,1,-1,1,0,0,2,
		0,1,0,0,0,1,0,-1,0,0,1,1,
		0,1,-1,0,0,0,-1,0,-1,0,0,3,
		0,1,0,0,-1,0,-1,0,1,0,0,2,
		0,-1,0,0,0,1,0,1,0,0,1,1,
		0,1,0,1,1,0,1,-1,1,0,0,1,
		0,1,1,0,0,1,1,-1,-1,0,0,2,
		1,0,0,1,0,1,-1,0,0,1,-1,2,
		0,-1,0,1,0,0,0,1,0,-1,0,2,
		0,1,0,0,1,-1,-1,0,-1,0,0,3,
		0,1,1,0,0,1,-1,0,1,0,0,1,
		-1,0,0,1,0,0,1,-1,0,1,-1,3,
	};

	for (int i = 0; i < 22; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			x * ma[i][8] + y * ma[i][9] + z * ma[i][10] +
			ma[i][11] - 1 >= 0
		);
	}
}

void Primitive::ONZ(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7) {
	//onz(x,y,z) = (x | ~y) ^ z = onx (z,x,y)= z ^ (x|~y)
	ONX(model, a4, a5, a0, a1, a2, a3, a6, a7);
}

void Primitive::ONZCut(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
	GRBVar& x, GRBVar& y, GRBVar& z) {
	//onz(x,y,z) = (x | ~y) ^ z = onx (z,x,y)= z ^ (x|~y)
	ONXCut(model, a4, a5, a0, a1, a2, a3, a6, a7, z, x, y);
}

//removing some low-probability difference transitions
void Primitive::IFZ(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7) {

	int ma[15][9] = {
		-1,0,-1,0,0,0,1,0,2,
		0,1,-1,0,0,0,1,-1,2,
		-1,0,0,1,0,0,1,-1,2,
		0,0,0,0,0,0,-1,1,1,
		-1,1,0,0,0,0,0,0,1,
		0,0,-1,1,0,0,0,0,1,
		0,0,0,0,-1,1,0,0,1,
		0,1,0,1,0,1,0,-1,1,
		1,-1,1,-1,0,0,0,1,2,
		1,-1,0,0,0,-1,-1,0,3,
		0,0,1,-1,0,-1,-1,0,3,
		1,0,1,0,0,1,-1,0,1,
		0,0,-1,0,0,-1,1,-1,3,
		-1,0,0,0,0,-1,1,-1,3,
		0,-1,0,-1,0,1,0,1,2,
	};

	for (int i = 0; i < 15; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			ma[i][8] - 1 >= 0
		);
	}
}

void Primitive::IFZCut(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
	GRBVar& x, GRBVar& y, GRBVar& z) {

	int ma[26][12] = {
		0,-1,0,0,0,1,0,1,0,0,-1,2,
		0,0,0,-1,0,1,0,1,0,0,1,1,
		0,1,0,0,0,1,0,-1,0,0,-1,2,
		0,0,0,1,0,1,0,-1,0,0,1,1,
		0,1,0,0,1,-1,-1,0,-1,0,0,3,
		0,0,0,1,1,-1,-1,0,0,1,0,2,
		0,1,0,0,-1,0,-1,0,1,0,0,2,
		0,0,0,1,-1,0,-1,0,0,-1,0,3,
		1,0,0,0,0,1,-1,0,0,0,-1,2,
		0,0,1,0,0,1,-1,0,0,0,1,1,
		-1,0,0,0,0,0,1,-1,0,0,-1,3,
		0,0,-1,0,0,0,1,-1,0,0,1,2,
		0,1,0,0,1,-1,1,-1,1,0,0,2,
		0,0,0,1,1,-1,1,-1,0,-1,0,3,
		0,1,0,0,-1,0,1,-1,-1,0,0,3,
		0,0,0,1,-1,0,1,-1,0,1,0,2,
		0,1,0,1,0,-1,0,1,1,-1,0,2,
		0,1,0,1,0,-1,0,1,-1,1,0,2,
		0,1,-1,0,1,-1,0,1,1,0,0,2,
		-1,0,0,1,1,-1,0,1,0,-1,0,3,
		0,1,-1,0,-1,0,0,1,-1,0,0,3,
		-1,0,0,1,-1,0,0,1,0,1,0,2,
		0,1,1,-1,1,-1,0,1,-1,0,0,3,
		1,-1,0,1,1,-1,0,1,0,1,0,2,
		0,1,1,-1,-1,0,0,1,1,0,0,2,
		1,-1,0,1,-1,0,0,1,0,-1,0,3,
	};

	for (int i = 0; i < 26; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			x * ma[i][8] + y * ma[i][9] + z * ma[i][10] +
			ma[i][11] - 1 >= 0
		);
	}
}

void Primitive::IFZFull(GRBModel& model,
	GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
	GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
	GRBVar& x, GRBVar& y, GRBVar& z) {

	int ma[35][12] = {
		0,-1,0,0,0,1,0,1,0,0,-1,2,
		0,0,-1,1,0,0,0,0,0,0,0,1,
		-1,1,0,0,0,0,0,0,0,0,0,1,
		0,0,0,-1,0,1,0,1,0,0,1,1,
		0,1,0,0,0,1,0,-1,0,0,-1,2,
		0,0,0,0,-1,1,0,0,0,0,0,1,
		0,0,0,1,0,1,0,-1,0,0,1,1,
		1,-1,1,-1,0,0,0,1,0,0,0,2,
		1,-1,0,0,0,-1,-1,0,0,0,0,3,
		0,0,1,-1,0,-1,-1,0,0,0,0,3,
		0,1,0,0,1,-1,-1,0,-1,0,0,3,
		0,0,0,1,1,-1,-1,0,0,1,0,2,
		0,1,0,0,-1,0,-1,0,1,0,0,2,
		0,0,0,1,-1,0,-1,0,0,-1,0,3,
		1,0,0,0,0,1,-1,0,0,0,-1,2,
		0,0,1,0,0,1,-1,0,0,0,1,1,
		-1,0,0,0,0,0,1,-1,0,0,-1,3,
		0,0,-1,0,0,0,1,-1,0,0,1,2,
		-1,0,-1,0,0,0,0,1,0,0,0,2,
		0,0,-1,0,0,-1,1,-1,0,0,0,3,
		-1,0,0,0,0,-1,1,-1,0,0,0,3,
		0,1,0,0,1,-1,1,-1,1,0,0,2,
		0,0,0,1,1,-1,1,-1,0,-1,0,3,
		0,1,0,0,-1,0,1,-1,-1,0,0,3,
		0,0,0,1,-1,0,1,-1,0,1,0,2,
		0,1,0,1,0,-1,0,1,1,-1,0,2,
		0,1,0,1,0,-1,0,1,-1,1,0,2,
		0,1,-1,0,1,-1,0,1,1,0,0,2,
		-1,0,0,1,1,-1,0,1,0,-1,0,3,
		0,1,-1,0,-1,0,0,1,-1,0,0,3,
		-1,0,0,1,-1,0,0,1,0,1,0,2,
		0,1,1,-1,1,-1,0,1,-1,0,0,3,
		1,-1,0,1,1,-1,0,1,0,1,0,2,
		0,1,1,-1,-1,0,0,1,1,0,0,2,
		1,-1,0,1,-1,0,0,1,0,-1,0,3,
	};

	for (int i = 0; i < 35; i++) {
		model.addConstr(
			a0 * ma[i][0] + a1 * ma[i][1] + a2 * ma[i][2] + a3 * ma[i][3] +
			a4 * ma[i][4] + a5 * ma[i][5] + a6 * ma[i][6] + a7 * ma[i][7] +
			x * ma[i][8] + y * ma[i][9] + z * ma[i][10] +
			ma[i][11] - 1 >= 0
		);
	}
}

void Primitive::loadConstraint(
	GRBModel& model, string diff,
	vector<GRBVar>& v, vector<GRBVar>& d) {
	for (int i = 0; i < 32; i++) {
		if (diff[31 - i] == 'u') {
			model.addConstr(v[i] == 1);
			model.addConstr(d[i] == 1);
		}
		else if (diff[31 - i] == 'n') {
			model.addConstr(v[i] == 0);
			model.addConstr(d[i] == 1);
		}
		else if(diff[31-i]=='='|| diff[31 - i] == '0'|| diff[31 - i] == '1'){
			model.addConstr(v[i] == 0);
			model.addConstr(d[i] == 0);
		}
	}
}

void Primitive::loadCondition(
	GRBModel& model, string diff,
	vector<GRBVar>& v) {
	for (int i = 0; i < 32; i++) {
		if (diff[31 - i] == '1') {
			model.addConstr(v[i] == 1);
		}
		else if (diff[31 - i] == '0') {
			model.addConstr(v[i] == 0);
		}
	}
}


/*
new load constraint
*/
void Primitive::loadConstraint(
	GRBModel& model, string& diff,
	vector<Var>& x) {
	for (int i = 0; i < 32; i++) {
		if (diff[31 - i] == 'u') {
			model.addConstr(x[i].v == 1);
			model.addConstr(x[i].d == 1);
		}
		else if (diff[31 - i] == 'n') {
			model.addConstr(x[i].v == 0);
			model.addConstr(x[i].d == 1);
		}
		else if (diff[31 - i] == '=' || diff[31 - i] == '0' || diff[31 - i] == '1') {
			model.addConstr(x[i].v == 0);
			model.addConstr(x[i].d == 0);
		}
	}
}

int Primitive::minus(int shift) {
	return (32 + shift) % 32;
}

void Primitive::setZero(GRBModel& model, Var& var) {
	model.addConstr(var.d == 0);
	model.addConstr(var.v == 0);
}

/*
start model the round function
*/
void Primitive::boolFastModel(GRBModel &model, string &name, vector<Var>& x, vector<Var>& y, vector<Var>& z, vector<Var>& w) {
	if (name=="XOR") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			XOR(model, x[i].v, x[i].d, y[i].v, y[i].d, z[j].v, z[j].d, w[i].v, w[i].d);
		}
	}
	else if (name=="ONX") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			ONX(model, x[i].v, x[i].d, y[i].v, y[i].d, z[j].v, z[j].d, w[i].v, w[i].d);
		}
	}
	else if (name=="IFZ") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			IFZ(model, x[i].v, x[i].d, y[i].v, y[i].d, z[j].v, z[j].d, w[i].v, w[i].d);
		}
	}
	else if (name=="ONZ") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			ONX(model, z[j].v, z[j].d, x[i].v, x[i].d, y[i].v, y[i].d, w[i].v, w[i].d);
		}
	}
	else if (name=="IFX") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			IFZ(model, y[i].v, y[i].d, z[j].v, z[j].d, x[i].v, x[i].d, w[i].v, w[i].d);
		}
	}
}

void Primitive::boolConditionModel(GRBModel& model, string &name, vector<Var>& x, vector<Var>& y, vector<Var>& z, vector<Var>& w,
	vector<GRBVar>& xv, vector<GRBVar>& yv, vector<GRBVar>& zv) {
	if (name=="XOR") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			XORCut(model, x[i].v, x[i].d, y[i].v, y[i].d, z[j].v, z[j].d, w[i].v, w[i].d, xv[i], yv[i], zv[j]);
		}
	}
	else if (name=="ONX") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			ONXCut(model, x[i].v, x[i].d, y[i].v, y[i].d, z[j].v, z[j].d, w[i].v, w[i].d, xv[i], yv[i], zv[j]);
		}
	}
	else if (name=="IFZ") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			IFZCut(model, x[i].v, x[i].d, y[i].v, y[i].d, z[j].v, z[j].d, w[i].v, w[i].d, xv[i], yv[i], zv[j]);
		}
	}
	else if (name=="ONZ") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			ONXCut(model, z[j].v, z[j].d, x[i].v, x[i].d, y[i].v, y[i].d, w[i].v, w[i].d, zv[j], xv[i], yv[i]);
		}
	}
	else if (name=="IFX") {
		for (int i = 0; i < 32; i++) {
			int j = minus(i - 10);
			IFZCut(model, y[i].v, y[i].d, z[j].v, z[j].d, x[i].v, x[i].d, w[i].v, w[i].d, yv[i], zv[j], xv[i]);
		}
	}
}

void Primitive::modAddition(GRBModel& model, vector<Var>& x, vector<Var>& y, vector<Var>& carry, vector<Var>& z, int sx, int sy) {
	for (int i = 0; i < 32; i++) {
		int j = minus(i - sx);
		int k = minus(i - sy);
		addition(model, x[j].v, x[j].d, y[k].v, y[k].d, carry[i].v, carry[i].d, z[i].v, z[i].d, carry[i + 1].v, carry[i + 1].d);
	}
}

//x is expanded to y (y = x <<< shift)
void Primitive::expandState(GRBModel& model, int shift, int isK,
	vector<Var>& x, vector<Var>& y, vector<Var>& carry) {

	if (isK == 0) {
		for (int i = 0; i < 32; i++) {
			expansion(model, x[i].v, x[i].d,
				carry[i].v, carry[i].d,
				y[(i + shift) % 32].v, y[(i + shift) % 32].d,
				carry[i + 1].v, carry[i + 1].d);
		}
	}
	else {
		for (int i = 0; i < 32; i++) {
			zeroAddition(model, x[i].v, x[i].d,
				y[(i + shift) % 32].v, y[(i + shift) % 32].d,
				carry[i].v, carry[i].d,
				carry[i + 1].v, carry[i + 1].d);
		}
	}
	//for the MSB, we can ignore c[32].
	//model.addConstr(y[(31 + shift) % 32].d >= y[(31 + shift) % 32].v);
	//model.addConstr(y[(31 + shift) % 32].d + x[31].d - carry[31].d >= 0);
	//model.addConstr(y[(31 + shift) % 32].d - x[31].d + carry[31].d >= 0);
	//model.addConstr(-y[(31 + shift) % 32].d + x[31].d + carry[31].d >= 0);
	//model.addConstr(-y[(31 + shift) % 32].d - x[31].d - carry[31].d + 2 >= 0);
}

/*
y = x <<< shift
*/
void Primitive::rotateState(GRBModel& model, int shift,
	vector<Var>& x, vector<Var>& y, Aux& aux) {

	//v[31:30]
	//v[31-shift,30-shift]
	for (int i = 0; i < 30 - shift; i++) {
		model.addConstr(x[i].v == y[(i + shift) % 32].v);
		model.addConstr(x[i].d == y[(i + shift) % 32].d);
	}
	rotation(model, x[31 - shift].v, x[31 - shift].d, x[30 - shift].v, x[30 - shift].d,
		y[31].v, y[31].d, y[30].v, y[30].d,
		aux.midV, aux.midD);

	for (int i = 32 - shift; i < 30; i++) {
		model.addConstr(x[i].v == y[(i + shift) % 32].v);
		model.addConstr(x[i].d == y[(i + shift) % 32].d);
	}
	//we do not need the carry for the MSB
	rotation(model, x[31].v, x[31].d, x[30].v, x[30].d,
		y[(31 + shift) % 32].v, y[(31 + shift) % 32].d, y[(30 + shift) % 32].v, y[(30 + shift) % 32].d, aux.highV, aux.highD);
}

void Primitive::rotateFirstStrategy(GRBModel& model,int shift,
	vector<Var>&b3, vector<Var>&a1, vector<Var>&b4, 
	vector<Var>& b5, vector<Var>&a5, vector<Var>& q,
	vector<Var>& carryAdd, vector<Var>& carryExp, vector<Var>& carryQ,
	Aux& aux, Par& par) {
	rotateState(model, shift, b3, b4, aux);//b4 = b3<<<shift
	//add b4 and a1<<<10, the sum is b5
	model.addConstr(carryAdd[0].v == aux.midV);
	model.addConstr(carryAdd[0].d == aux.midD);
	for (int i = 0; i < 32; i++) {
		modAddition(model, b4, a1, carryAdd, b5, 0, 10);//b5 = b4 + al<<<10
	}
	//expand b5 and get a5
	setZero(model, carryExp[0]);
	expandState(model, 0, par.isK, b5, a5, carryExp);//expand b5 and get a5 (a5=b5)

	if (par.isV) {
		//compute q = b4
		setZero(model,carryQ[0]);
		for (int j = 0; j <= shift; j++) {
			if (j == 0) {
				addition(model, b4[j].v, b4[j].d,
					aux.midV, aux.midD,
					carryQ[j].v, carryQ[j].d,
					q[j].v, q[j].d,
					carryQ[j + 1].v, carryQ[j + 1].d);
			}
			else {
				addZero(model, b4[j].v, b4[j].d,
					0, 0,
					carryQ[j].v, carryQ[j].d,
					q[j].v, q[j].d,
					carryQ[j + 1].v, carryQ[j + 1].d);
			}
		}
	}
}

void Primitive::rotateSecondStrategy(GRBModel &model, int shift,
	vector<Var>& a1, vector<Var>&b3, vector<Var>& b4, vector<Var>& a5,
	vector<Var>& c1, vector<Var>& c2) {
	//b4 = expand(b3)
	setZero(model, c1[0]);
	zeroState(model, b3, b4, c1);
	//a5 = a1<<<10 + b4<<<shift
	setZero(model, c2[0]);
	addExpandState(model, 10, shift, a1, b4, a5, c2);
}

//yv = x<<<10 + qv
void Primitive::detectContraditions(GRBModel& model, int shift,
	vector<Var>& a1, vector<Var>& a5,
	vector<Var>& b3, vector<Var>& q,
	vector<GRBVar>& a1v, vector<GRBVar>& a5v, 
	vector<GRBVar>& qv,vector<GRBVar>& outCarry,
	vector<GRBVar>& value, vector<Var>& innerCarry, vector<Var>& outerCarry) {

	outCarry.resize(33);
	qv.resize(32);
	value.resize(32);
	innerCarry.resize(32 - shift + 2);
	outerCarry.resize(shift + 2);

	initializeVar(model, qv, 32);
	initializeVar(model, value, 32);
	initializeVar(model, outCarry, 33);
	initializeVar(model, innerCarry, 32 - shift + 2);
	initializeVar(model, outerCarry, shift + 2);

	//first compute yv=xv<<<10 + qv
	//convert the signed diff on (x,y) to conditions
	for (int j = 0; j < 32; j++) {
		model.addConstr(-a5[j].v + a5v[j] >= 0);
		model.addConstr(a5[j].v - a5[j].d - a5v[j] + 1 >= 0);
		model.addConstr(-a1[j].v + a1v[j] >= 0);
		model.addConstr(a1[j].v - a1[j].d - a1v[j] + 1 >= 0);
	}
	model.addConstr(outCarry[0] == 0);
	valueAddition(model, 10, a1v, qv, a5v, outCarry);

	//compute value=qv>>>shift + b3, until 32-shift
	setZero(model,innerCarry[0]);
	for (int j = 0; j <= 32 - shift; j++) {
		model.addConstr(
			2 * (innerCarry[j + 1].d - 2 * innerCarry[j + 1].v) + value[j]
			== (innerCarry[j].d - 2 * innerCarry[j].d) + (b3[j].d - 2 * b3[j].v)
			+ qv[(j + shift) % 32]
		);
		model.addConstr(innerCarry[j + 1].d >= innerCarry[j + 1].v);
	}

	//compute value<<<shift = qv + q, until shift:   
	//( qv + q ) = ( qv>>> shift + b3 )<<<shift  ----> (qv + q) = value <<< shift
	setZero(model, outerCarry[0]);
	for (int j = 0; j <= shift; j++) {
		model.addConstr(
			2 * (outerCarry[j + 1].d - 2 * outerCarry[j + 1].v) + value[(j + 32 - shift) % 32]
			== (outerCarry[j].d - 2 * outerCarry[j].v) + (q[j].d - 2 * q[j].v)
			+ qv[j]
		);
		model.addConstr(outerCarry[j + 1].d >= outerCarry[j + 1].v);
	}
}

void Primitive::modelRoundFunction(
	GRBModel& model, int step,
	Par& par, int shift,
	vector<Var>& m, vector<Var>& a0, vector<Var>& a1, vector<Var>& a2, vector<Var>& a3, vector<Var>& a4, vector<Var>& a5,
	vector<Var>& b1, vector<Var>& b2, vector<Var>& b3, vector<Var>& b4, vector<Var>& b5, 
	vector<Var>& c1, vector<Var>& c2, vector<Var>& c3, vector<Var>& c4, vector<Var>& c5,
	vector<GRBVar>& v1, vector<GRBVar>& v2, vector<GRBVar>& v3, vector<GRBVar>& v4, vector<GRBVar>& v5,
	vector<GRBVar>& qv, vector<GRBVar>& outCarry, vector<GRBVar>& value, vector<Var>& innerCarry, vector<Var>& outerCarry,
	vector<Var>& q, Aux& aux
	) {

	boolFastModel(model, par.fNa, a4, a3, a2, b1);//b1 = F(a4,a3,a2<<<10)
	if (par.isC == 1) {
		boolConditionModel(model, par.fNa, a4, a3, a2, b1, v4, v3, v2);
	}

	//if (step == 1) {
		//string str = "=====nn=nununuuunuuunu=uu=======";
		//loadConstraint(model, str, b1);
	//}

	setZero(model, c1[0]);
	setZero(model, c2[0]);

	modAddition(model, m, b1, c1, b2, 0, 0);//b2 = b1+m
	modAddition(model, b2, a0, c2, b3, 0, 10);//b3 = a0<<<10 + b2

	
	//b4 = b3<<< sr + a1 [shift]
	if (par.isF == 1) {//use the first strategy
		rotateFirstStrategy(model, shift, b3, a1, b4, b5, a5, q, c3, c4, c5, aux, par);
	}
	else {//use the second strategy
		rotateSecondStrategy(model, shift, a1, b3, b4, a5, c3, c4);
	}

	if (par.isV) {//detect more contraditions
		detectContraditions(model, shift, a1, a5, b3, q, v1, v5, qv, outCarry, value, innerCarry, outerCarry);
	}
	
}




bool Primitive::buildModel(int start, int end, int isRight,
	int *isC, int *isF, int *isV,
	vector<string>& diff,
	string mPattern[],
	vector<string> &bfDiff){

	int varLen = diff.size();
	int size = diff.size() - 5;
	vector<Par> par(size);
	for (int i = 0; i < size; i++) {
		par[i].isC = isC[i];
		par[i].isF = isF[i];
		par[i].isV = isV[i];
	}

	vector<vector<Var>> a(varLen);//diff
	vector<vector<GRBVar>> av(varLen);//value
	vector<vector<Var>> m(16);

	vector<vector<Var>> b(size * 5);//tmp diff
	vector<vector<Var>> q(size);//diff
	vector<vector<Var>> cMod1(size);//carry diff
	vector<vector<Var>> cMod2(size);//carry diff
	vector<vector<Var>> cMod3(size);//carry diff
	vector<vector<Var>> cExp(size);//carry diff for expand
	vector<vector<Var>> cQ(size);//carry diff for computing q

	int detectLength = 0;
	for (int i = 0; i < size; i++) {
		if (isV[i] == 1)
			detectLength++;
	}
	vector<vector<GRBVar> > qv(detectLength); 
	vector<vector<GRBVar> > outCarry(detectLength);
	vector<vector<GRBVar> > value(detectLength);
	vector<vector<Var> > innerCarry(detectLength);
	vector<vector<Var> >outerCarry(detectLength);
	int currentDetect = 0;

	vector<Aux> aux(size);

	GRBEnv env = GRBEnv();
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);

	initializeVar(model, a, 32);
	initializeVar(model, av, 32);
	initializeVar(model, m, 32);
	initializeVar(model, b, 32);
	initializeVar(model, q, 32);

	initializeVar(model, cMod1, 33);
	initializeVar(model, cMod2, 33);
	initializeVar(model, cMod3, 33);
	initializeVar(model, cExp, 33);
	initializeVar(model, cQ, 33);

	initializeVar(model, aux);

	for (int i = 0; i < 16; i++) {
		loadConstraint(model, mPattern[i], m[i]);
	}

	for (int i = 0; i < varLen; i++) {
		loadConstraint(model, diff[i], a[i]);
	}

	for (int i = start; i <size+start; i++) {
		int k = i - start;
		//cout << k << "---" << i<<" -- "<<size << endl;
		if (isRight) {
			if (i <= end)
				par[k].isK = 0;
			else
				par[k].isK = 1;//use another expansion
			par[k].fNa = fNameR[i / 16];
			modelRoundFunction(model, k, par[k], sr[i], m[pr[i]],
				a[k], a[k+1], a[k+2], a[k+3], a[k+4], a[k+5],
				b[5 * k], b[5 * k + 1], b[5 * k + 2], b[5 * k + 3], b[5 * k + 4],
				cMod1[k], cMod2[k], cMod3[k], cExp[k], cQ[k],
				av[k+1], av[k+2], av[k+3], av[k+4], av[k+5],
				qv[currentDetect], outCarry[currentDetect], value[currentDetect], innerCarry[currentDetect], outerCarry[currentDetect],
				q[k], aux[k]);

			if (par[k].isV == 1)
				currentDetect++;
			//cout << "#current detected points:" << currentDetect << endl;
		}
		else {
			if (i <= end)
				par[k].isK = 0;
			else
				par[k].isK = 1;//use another expansion
			par[k].fNa = fNameL[i / 16];
			modelRoundFunction(model, k, par[k], sl[i], m[pl[i]],
				a[k], a[k + 1], a[k + 2], a[k + 3], a[k + 4], a[k + 5],
				b[5 * k], b[5 * k + 1], b[5 * k + 2], b[5 * k + 3], b[5 * k + 4],
				cMod1[k], cMod2[k], cMod3[k], cExp[k], cQ[k],
				av[k + 1], av[k + 2], av[k + 3], av[k + 4], av[k + 5],
				qv[currentDetect], outCarry[currentDetect], value[currentDetect], innerCarry[currentDetect], outerCarry[currentDetect],
				q[k], aux[k]);

			if (par[k].isV == 1)
				currentDetect++;
		}
	}

	model.optimize();

	if (model.get(GRB_IntAttr_Status) == 3) {
		cout << "The model is infeasible" << endl;
		return 0;
	}

	string str;
	for (int i = 0; i < varLen; i++) {
		getSignedString(a[i], diff[i]);
		//cout << diff[i] << endl;
		//resDiff[i - start] = str;
	}
	for (int i = 0; i < size; i++) {
		getSignedString(b[5*i], bfDiff[i]);
		//bfDiff[i] = str;
	}
	return 1;
}


bool Primitive::buildModel(int start, int end, int opStart, int opSteps, int isRight,
	vector<int>& isC, vector<int>& isF, vector<int>& isV,
	vector<string>& diff,
	string mPattern[],
	vector<string>& bfDiff) {

	int varLen = diff.size();
	int size = diff.size() - 5;
	vector<Par> par(size);
	for (int i = 0; i < size; i++) {
		par[i].isC = isC[i];
		par[i].isF = isF[i];
		par[i].isV = isV[i];
	}

	vector<vector<Var>> a(varLen);//diff
	vector<vector<GRBVar>> av(varLen);//value
	vector<vector<Var>> m(16);

	vector<vector<Var>> b(size * 5);//tmp diff
	vector<vector<Var>> q(size);//diff
	vector<vector<Var>> cMod1(size);//carry diff
	vector<vector<Var>> cMod2(size);//carry diff
	vector<vector<Var>> cMod3(size);//carry diff
	vector<vector<Var>> cExp(size);//carry diff for expand
	vector<vector<Var>> cQ(size);//carry diff for computing q

	int detectLength = 0;
	for (int i = 0; i < size; i++) {
		if (isV[i] == 1)
			detectLength++;
	}
	vector<vector<GRBVar> > qv(detectLength);
	vector<vector<GRBVar> > outCarry(detectLength);
	vector<vector<GRBVar> > value(detectLength);
	vector<vector<Var> > innerCarry(detectLength);
	vector<vector<Var> >outerCarry(detectLength);
	int currentDetect = 0;

	vector<Aux> aux(size);

	cout << "000000" << endl;
	GRBEnv env = GRBEnv();
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);

	cout << "111111" << endl;

	initializeVar(model, a, 32);
	initializeVar(model, av, 32);
	initializeVar(model, m, 32);
	initializeVar(model, b, 32);
	initializeVar(model, q, 32);

	initializeVar(model, cMod1, 33);
	initializeVar(model, cMod2, 33);
	initializeVar(model, cMod3, 33);
	initializeVar(model, cExp, 33);
	initializeVar(model, cQ, 33);

	initializeVar(model, aux);

	for (int i = 0; i < 16; i++) {
		loadConstraint(model, mPattern[i], m[i]);
	}

	for (int i = 0; i < varLen; i++) {
		loadConstraint(model, diff[i], a[i]);
	}

	for (int i = start; i < size + start; i++) {
		int k = i - start;
		//cout << k << "---" << i<<" -- "<<size << endl;
		if (isRight) {
			if (i <= end)
				par[k].isK = 0;
			else
				par[k].isK = 1;//use another expansion
			par[k].fNa = fNameR[i / 16];
			modelRoundFunction(model, k, par[k], sr[i], m[pr[i]],
				a[k], a[k + 1], a[k + 2], a[k + 3], a[k + 4], a[k + 5],
				b[5 * k], b[5 * k + 1], b[5 * k + 2], b[5 * k + 3], b[5 * k + 4],
				cMod1[k], cMod2[k], cMod3[k], cExp[k], cQ[k],
				av[k + 1], av[k + 2], av[k + 3], av[k + 4], av[k + 5],
				qv[currentDetect], outCarry[currentDetect], value[currentDetect], innerCarry[currentDetect], outerCarry[currentDetect],
				q[k], aux[k]);

			if (par[k].isV == 1)
				currentDetect++;
			//cout << "#current detected points:" << currentDetect << endl;
		}
		else {
			if (i <= end)
				par[k].isK = 0;
			else
				par[k].isK = 1;//use another expansion
			par[k].fNa = fNameL[i / 16];
			modelRoundFunction(model, k, par[k], sl[i], m[pl[i]],
				a[k], a[k + 1], a[k + 2], a[k + 3], a[k + 4], a[k + 5],
				b[5 * k], b[5 * k + 1], b[5 * k + 2], b[5 * k + 3], b[5 * k + 4],
				cMod1[k], cMod2[k], cMod3[k], cExp[k], cQ[k],
				av[k + 1], av[k + 2], av[k + 3], av[k + 4], av[k + 5],
				qv[currentDetect], outCarry[currentDetect], value[currentDetect], innerCarry[currentDetect], outerCarry[currentDetect],
				q[k], aux[k]);

			if (par[k].isV == 1)
				currentDetect++;
		}
	}

	//minimize the hamming weight
	//int opStart = 6, opSteps = 5;
	if (opSteps > 0) {
		GRBLinExpr hwSum = 0;
		for (int i = 0; i < opSteps; i++) {
			int k = i + (opStart - start) + 5;
			cout << k - 5 << " ";
			for (int j = 0; j < 32; j++) {
				hwSum = hwSum + a[k][j].d;
			}
		}
		cout << endl;
		model.setObjective(hwSum, GRB_MINIMIZE);
	}

	model.optimize();

	if (model.get(GRB_IntAttr_Status) == 3) {
		cout << "The model is infeasible" << endl;
		return 0;
	}

	string str;
	for (int i = 0; i < varLen; i++) {
		getSignedString(a[i], diff[i]);
		//cout << diff[i] << endl;
		//resDiff[i - start] = str;
	}
	for (int i = 0; i < size; i++) {
		getSignedString(b[5 * i], bfDiff[i]);
		//bfDiff[i] = str;
	}
	return 1;
}


bool Primitive::autoCheck(int start, int isRight, 
	string msgDiff[], vector<string> &input, vector<string>& boolOut) {
	GRBEnv env;
	env.set(GRB_IntParam_Threads, 4);
	env.set(GRB_IntParam_OutputFlag, 0);
	GRBModel model = GRBModel(env);

	int inputSize = input.size();
	int bOutSize = boolOut.size();

	vector<int> s(80), p(80);
	for (int i = 0; i < 80; i++) {
		if (isRight == 1) {
			s[i] = sr[i];
			p[i] = pr[i];
		}
		else {
			s[i] = sl[i];
			p[i] = pl[i];
		}
	}

	vector<vector<Var> > x(inputSize);
	vector<vector<GRBVar> > con(inputSize);
	vector<vector<Var> > bOut(bOutSize);
	vector<vector<Var> > m(16);
	initializeVar(model, x, 32);
	initializeVar(model, con, 32);
	initializeVar(model, bOut, 32);
	initializeVar(model, m, 32);

	for (int i = 0; i < inputSize; i++) {
		loadConstraint(model, input[i], x[i]);
	}

	for (int i = 0; i < bOutSize; i++) {
		loadConstraint(model, boolOut[i], bOut[i]);
	}

	for (int i = 0; i < 16; i++) {
		loadConstraint(model, msgDiff[i], m[i]);
	}

	//convert the signed diff to conditions
	for (int i = 0; i < inputSize; i++) {
		for (int j = 0; j < 32; j++) {
			model.addConstr(-x[i][j].v + con[i][j] >= 0);
			model.addConstr(x[i][j].v - x[i][j].d - con[i][j] + 1 >= 0);
		}
	}

	for (int i = 5; i < inputSize; i++) {
		int k = i - 5 + start;
		int u = i - 5;
		string name = fNameL[k / 16];
		if (isRight==1) 
			name = fNameR[k / 16];
		boolConditionModel(model, name, x[i - 1], x[i - 2], x[i - 3], bOut[u], con[i - 1], con[i - 2], con[i - 3]);
	}

	vector<u32> inner(inputSize), outer(inputSize);
	string inStr, outStr;
	for (int i = 5; i < inputSize; i++) {
		int u = i - 5;
		int k = i - 5 + start;
		outer[i] = getValFromSignedDiff(input[i], 0) - getValFromSignedDiff(input[i - 4], 10);
		inner[i] = getValFromSignedDiff(boolOut[u], 0) + getValFromSignedDiff(input[i - 5], 10);
		inner[i] = inner[i] + getValFromSignedDiff(msgDiff[p[k]], 0);
	}

	//satisfy the modular difference
	vector<vector<GRBVar> > tmp(inputSize * 2);
	initializeVar(model, tmp, 32);
	vector<vector<GRBVar> > q(inputSize);
	initializeVar(model, q, 32);
	vector<vector<GRBVar> > c(inputSize * 3);
	initializeVar(model, c, 33);
	vector<vector<GRBVar> > inVar(inputSize);
	initializeVar(model, inVar, 32);
	vector<vector<GRBVar>> outVar(inputSize);
	initializeVar(model, outVar, 32);

	for (int i = 5; i < inputSize; i++) {
		loadConditionFromValue(model, inner[i], inVar[i]);
		loadConditionFromValue(model, outer[i], outVar[i]);
	}

	for (int i = 0; i < inputSize; i++) {
		toConditionFromSignedDiff(model, input[i], con[i]);
	}

	int testLen = inputSize;
	for (int i = 5; i < testLen; i++) {
		int k = i - 5 + start;
		model.addConstr(c[3 * i][0] == 0);
		model.addConstr(c[3 * i + 1][0] == 0);
		model.addConstr(c[3 * i + 2][0] == 0);
		valueAddition(model, 32 - s[k], q[i], inVar[i],
			tmp[2 * i], c[3 * i]);// q[i]<<< s[k] + inVar[i] = tmp[2i]
		valueAddition(model, 0, q[i], outVar[i],
			tmp[2 * i + 1], c[3 * i + 1]);// q[i] + outVar[i] = tmp[2i+1]
		valueAddition(model, 10, con[i - 4], q[i],
			con[i], c[3 * i + 2]);//con[i-4]<<<10 + q[i] = con[i]

		for (int j = 0; j < 32; j++) {
			model.addConstr(tmp[2 * i][(32 - s[k] + j) % 32] == tmp[2 * i + 1][j]);//tmp[2i]=tmp[2i+1]<<<s[k]
		}
	}

	model.optimize();
	if (model.get(GRB_IntAttr_Status) == 3) {
		cout << "The model is infeasible" << endl;
		return 0;
	}
	else {
		cout << "feasible" << endl;
		for (int i = 5; i < testLen; i++) {
			int k = i - 5 + start;
			getBinaryString(inner[i], inStr);
			getBinaryString(outer[i], outStr);
			cout << "step : " << dec << k << endl;
			cout << "shift: " << dec << s[k] << endl;
			cout << "inner: " << inStr << endl;
			cout << "outer: " << outStr << endl;
			cout << "inner: 0x" << hex << inner[i] << "   ";
			cout << "outer: 0x" << hex << outer[i] << endl << endl;
		}
	}


	vector<string> cx(inputSize);
	for (int i = 0; i < inputSize; i++) {
		cx[i] = input[i];
	}
	for (int i = 5; i < inputSize; i++) {
		int k = i - 5 + start;
		int u = i - 5;
		string name = fNameL[k / 16];
		if (isRight == 1)
			name = fNameR[k / 16];

		if (name=="ONX") {
			getConditonONX(input[i - 1], input[i - 2], input[i - 3], boolOut[u],
				cx[i - 1], cx[i - 2], cx[i - 3]);
		}
		else if (name=="IFZ") {
			getConditonIFZ(input[i - 1], input[i - 2], input[i - 3], boolOut[u],
				cx[i - 1], cx[i - 2], cx[i - 3]);
		}
		else if (name=="ONZ") {
			//getConditonONZ(input[i - 1], input[i - 2], input[i - 3], boolOut[u],
				//cx[i - 1], cx[i - 2], cx[i - 3]);
		}
		else if (name=="IFX") {
			getConditonIFX(input[i - 1], input[i - 2], input[i - 3], boolOut[u],
				cx[i - 1], cx[i - 2], cx[i - 3]);
		}
	}

	cout << "after adding bit conditions:" << endl;
	for (int i = 0; i < inputSize; i++) {
		int k = i - 5 + start;
		cout << setw(2) <<dec<< k << ": " << cx[i] << endl;
	}

	return 1;
}


















void Primitive::model_RIPEMD160_FirstRound_Left() {
	//X[i+5] = X[i+1]<<<10 + (XOR(X[i+4],X[i+3],X[i+2]<<<10) + X[i]<<<10 + m)<<s
	int s[32] = {
		11,14,15,12,5,8,7,9,11,13,14,15,6,7,9,8,
		7,6,8,13,11,9,7,15,7,12,15,9,11,7,13,12
	};

	string m7 = "=n==============n==========n====";
	string m7Inv = "=u==============u==========u====";

	//string m7 =    "=======u=======================n";
	//string m7Inv = "=======n=======================u";

	string input[10] = {
"uuuuuuuu=nuuuuuuuu=uuuuuu=======",
"uun===u=u=n=un==nuuuuu==u=nuu=nu",
"nu=nnn====n=un=n=un=====u=nu=uu=",
"=nu=========n=uuuuuu============",
"nnnnnn====un=unnnnnnnnnnnnnunnnn",
	};

	GRBEnv env = GRBEnv();
	cout << "starting" << endl;
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);
	vector<vector<GRBVar> > xv, xd, tmpv, tmpd, rxv, rxd;
	vector<GRBVar> mv, md,mvInv,mdInv;
	vector<GRBVar> midXV, midXD, highXV, highXD;

	vector<GRBVar> zero;
	zero.resize(32);
	for (int i = 0; i < 32; i++) {
		zero[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(zero[i] == 0);
	}

	xv.resize(10);
	xd.resize(10);
	rxv.resize(10);
	rxd.resize(10);

	midXD.resize(10);
	midXV.resize(10);
	highXD.resize(10);
	highXV.resize(10);

	tmpv.resize(30);
	tmpd.resize(30);

	mv.reserve(32);
	md.resize(32);
	mvInv.resize(32);
	mdInv.resize(32);

	for (int i = 0; i < 10; i++) {
		xv[i].resize(32);
		xd[i].resize(32);
		rxv[i].resize(32);
		rxd[i].resize(32);
		
		midXV[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		midXD[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		highXV[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		highXD[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

		for (int j = 0; j < 32; j++) {
			xv[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			xd[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			rxv[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			rxd[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		for (int k = 0; k < 3; k++) {
			tmpv[3 * i + k].resize(33);
			tmpd[3 * i + k].resize(33);
		}

		for (int j = 0; j < 33; j++) {
			for (int k = 0; k < 3; k++) {
				tmpv[3 * i + k][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				tmpd[3 * i + k][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
	}

	for (int i = 0; i < 32; i++) {
		mv[i]=model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		md[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		mvInv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		mdInv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	}
	loadConstraint(model, m7, mv, md);
	loadConstraint(model, m7Inv, mvInv, mdInv);

	//x7 = (m7+0)<<s[7]
	//first, rotate m7; second, expand it and assign it to x7
	rotateState(model, s[7], mv, md, rxv[0], rxd[0], 
		tmpv[0][0], tmpd[0][0], highXV[0], highXD[0]);
	expandState(model, 0, rxv[0], rxd[0], xv[0], xd[0], tmpv[0], tmpd[0]);

	//determine x8,x9,x10 by propagating signed differences
	//x8 = x4<<<10 + (xor(x7,x6,x5<<<10) + x3<<<10 )<<s[8]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[0][i], xd[0][i], zero[i], zero[i], zero[i], zero[i],
			tmpv[3][i], tmpd[3][i]);
	}
	rotateState(model, s[8], tmpv[3], tmpd[3], rxv[1], rxd[1],
		tmpv[4][0], tmpd[4][0], highXV[1], highXD[1]);
	expandState(model, 0, rxv[1], rxd[1], xv[1], xd[1], tmpv[4], tmpd[4]);

	//x9 = x5<<<10 + (xor(x8,x7,x6<<<10) + x4<<<10 )<<s[9]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[1][i], xd[1][i], xv[0][i], xd[0][i], zero[i], zero[i],
			tmpv[6][i], tmpd[6][i]);
	}
	rotateState(model, s[9], tmpv[6], tmpd[6], rxv[2], rxd[2],
		tmpv[7][0], tmpd[7][0], highXV[2], highXD[2]);
	expandState(model, 0, rxv[2], rxd[2], xv[2], xd[2], tmpv[7], tmpd[7]);
		
	//x10 = x6<<<10 + (xor(x9,x8,x7<<<10) + x5<<<10 )<<s[10]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[2][i], xd[2][i], xv[1][i], xd[1][i], xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			tmpv[9][i], tmpd[9][i]);
	}
	rotateState(model, s[10], tmpv[9], tmpd[9], rxv[3], rxd[3],
		tmpv[10][0], tmpd[10][0], highXV[3], highXD[3]);
	expandState(model, 0, rxv[3], rxd[3], xv[3], xd[3], tmpv[10], tmpd[10]);

	//x11 = x7<<<10 + (xor(x10,x9,x8<<<10) + x6<<<10 )<<s[11]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[3][i], xd[3][i], xv[2][i], xd[2][i], xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			tmpv[12][i], tmpd[12][i]);
	}
	rotateState(model, s[11], tmpv[12], tmpd[12], rxv[4], rxd[4],
		tmpv[13][0], tmpd[13][0], highXV[4], highXD[4]);
	//add rxv[4] and xv[0]<<<10, the sum is tmpv[14]
	for (int i = 0; i < 32; i++) {
		addition(model, rxv[4][i], rxd[4][i],
			xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			tmpv[13][i], tmpd[13][i],
			tmpv[14][i], tmpd[14][i],
			tmpv[13][i + 1], tmpd[13][i + 1]);
	}
	//expand tmpv[14] and get xv[4]
	model.addConstr(tmpv[15][0] == 0);
	model.addConstr(tmpd[15][0] == 0);//no carry for LSB
	expandState(model, 0, tmpv[14], tmpd[14], xv[4], xd[4], tmpv[15], tmpd[15]);

	//x12 = x8<<<10 + (xor(x11,x10,x9<<<10) + x7<<<10 )<<s[12]
	vector<vector<GRBVar> > interV, interD;
	interV.resize(5);
	interD.resize(5);
	for (int i = 0; i < 5; i++) {
		interV[i].resize(33);
		interD[i].resize(33);
		for (int j = 0; j < 33; j++) {
			interV[i][j]=model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			interD[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[4][i], xd[4][i], xv[3][i], xd[3][i], xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			interV[0][i], interD[0][i]);
	}
	model.addConstr(interV[1][0] == 0);
	model.addConstr(interD[1][0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, interV[0][i], interD[0][i],
			xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			interV[1][i], interD[1][i],
			interV[2][i], interD[2][i],
			interV[1][i + 1], interD[1][i + 1]);
	}
	//rotate interV[2] to interV[3]
	rotateState(model, s[12],
		interV[2], interD[2],
		interV[3], interD[3],
		interV[4][0], interD[4][0],
		highXV[5], highXD[5]);
	//addition
	for (int i = 0; i < 32; i++) {
		addition(model, interV[3][i], interD[3][i],
			xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			interV[4][i], interD[4][i],
			zero[i], zero[i],
			interV[4][i + 1], interD[4][i + 1]);
	}

	//x13 = x9<<<10 + (xor(x12,x11,x10<<<10) + x8<<<10 )<<s[13]
	for (int i = 0; i < 32; i++) {
		XOR(model, zero[i], zero[i], xv[4][i], xd[4][i], xv[3][(i + 22) % 32], xd[3][(i + 22) % 32],
			tmpv[17][i], tmpd[17][i]);
	}
	model.addConstr(tmpv[18][0] == 0);
	model.addConstr(tmpd[18][0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[17][i], tmpd[17][i],
			xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			tmpv[18][i], tmpd[18][i],
			tmpv[19][i], tmpd[19][i],
			tmpv[18][i + 1], tmpd[18][i + 1]);
	}
	//rotate tmpv[19] to tmpv[20]
	rotateState(model, s[13],
		tmpv[19], tmpd[19],
		tmpv[20], tmpd[20],
		tmpv[21][0], tmpd[21][0],
		highXV[6], highXD[6]);
	//addition
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[20][i], tmpd[20][i],
			xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			tmpv[21][i], tmpd[21][i],
			zero[i], zero[i],
			tmpv[21][i + 1], tmpd[21][i + 1]);
	}

	//x14 = x10<<<10 + (xor(x13,x12,x11<<<10) + x9<<<10 )<<s[14]
	for (int i = 0; i < 32; i++) {
		XOR(model, zero[i], zero[i], zero[i], zero[i], xv[4][(i + 22) % 32], xd[4][(i + 22) % 32],
			tmpv[22][i], tmpd[22][i]);
	}
	model.addConstr(tmpv[23][0] == 0);
	model.addConstr(tmpd[23][0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[22][i], tmpd[22][i],
			xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			tmpv[23][i], tmpd[23][i],
			tmpv[24][i], tmpd[24][i],
			tmpv[23][i + 1], tmpd[23][i + 1]);
	}
	//rotate tmpv[24] to tmpv[25]
	rotateState(model, s[14],
		tmpv[24], tmpd[24],
		tmpv[25], tmpd[25],
		tmpv[26][0], tmpd[26][0],
		highXV[7], highXD[7]);
	//addition
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[25][i],tmpd[25][i],
			xv[3][(i + 22) % 32], xd[3][(i + 22) % 32],
			tmpv[26][i], tmpd[26][i],
			zero[i], zero[i],
			tmpv[26][i + 1], tmpd[26][i + 1]);
	}

	//x15 = x11<<<10 + (xor(x14,x13,x12<<<10) + x10<<<10 )<<s[15]
	rotateStateShift10(model, s[15], xv[3], xd[3], rxv[8], rxd[8],
		tmpv[27][0], tmpd[27][0], highXV[8], highXD[8]);
	for (int i = 0; i < 32; i++) {
		addition(model, rxv[8][i], rxd[8][i],
			xv[4][(i+22)%32], xd[4][(i+22)%32],
			tmpv[27][i], tmpd[27][i],
			tmpv[28][i], tmpd[28][i],
			tmpv[27][i + 1], tmpd[27][i + 1]);
	}
	for (int i = 0; i < 32; i++) {
		model.addConstr(tmpv[28][i] == 0);
		model.addConstr(tmpd[28][i] == 0);
	}


	//x11<<<10 + m7 = 0 (x16)
	//expan m7 and then assign it to x11
	model.addConstr(tmpv[29][0] == 0);
	model.addConstr(tmpd[29][0] == 0);//no carry for LSB
	expandState(model, 22, mvInv, mdInv, xv[4], xd[4],
		tmpv[29], tmpd[29]);

	GRBLinExpr wei = 0;
	for (int j = 1; j < 4; j++) {
		for (int i = 0; i < 32; i++) {
			wei += xd[j][i];
		}
	}

	loadConstraint(model, input[0], xv[0], xd[0]);
	loadConstraint(model, input[1], xv[1], xd[1]);
	loadConstraint(model, input[2], xv[2], xd[2]);
	loadConstraint(model, input[3], xv[3], xd[3]);
	loadConstraint(model, input[4], xv[4], xd[4]);

	//model.setObjective(wei, GRB_MINIMIZE);
	model.optimize();

	string str;
	for (int i = 0; i < 5; i++) {
		getSignedString(xv[i], xd[i], str);
		cout << "X" << (i + 7) << ": " << str << endl;
	}
	
	getSignedString(tmpv[3], tmpd[3], str);
	cout << "xor0:" << str << endl;
	getSignedString(tmpv[6], tmpd[6], str);
	cout << "xor1:" << str << endl;
	getSignedString(tmpv[9], tmpd[9], str);
	cout << "xor2:" << str << endl;
	getSignedString(tmpv[12], tmpd[12], str);
	cout << "xor3:" << str << endl;
	getSignedString(interV[0], interD[0], str);
	cout << "xor4:" << str << endl;
	getSignedString(tmpv[17], tmpd[17], str);
	cout << "xor5:" << str << endl;
	getSignedString(tmpv[22], tmpd[22], str);
	cout << "xor6:" << str << endl;

	//check
	cout << endl;
	getSignedString(rxv[4], rxd[4], str);
	cout << "rot11:" << str << endl;
	getSignedString(tmpv[14], tmpd[14], str);
	cout << "sum11:" << str << endl;
}

void Primitive::model_RIPEMD160_FirstRound_Left_36Collision_All() {
	//X[i+5] = X[i+1]<<<10 + (XOR(X[i+4],X[i+3],X[i+2]<<<10) + X[i]<<<10 + m)<<s
	int s[32] = {
		11,14,15,12,5,8,7,9,11,13,14,15,6,7,9,8,
		7,6,8,13,11,9,7,15,7,12,15,9,11,7,13,12
	};

	string m0 =    "===========u====================";
	string m6 =    "=======================n========";
	string m9Inv = "===============================n";//n
	
	string input[10] = {
		"????????????????????????????????",
		"????????????????????????????????",
		"????????????????????????????????",
		"????????????????????????????????",
		"????????????????????????????????",
	};

	GRBEnv env = GRBEnv();
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);
	vector<vector<GRBVar> > xv, xd, tmpv, tmpd, rxv, rxd;
	vector<GRBVar> mv, md, mvInv, mdInv,m6v,m6d;
	vector<GRBVar> midXV, midXD, highXV, highXD;

	vector<GRBVar> zero;
	zero.resize(32);
	for (int i = 0; i < 32; i++) {
		zero[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(zero[i] == 0);
	}

	xv.resize(10);
	xd.resize(10);
	rxv.resize(10);
	rxd.resize(10);

	midXD.resize(10);
	midXV.resize(10);
	highXD.resize(10);
	highXV.resize(10);

	tmpv.resize(30);
	tmpd.resize(30);

	mv.reserve(32);
	md.resize(32);
	mvInv.resize(32);
	mdInv.resize(32);
	m6v.reserve(32);
	m6d.resize(32);

	for (int i = 0; i < 10; i++) {
		xv[i].resize(32);
		xd[i].resize(32);
		rxv[i].resize(32);
		rxd[i].resize(32);

		midXV[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		midXD[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		highXV[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		highXD[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

		for (int j = 0; j < 32; j++) {
			xv[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			xd[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			rxv[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			rxd[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		for (int k = 0; k < 3; k++) {
			tmpv[3 * i + k].resize(33);
			tmpd[3 * i + k].resize(33);
		}

		for (int j = 0; j < 33; j++) {
			for (int k = 0; k < 3; k++) {
				tmpv[3 * i + k][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				tmpd[3 * i + k][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
	}

	for (int i = 0; i < 32; i++) {
		mv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		md[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		m6v[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		m6d[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		mvInv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		mdInv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	}
	loadConstraint(model, m0, mv, md);
	loadConstraint(model, m6, m6v, m6d);
	loadConstraint(model, m9Inv, mvInv, mdInv);

	cout << "after loading" << endl;

	//x7 = (m7+0)<<s[7]
	//first, rotate m7; second, expand it and assign it to x7
	rotateState(model, s[0], mv, md, rxv[0], rxd[0],
		tmpv[0][0], tmpd[0][0], highXV[0], highXD[0]);
	expandState(model, 0, rxv[0], rxd[0], xv[0], xd[0], tmpv[0], tmpd[0]);

	//determine x8,x9,x10 by propagating signed differences
	//x8 = x4<<<10 + (xor(x7,x6,x5<<<10) + x3<<<10 )<<s[8]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[0][i], xd[0][i], zero[i], zero[i], zero[i], zero[i],
			tmpv[3][i], tmpd[3][i]);
	}
	rotateState(model, s[1], tmpv[3], tmpd[3], rxv[1], rxd[1],
		tmpv[4][0], tmpd[4][0], highXV[1], highXD[1]);
	expandState(model, 0, rxv[1], rxd[1], xv[1], xd[1], tmpv[4], tmpd[4]);

	

	//x9 = x5<<<10 + (xor(x8,x7,x6<<<10) + x4<<<10 )<<s[9]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[1][i], xd[1][i], xv[0][i], xd[0][i], zero[i], zero[i],
			tmpv[6][i], tmpd[6][i]);
	}
	rotateState(model, s[2], tmpv[6], tmpd[6], rxv[2], rxd[2],
		tmpv[7][0], tmpd[7][0], highXV[2], highXD[2]);
	expandState(model, 0, rxv[2], rxd[2], xv[2], xd[2], tmpv[7], tmpd[7]);

	//x10 = x6<<<10 + (xor(x9,x8,x7<<<10) + x5<<<10 )<<s[10]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[2][i], xd[2][i], xv[1][i], xd[1][i], xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			tmpv[9][i], tmpd[9][i]);
	}
	rotateState(model, s[3], tmpv[9], tmpd[9], rxv[3], rxd[3],
		tmpv[10][0], tmpd[10][0], highXV[3], highXD[3]);
	expandState(model, 0, rxv[3], rxd[3], xv[3], xd[3], tmpv[10], tmpd[10]);
	

	//x11 = x7<<<10 + (xor(x10,x9,x8<<<10) + x6<<<10 )<<s[11]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[3][i], xd[3][i], xv[2][i], xd[2][i], xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			tmpv[12][i], tmpd[12][i]);
	}

	//new way
	/*model.addConstr(tmpv[15][0] == 0);
	model.addConstr(tmpd[15][0] == 0);//no carry for LSB
	expandState(model, 0, tmpv[12], tmpd[12], rxv[4], rxd[4], tmpv[15], tmpd[15]);
	model.addConstr(tmpv[13][0] == 0);
	model.addConstr(tmpd[13][0] == 0);//no carry for LSB
	cout << "test" << endl;
	for (int i = 0; i < 32; i++) {
		addition(model, rxv[4][(i+32-s[4])%32], rxd[4][(i + 32 - s[4]) % 32],
			xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			tmpv[13][i], tmpd[13][i],
			xv[4][i], xd[4][i],
			tmpv[13][i + 1], tmpd[13][i + 1]);
	}*/

	rotateState(model, s[4], tmpv[12], tmpd[12], rxv[4], rxd[4],
		tmpv[13][0], tmpd[13][0], highXV[4], highXD[4]);
	model.addConstr(tmpd[13][0] == 1);
	model.addConstr(highXD[4] == 1);
	//add rxv[4] and xv[0]<<<10, the sum is tmpv[14]
	for (int i = 0; i < 32; i++) {
		addition(model, rxv[4][i], rxd[4][i],
			xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			tmpv[13][i], tmpd[13][i],
			tmpv[14][i], tmpd[14][i],
			tmpv[13][i + 1], tmpd[13][i + 1]);
	}
	//expand tmpv[14] and get xv[4]
	model.addConstr(tmpv[15][0] == 0);
	model.addConstr(tmpd[15][0] == 0);//no carry for LSB
	expandState(model, 0, tmpv[14], tmpd[14], xv[4], xd[4], tmpv[15], tmpd[15]);

	//x12 = x8<<<10 + (xor(x11,x10,x9<<<10) + x7<<<10 )<<s[12]
	vector<vector<GRBVar> > interV, interD;
	interV.resize(5);
	interD.resize(5);
	for (int i = 0; i < 5; i++) {
		interV[i].resize(33);
		interD[i].resize(33);
		for (int j = 0; j < 33; j++) {
			interV[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			interD[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[4][i], xd[4][i], xv[3][i], xd[3][i], xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			interV[0][i], interD[0][i]);
	}
	model.addConstr(interV[1][0] == 0);
	model.addConstr(interD[1][0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, interV[0][i], interD[0][i],
			xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			interV[1][i], interD[1][i],
			interV[2][i], interD[2][i],
			interV[1][i + 1], interD[1][i + 1]);
	}
	//rotate interV[2] to interV[3]
	rotateState(model, s[5],
		interV[2], interD[2],
		interV[3], interD[3],
		interV[4][0], interD[4][0],
		highXV[5], highXD[5]);
	//addition
	for (int i = 0; i < 32; i++) {
		addition(model, interV[3][i], interD[3][i],
			xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			interV[4][i], interD[4][i],
			zero[i], zero[i],
			interV[4][i + 1], interD[4][i + 1]);
	}

	//x13 = x9<<<10 + (xor(x12,x11,x10<<<10) + x8<<<10 )<<s[13]
	for (int i = 0; i < 32; i++) {
		XOR(model, zero[i], zero[i], xv[4][i], xd[4][i], xv[3][(i + 22) % 32], xd[3][(i + 22) % 32],
			tmpv[17][i], tmpd[17][i]);
	}

	vector<GRBVar> extraV, extraD, extraCV, extraCD;
	extraV.resize(32);
	extraD.resize(32);
	extraCV.resize(33);
	extraCD.resize(33);
	initializeVar(model, extraV, 32);
	initializeVar(model, extraD, 32);
	initializeVar(model, extraCV, 33);
	initializeVar(model, extraCD, 33);

	model.addConstr(extraCV[0] == 0);
	model.addConstr(extraCD[0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[17][i], tmpd[17][i],
			m6v[i], m6d[i],
			extraCV[i], extraCD[i],
			extraV[i], extraD[i],
			extraCV[i + 1], extraCD[i + 1]);
	}

	model.addConstr(tmpv[18][0] == 0);
	model.addConstr(tmpd[18][0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, extraV[i], extraD[i],
			xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			tmpv[18][i], tmpd[18][i],
			tmpv[19][i], tmpd[19][i],
			tmpv[18][i + 1], tmpd[18][i + 1]);
	}
	//rotate tmpv[19] to tmpv[20]
	rotateState(model, s[6],
		tmpv[19], tmpd[19],
		tmpv[20], tmpd[20],
		tmpv[21][0], tmpd[21][0],
		highXV[6], highXD[6]);
	//addition
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[20][i], tmpd[20][i],
			xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			tmpv[21][i], tmpd[21][i],
			zero[i], zero[i],
			tmpv[21][i + 1], tmpd[21][i + 1]);
	}

	//x14 = x10<<<10 + (xor(x13,x12,x11<<<10) + x9<<<10 )<<s[14]
	for (int i = 0; i < 32; i++) {
		XOR(model, zero[i], zero[i], zero[i], zero[i], xv[4][(i + 22) % 32], xd[4][(i + 22) % 32],
			tmpv[22][i], tmpd[22][i]);
	}
	model.addConstr(tmpv[23][0] == 0);
	model.addConstr(tmpd[23][0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[22][i], tmpd[22][i],
			xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			tmpv[23][i], tmpd[23][i],
			tmpv[24][i], tmpd[24][i],
			tmpv[23][i + 1], tmpd[23][i + 1]);
	}
	//rotate tmpv[24] to tmpv[25]
	rotateState(model, s[7],
		tmpv[24], tmpd[24],
		tmpv[25], tmpd[25],
		tmpv[26][0], tmpd[26][0],
		highXV[7], highXD[7]);
	//addition
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[25][i], tmpd[25][i],
			xv[3][(i + 22) % 32], xd[3][(i + 22) % 32],
			tmpv[26][i], tmpd[26][i],
			zero[i], zero[i],
			tmpv[26][i + 1], tmpd[26][i + 1]);
	}

	//x15 = x11<<<10 + (xor(x14,x13,x12<<<10) + x10<<<10 )<<s[15]
	rotateStateShift10(model, s[8], xv[3], xd[3], rxv[8], rxd[8],
		tmpv[27][0], tmpd[27][0], highXV[8], highXD[8]);
	for (int i = 0; i < 32; i++) {
		addition(model, rxv[8][i], rxd[8][i],
			xv[4][(i + 22) % 32], xd[4][(i + 22) % 32],
			tmpv[27][i], tmpd[27][i],
			tmpv[28][i], tmpd[28][i],
			tmpv[27][i + 1], tmpd[27][i + 1]);
	}
	for (int i = 0; i < 32; i++) {
		model.addConstr(tmpv[28][i] == 0);
		model.addConstr(tmpd[28][i] == 0);
	}

	//x11<<<10 + m7 = 0 (x16)
	//expan m7 and then assign it to x11
	model.addConstr(tmpv[29][0] == 0);
	model.addConstr(tmpd[29][0] == 0);//no carry for LSB
	expandState(model, 22, mvInv, mdInv, xv[4], xd[4],
		tmpv[29], tmpd[29]);

	GRBLinExpr wei = 0;
	for (int j = 0; j < 5; j++) {
		if (j == 0 || j == 4) {
			for (int i = 0; i < 32; i++) {
				wei += xd[j][i];
			}
		}
	}

	model.setObjective(wei, GRB_MINIMIZE);
	model.optimize();

	string str;
	for (int i = 0; i < 5; i++) {
		getSignedString(xv[i], xd[i], str);
		cout << "X" << (i + 7) << ": " << str << endl;
	}

	getSignedString(tmpv[3], tmpd[3], str);
	cout << "xor0:" << str << endl;
	getSignedString(tmpv[6], tmpd[6], str);
	cout << "xor1:" << str << endl;
	getSignedString(tmpv[9], tmpd[9], str);
	cout << "xor2:" << str << endl;
	getSignedString(tmpv[12], tmpd[12], str);
	cout << "xor3:" << str << endl;
	//getSignedString(interV[0], interD[0], str);
	//cout << "xor4:" << str << endl;
	getSignedString(tmpv[17], tmpd[17], str);
	cout << "xor5:" << str << endl;
	getSignedString(tmpv[22], tmpd[22], str);
	cout << "xor6:" << str << endl;
}

void Primitive::model_RIPEMD160_FirstRound_Left_36Collision(int sh[],int length) {
	//X[i+5] = X[i+1]<<<10 + (XOR(X[i+4],X[i+3],X[i+2]<<<10) + X[i]<<<10 + m)<<s
	int s[32] = {
		11,14,15,12,5,8,7,9,11,13,14,15,6,7,9,8,
		7,6,8,13,11,9,7,15,7,12,15,9,11,7,13,12
	};

	string input[5] = {
		"uuuuuuuuuuuuuuuuuunuuuuuuuuuuuu=",
		"=====u==nn===ununnuu=========n=u",
		"uuuuu=un==nn=ununn==u=uunnuu=n==",
		"===nuu=======nuuuuu=============",
		"nnnnnnnn==unnnnnnnnnnnnnnnnunnnn",
	};

	GRBEnv env = GRBEnv();
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);
	vector<vector<GRBVar> > xv, xd, tmpv, tmpd, rxv, rxd;
	vector<GRBVar> mv, md, mvInv, mdInv, m6v, m6d;
	vector<GRBVar> midXV, midXD, highXV, highXD;

	vector<GRBVar> zero;
	zero.resize(32);
	for (int i = 0; i < 32; i++) {
		zero[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(zero[i] == 0);
	}

	xv.resize(10);
	xd.resize(10);
	rxv.resize(10);
	rxd.resize(10);

	midXD.resize(10);
	midXV.resize(10);
	highXD.resize(10);
	highXV.resize(10);

	tmpv.resize(30);
	tmpd.resize(30);

	mv.reserve(32);
	md.resize(32);
	mvInv.resize(32);
	mdInv.resize(32);
	m6v.reserve(32);
	m6d.resize(32);

	vector<vector<GRBVar> > xc(6);//condition
	initializeVar(model, xc, 32);

	for (int i = 0; i < 10; i++) {
		xv[i].resize(32);
		xd[i].resize(32);
		rxv[i].resize(32);
		rxd[i].resize(32);

		midXV[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		midXD[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		highXV[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		highXD[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

		for (int j = 0; j < 32; j++) {
			xv[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			xd[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			rxv[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			rxd[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		for (int k = 0; k < 3; k++) {
			tmpv[3 * i + k].resize(33);
			tmpd[3 * i + k].resize(33);
		}

		for (int j = 0; j < 33; j++) {
			for (int k = 0; k < 3; k++) {
				tmpv[3 * i + k][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				tmpd[3 * i + k][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
	}

	for (int i = 0; i < 32; i++) {
		mv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		md[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		m6v[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		m6d[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		mvInv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		mdInv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	}

	//loadConstraint(model, m0, mv, md);
	//loadConstraint(model, m6, m6v, m6d);
	//loadConstraint(model, m9Inv, mvInv, mdInv);

	for (int i = 0; i < 5; i++) {
		loadConstraint(model, input[i], xv[i], xd[i]);
	}
	
	bool isUsed[32];
	for (int i = 0; i < 32; i++) {
		isUsed[i] = false;
	}
	for (int i = 0; i < length; i++) {
		isUsed[sh[i]] = 1;
	}
	for (int i = 0; i < 32; i++) {
		if (isUsed[i]) {
			model.addConstr(mv[i] == 0);
			model.addConstr(md[i] == 1);
			model.addConstr(m6v[(i + 25) % 32] == 1);
			model.addConstr(m6d[(i + 25) % 32] == 1);
			model.addConstr(mvInv[(i + 12) % 32] == 1);
			model.addConstr(mdInv[(i + 12) % 32] == 1);
		}
		else {
			model.addConstr(mv[i] == 0);
			model.addConstr(md[i] == 0);
			model.addConstr(m6v[(i + 25) % 32] == 0);
			model.addConstr(m6d[(i + 25) % 32] == 0);
			model.addConstr(mvInv[(i + 12) % 32] == 0);
			model.addConstr(mdInv[(i + 12) % 32] == 0);
		}
	}

	//x7 = (m7+0)<<s[7]
	//expand it and assign it to x7
	model.addConstr(tmpv[0][0] == 0);
	model.addConstr(tmpd[0][0] == 0);
	//expandState(model, 0, mv, md, rxv[0], rxd[0], tmpv[0], tmpd[0]);
	zeroState(model, mv, md, rxv[0], rxd[0], tmpv[0], tmpd[0]);
	model.addConstr(0 == tmpv[1][0]);
	model.addConstr(0 == tmpd[1][0]);
	addExpandState(model, 0, s[0], zero, zero, 
		rxv[0], rxd[0], xv[0], xd[0], tmpv[1], tmpd[1]);


	//determine x8,x9,x10 by propagating signed differences
	//x8 = x4<<<10 + (xor(x7,x6,x5<<<10) + x3<<<10 )<<s[8]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[0][i], xd[0][i], zero[i], zero[i], zero[i], zero[i],
			tmpv[3][i], tmpd[3][i]);
	}
	model.addConstr(tmpv[4][0] == 0);
	model.addConstr(tmpd[4][0] == 0);
	//expandState(model, 0, tmpv[3], tmpd[3], rxv[1], rxd[1], tmpv[4], tmpd[4]);
	zeroState(model, tmpv[3], tmpd[3], rxv[1], rxd[1], tmpv[4], tmpd[4]);
	model.addConstr(0 == tmpv[5][0]);
	model.addConstr(0 == tmpd[5][0]);
	addExpandState(model, 0, s[1], zero, zero,
		rxv[1], rxd[1], xv[1], xd[1], tmpv[5], tmpd[5]);


	/*//x9 = x5<<<10 + (xor(x8,x7,x6<<<10) + x4<<<10 )<<s[9]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[1][i], xd[1][i], xv[0][i], xd[0][i], zero[i], zero[i],
			tmpv[6][i], tmpd[6][i]);
	}
	rotateState(model, s[2], tmpv[6], tmpd[6], rxv[2], rxd[2],
		tmpv[7][0], tmpd[7][0], highXV[2], highXD[2]);
	expandState(model, 0, rxv[2], rxd[2], xv[2], xd[2], tmpv[7], tmpd[7]);

	//x10 = x6<<<10 + (xor(x9,x8,x7<<<10) + x5<<<10 )<<s[10]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[2][i], xd[2][i], xv[1][i], xd[1][i], xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			tmpv[9][i], tmpd[9][i]);
	}
	rotateState(model, s[3], tmpv[9], tmpd[9], rxv[3], rxd[3],
		tmpv[10][0], tmpd[10][0], highXV[3], highXD[3]);
	expandState(model, 0, rxv[3], rxd[3], xv[3], xd[3], tmpv[10], tmpd[10]);
	*/

	//x9 = x5<<<10 + (xor(x8,x7,x6<<<10) + x4<<<10 )<<s[9]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[1][i], xd[1][i], xv[0][i], xd[0][i], zero[i], zero[i],
			tmpv[6][i], tmpd[6][i]);
	}
	model.addConstr(tmpv[7][0] == 0);
	model.addConstr(tmpd[7][0] == 0);
	//expandState(model, 0, tmpv[6], tmpd[6], rxv[2], rxd[2], tmpv[7], tmpd[7]);
	zeroState(model, tmpv[6], tmpd[6], rxv[2], rxd[2], tmpv[7], tmpd[7]);
	model.addConstr(0 == tmpv[8][0]);
	model.addConstr(0 == tmpd[8][0]);
	addExpandState(model, 0, s[2], zero, zero,
		rxv[2], rxd[2], xv[2], xd[2], tmpv[8], tmpd[8]);

	//x10 = x6<<<10 + (xor(x9,x8,x7<<<10) + x5<<<10 )<<s[10]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[2][i], xd[2][i], xv[1][i], xd[1][i], xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			tmpv[9][i], tmpd[9][i]);
		//XORCut
		XORCut(model, xv[2][i], xd[2][i], xv[1][i], xd[1][i], xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			tmpv[9][i], tmpd[9][i],
			xc[2][i],xc[1][i], xc[0][(i + 22) % 32]);
	}
	model.addConstr(tmpv[10][0] == 0);
	model.addConstr(tmpd[10][0] == 0);
	//expandState(model, 0, tmpv[9], tmpd[9], rxv[3], rxd[3], tmpv[10], tmpd[10]);
	zeroState(model, tmpv[9], tmpd[9], rxv[3], rxd[3], tmpv[10], tmpd[10]);
	model.addConstr(0 == tmpv[11][0]);
	model.addConstr(0 == tmpd[11][0]);
	addExpandState(model, 0, s[3], zero, zero,
		rxv[3], rxd[3], xv[3], xd[3], tmpv[11], tmpd[11]);

	//x11 = x7<<<10 + (xor(x10,x9,x8<<<10) + x6<<<10 )<<s[11]
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[3][i], xd[3][i], xv[2][i], xd[2][i], xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			tmpv[12][i], tmpd[12][i]);
		//XORCut
		XORCut(model, xv[3][i], xd[3][i], xv[2][i], xd[2][i], xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			tmpv[12][i], tmpd[12][i],
			xc[3][i], xc[2][i], xc[1][(i + 22) % 32]);
	}

	//new way
	model.addConstr(tmpv[15][0] == 0);
	model.addConstr(tmpd[15][0] == 0);//no carry for LSB
	//expandState(model, 0, tmpv[12], tmpd[12], rxv[4], rxd[4], tmpv[15], tmpd[15]);
	zeroState(model, tmpv[12], tmpd[12], rxv[4], rxd[4], tmpv[15], tmpd[15]);
	model.addConstr(tmpv[13][0] == 0);
	model.addConstr(tmpd[13][0] == 0);//no carry for LSB
	addExpandState(model, 10, s[4], xv[0], xd[0],
		rxv[4], rxd[4], xv[4], xd[4], tmpv[13], tmpd[13]);

	//x12 = x8<<<10 + (xor(x11,x10,x9<<<10) + x7<<<10 )<<s[12]
	vector<vector<GRBVar> > interV, interD;
	interV.resize(5);
	interD.resize(5);
	for (int i = 0; i < 5; i++) {
		interV[i].resize(33);
		interD[i].resize(33);
		for (int j = 0; j < 33; j++) {
			interV[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			interD[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
	for (int i = 0; i < 32; i++) {
		XOR(model, xv[4][i], xd[4][i], xv[3][i], xd[3][i], xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			interV[0][i], interD[0][i]);
		XORCut(model, xv[4][i], xd[4][i], xv[3][i], xd[3][i], xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			interV[0][i], interD[0][i],
			xc[4][i],xc[3][i],xc[2][(i + 22) % 32]);
	}
	model.addConstr(interV[1][0] == 0);
	model.addConstr(interD[1][0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, interV[0][i], interD[0][i],
			xv[0][(i + 22) % 32], xd[0][(i + 22) % 32],
			interV[1][i], interD[1][i],
			interV[2][i], interD[2][i],
			interV[1][i + 1], interD[1][i + 1]);
	}
	//rotate interV[2] to interV[3]
	rotateState(model, s[5],
		interV[2], interD[2],
		interV[3], interD[3],
		interV[4][0], interD[4][0],
		highXV[5], highXD[5]);
	//addition
	for (int i = 0; i < 32; i++) {
		addition(model, interV[3][i], interD[3][i],
			xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			interV[4][i], interD[4][i],
			zero[i], zero[i],
			interV[4][i + 1], interD[4][i + 1]);
	}

	//x13 = x9<<<10 + (xor(x12,x11,x10<<<10) + x8<<<10 )<<s[13]
	for (int i = 0; i < 32; i++) {
		XOR(model, zero[i], zero[i], xv[4][i], xd[4][i], xv[3][(i + 22) % 32], xd[3][(i + 22) % 32],
			tmpv[17][i], tmpd[17][i]);
	}

	vector<GRBVar> extraV, extraD, extraCV, extraCD;
	extraV.resize(32);
	extraD.resize(32);
	extraCV.resize(33);
	extraCD.resize(33);
	initializeVar(model, extraV, 32);
	initializeVar(model, extraD, 32);
	initializeVar(model, extraCV, 33);
	initializeVar(model, extraCD, 33);

	model.addConstr(extraCV[0] == 0);
	model.addConstr(extraCD[0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[17][i], tmpd[17][i],
			m6v[i], m6d[i],
			extraCV[i], extraCD[i],
			extraV[i], extraD[i],
			extraCV[i + 1], extraCD[i + 1]);
	}

	model.addConstr(tmpv[18][0] == 0);
	model.addConstr(tmpd[18][0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, extraV[i], extraD[i],
			xv[1][(i + 22) % 32], xd[1][(i + 22) % 32],
			tmpv[18][i], tmpd[18][i],
			tmpv[19][i], tmpd[19][i],
			tmpv[18][i + 1], tmpd[18][i + 1]);
	}
	//rotate tmpv[19] to tmpv[20]
	rotateState(model, s[6],
		tmpv[19], tmpd[19],
		tmpv[20], tmpd[20],
		tmpv[21][0], tmpd[21][0],
		highXV[6], highXD[6]);
	//addition
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[20][i], tmpd[20][i],
			xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			tmpv[21][i], tmpd[21][i],
			zero[i], zero[i],
			tmpv[21][i + 1], tmpd[21][i + 1]);
	}

	//x14 = x10<<<10 + (xor(x13,x12,x11<<<10) + x9<<<10 )<<s[14]
	for (int i = 0; i < 32; i++) {
		XOR(model, zero[i], zero[i], zero[i], zero[i], xv[4][(i + 22) % 32], xd[4][(i + 22) % 32],
			tmpv[22][i], tmpd[22][i]);
	}
	model.addConstr(tmpv[23][0] == 0);
	model.addConstr(tmpd[23][0] == 0);
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[22][i], tmpd[22][i],
			xv[2][(i + 22) % 32], xd[2][(i + 22) % 32],
			tmpv[23][i], tmpd[23][i],
			tmpv[24][i], tmpd[24][i],
			tmpv[23][i + 1], tmpd[23][i + 1]);
	}
	//rotate tmpv[24] to tmpv[25]
	rotateState(model, s[7],
		tmpv[24], tmpd[24],
		tmpv[25], tmpd[25],
		tmpv[26][0], tmpd[26][0],
		highXV[7], highXD[7]);
	//addition
	for (int i = 0; i < 32; i++) {
		addition(model, tmpv[25][i], tmpd[25][i],
			xv[3][(i + 22) % 32], xd[3][(i + 22) % 32],
			tmpv[26][i], tmpd[26][i],
			zero[i], zero[i],
			tmpv[26][i + 1], tmpd[26][i + 1]);
	}

	//x15 = x11<<<10 + (xor(x14,x13,x12<<<10) + x10<<<10 )<<s[15]
	rotateStateShift10(model, s[8], xv[3], xd[3], rxv[8], rxd[8],
		tmpv[27][0], tmpd[27][0], highXV[8], highXD[8]);
	for (int i = 0; i < 32; i++) {
		addition(model, rxv[8][i], rxd[8][i],
			xv[4][(i + 22) % 32], xd[4][(i + 22) % 32],
			tmpv[27][i], tmpd[27][i],
			tmpv[28][i], tmpd[28][i],
			tmpv[27][i + 1], tmpd[27][i + 1]);
	}
	for (int i = 0; i < 32; i++) {
		model.addConstr(tmpv[28][i] == 0);
		model.addConstr(tmpd[28][i] == 0);
	}

	//x11<<<10 + m7 = 0 (x16)
	//expan m7 and then assign it to x11
	model.addConstr(tmpv[29][0] == 0);
	model.addConstr(tmpd[29][0] == 0);//no carry for LSB
	expandState(model, 22, mvInv, mdInv, xv[4], xd[4],
		tmpv[29], tmpd[29]);
		

	GRBLinExpr wei = 0;
	for (int j = 0; j < 5; j++) {
		if (j == 0 || j == 4) {
			for (int i = 0; i < 32; i++) {
				wei += xd[j][i];
			}
		}
	}

	//string x0 = "u==nuuuuuuuuuuuuuuuuuuuuuuuuuuuu";
	//string x4 = "nnnnnnnnnn=unnnnnnnnnnnnnnnnnnnn";
	//loadConstraint(model, x0, xv[0], xd[0]);
	//loadConstraint(model, x4, xv[4], xd[4]);

	string xString[5] = {
		"n==========unnnnnnnnn=nuuuuuuuuu",
		"n====n=u=uu===un=uu========un==n",
		"=nun=n=u=uunun==n==nuuu==uu==nu=",
		"==============nuuuuuu===========",
		"nnnnnnnnnn=unnnnnnnnnnnnnnnnnnnn",
	};
	for (int i = 0; i < 5; i++) {
		//setHintValue(model, xString[i], xv[i], xd[i]);
		//if (i == 0 || i == 4) {
			//loadConstraint(model, xString[i], xv[i], xd[i]);
		//}
	}
	

	//model.setObjective(wei, GRB_MINIMIZE);
	model.optimize();

	string str;
	for (int i = 0; i < 5; i++) {
		getSignedString(xv[i], xd[i], str);
		cout << "X" << (i + 7) << ": \"" << str <<"\","<< endl;
	}

	getSignedString(tmpv[3], tmpd[3], str);
	cout << "\"" << str << "\"," << endl;
	getSignedString(tmpv[6], tmpd[6], str);
	cout << "\"" << str << "\"" << endl;
	getSignedString(tmpv[9], tmpd[9], str);
	cout << "\"" << str << "\"," << endl;
	getSignedString(tmpv[12], tmpd[12], str);
	cout << "\"" << str << "\"," << endl;
	getSignedString(interV[0], interD[0], str);
	cout << "\"" << str << "\"," << endl;
	getSignedString(tmpv[17], tmpd[17], str);
	cout << "\"" << str << "\"," << endl;
	getSignedString(tmpv[22], tmpd[22], str);
	cout << "\"" << str << "\"," << endl;

	getSignedString(m6v, m6d, str);
	cout << "m6:" << str << endl;
}

int Primitive::model_RIPEMD160_FirstRound_Right_36Collision(int sh[],int length,GRBEnv& env) {
	int s[32] = {
		8,9,9,11,13,15,15,5,7,7,8,11,14,14,12,6,
		9,13,15,7,12,8,9,11,7,7,12,7,6,15,13,11,
	};

	int p[32] = {
		5,14,7,0,9,2,11,4,13,6,15,8,1,10,3,12,
		6,11,3,7,0,13,5,10,14,15,8,12,4,9,1,2
	};

	string msgDiff[16] = {
		"===========n====================",//m0[+20]
		"================================",//m1
		"================================",//m2
		"================================",//m3
		"================================",//m4
		"================================",//m5
		"==================u=============",//m6[-13]
		"================================",//m7
		"================================",//m8
		"===============================n",//m9[+0]
		"================================",//m10
		"================================",//m11
		"================================",//m12
		"================================",//m13
		"================================",//m14
		"================================",//m15
	};

	const int inputSize = 19;

	string input[inputSize] = {
		/*"===n=========n============n=====",//y11
		"=====================n==========",//y12
		"========u========n======u===u===",//y13
		"==u======n===n===========n======",//y14
		"======u============uu======n====",//y15*/

		/*"===============n================",
		"=======u========================",
		"================================",
		"==========n=============n=======",
		"===nu===========================",
		"================================",
		"====================u===========",
		"n===========nu==================",
		"=======u==================u=====",*/

		/*"????????????????????????????????",
		"????????????????????????????????",
		"????????????????????????????????",
		"????????????????????????????????",*/
		"????????????????????????????????",
		"????????????????????????????????",
		"????????????????????????????????",
		"????????????????????????????????",
		"????????????????????????????????",

		"================================",
		"================================",
		"=u============================n=",
		"==========u============n========",
		"====n============u==============",
		"================================",
		"===================nu===========",
		"n============u==================",
		"=======n==================u=====",

		"================================",//y25
		"================================",//y26
		"================================",//y27
		"================================",//y28
		"================================",//y29
	};

	//GRBEnv env;
	//env.set(GRB_IntParam_OutputFlag, 0);
	//env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);

	//condition
	vector<vector<GRBVar> > condition;//x0-x15
	condition.resize(inputSize);
	initializeVar(model, condition, 32);

	//zero
	vector<GRBVar> zerov, zerod;
	zerov.resize(32);
	zerod.resize(32);
	initializeVar(model, zerov, 32);
	initializeVar(model, zerod, 32);

	//state
	vector<vector<GRBVar> > xv, xd;
	xv.resize(inputSize);
	xd.resize(inputSize);
	initializeVar(model, xv, 32);
	initializeVar(model, xd, 32);

	//intermediate state (carry, rotate)
	vector<vector<GRBVar> > tv, td;
	tv.resize(inputSize * 6);
	td.resize(inputSize * 6);
	initializeVar(model, tv, 33);
	initializeVar(model, td, 33);

	//output of boolean functions
	vector<vector<GRBVar> > bOutV, bOutD;
	bOutV.resize(inputSize);
	bOutD.resize(inputSize);
	initializeVar(model, bOutV, 32);
	initializeVar(model, bOutD, 32);

	//diff after adding the msg difference
	vector<vector<GRBVar> > msgAV(inputSize), msgAD(inputSize);
	vector<vector<GRBVar> > msgCV(inputSize), msgCD(inputSize);
	vector<vector<GRBVar> > msgOV(inputSize), msgOD(inputSize);
	initializeVar(model, msgAV, 32);
	initializeVar(model, msgAD, 32);
	initializeVar(model, msgCV, 33);
	initializeVar(model, msgCD, 33);
	initializeVar(model, msgOV, 32);
	initializeVar(model, msgOD, 32);

	//diff of the message difference
	vector<vector<GRBVar> > msgV(16), msgD(16);
	initializeVar(model, msgV, 32);
	initializeVar(model, msgD, 32);

	vector<GRBVar> highXD, highXV;
	highXD.resize(32);
	highXV.resize(32);
	initializeVar(model, highXD, 32);
	initializeVar(model, highXV, 32);

	//loading constraints
	for (int i = 0; i < inputSize; i++) {
		loadConstraint(model, input[i], xv[i], xd[i]);
	}

	for (int i = 0; i < 16; i++) {
		if ((i != 0) && (i != 6) && (i != 9)) {
			loadConstraint(model, msgDiff[i], msgV[i], msgD[i]);
		}
	}

	bool isUsed[32];
	for (int i = 0; i < 32; i++) {
		isUsed[i] = false;
	}
	for (int i = 0; i < length; i++) {
		isUsed[sh[i]] = 1;
	}
	bool first = true;
	for (int i = 0; i < 32; i++) {
		if (isUsed[i] && first) {
			model.addConstr(msgV[0][i] == 0);
			model.addConstr(msgD[0][i] == 1);
			model.addConstr(msgV[6][(i + 25) % 32] == 1);
			model.addConstr(msgD[6][(i + 25) % 32] == 1);
			model.addConstr(msgV[9][(i + 12) % 32] == 0);
			model.addConstr(msgD[9][(i + 12) % 32] == 1);
			first = false;
		}
		else if (isUsed[i] && !first) {
			model.addConstr(msgV[0][i] == 1);
			model.addConstr(msgD[0][i] == 1);
			model.addConstr(msgV[6][(i + 25) % 32] == 0);
			model.addConstr(msgD[6][(i + 25) % 32] == 1);
			model.addConstr(msgV[9][(i + 12) % 32] == 1);
			model.addConstr(msgD[9][(i + 12) % 32] == 1);
		}
		else {
			model.addConstr(msgV[0][i] == 0);
			model.addConstr(msgD[0][i] == 0);
			model.addConstr(msgV[6][(i + 25) % 32] == 0);
			model.addConstr(msgD[6][(i + 25) % 32] == 0);
			model.addConstr(msgV[9][(i + 12) % 32] == 0);
			model.addConstr(msgD[9][(i + 12) % 32] == 0);
		}
	}


	//start from i=5 (in the second round)
	for (int i = 5; i < inputSize; i++) {
		//x[i] = x[i-4]<<<10 + (F(x[i-1],x[i-2],x[i-3]<<<10) + x[i-5]<<<10 )<<<s[i]
		//compute F(x[i-1],x[i-2],x[i-3])
		for (int j = 0; j < 32; j++) {
			//ONX(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
				//bOutV[i][j], bOutD[i][j]);
			//condition
			IFZFull(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
				bOutV[i][j], bOutD[i][j], condition[i - 1][j], condition[i - 2][j], condition[i - 3][(j + 22) % 32]);
		}

		//compute bOut + mdiff
		model.addConstr(msgCV[i][0] == 0);
		model.addConstr(msgCD[i][0] == 0);
		for (int j = 0; j < 32; j++) {
			addition(model, bOutV[i][j], bOutD[i][j],
				msgV[p[i+16-5]][j], msgD[p[i+16-5]][j],
				msgCV[i][j], msgCD[i][j],
				msgOV[i][j], msgOD[i][j],
				msgCV[i][j + 1], msgCD[i][j + 1]);
		}

		//compute msgDiff + x[i-5]<<<10
		model.addConstr(tv[i * 6][0] == 0);
		model.addConstr(td[i * 6][0] == 0);
		for (int j = 0; j < 32; j++) {
			addition(model, msgOV[i][j], msgOD[i][j],
				xv[i - 5][(j + 22) % 32], xd[i - 5][(j + 22) % 32],
				tv[i * 6][j], td[i * 6][j],
				tv[i * 6 + 1][j], td[i * 6 + 1][j],
				tv[i * 6][j + 1], td[i * 6][j + 1]);;
		}
		//rotate tv[i*6+1] to tv[i*6+2], carry is tv[i*6+3]
		rotateState(model, s[i+16-5], tv[i * 6 + 1], td[i * 6 + 1], tv[i * 6 + 2], td[i * 6 + 2],
			tv[i * 6 + 3][0], td[i * 6 + 3][0], highXV[i], highXD[i]);
		
		
		//add tv[i*6+2] + x[i-4]<<<10
		//addition
		for (int j = 0; j < 32; j++) {
			addition(model, tv[i * 6 + 2][j], td[i * 6 + 2][j],
				xv[i - 4][(j + 22) % 32], xd[i - 4][(j + 22) % 32],
				tv[i * 6 + 3][j], td[i * 6 + 3][j],
				tv[i * 6 + 4][j], td[i * 6 + 4][j],
				tv[i * 6 + 3][j + 1], td[i * 6 + 3][j + 1]);
		}
		//expansion (expand tv[i*6+4] to xv[i])
		model.addConstr(tv[i * 6 + 5][0] == 0);
		model.addConstr(td[i * 6 + 5][0] == 0);//no carry for LSB

		//if (i < 10) {
			expandState(model, 0,
				tv[i * 6 + 4], td[i * 6 + 4],
				xv[i], xd[i],
				tv[i * 6 + 5], td[i * 6 + 5]);
		//}
		/*else {
			for (int j = 0; j < 32; j++) {
				zeroAddition(model,
					tv[i * 6 + 4][j], td[i * 6 + 4][j],
					xv[i][j], xd[i][j],
					tv[i * 6 + 5][j], td[i * 6 + 5][j],
					tv[i * 6 + 5][j + 1], td[i * 6 + 5][j + 1]);
			}
		}*/
	}

	GRBLinExpr wei = 0;
	for (int j = 0; j < 5; j++) {
		for (int i = 0; i < 32; i++) {
			wei += xd[j][i];
		}
	}
	//model.addConstr(wei == 15);

	/*GRBLinExpr xwei[5];
	for (int j = 3; j < 5; j++) {
		xwei[j-3].clear();
		for (int i = 0; i < 32; i++) {
			xwei[j-3] += xd[j][i];
		}
		model.addConstr(xwei[j-3] <= 6);
	}*/

	model.setObjective(wei, GRB_MINIMIZE);

	model.optimize();

	int res = model.getObjective().getValue();

	//the last N3 steps (N3=5) [---sparse]


	//output
	string str[inputSize];
	for (int i = 0; i < inputSize; i++) {
		//getSignedString(xv[i], xd[i], condition[i], str[i]);
		getSignedString(xv[i], xd[i], str[i]);
		cout << "X" << setw(2) << i << " : " << str[i]<<" : ";
		outputSignedPos(str[i]);
	}

	for (int i = 0; i < inputSize; i++) {
		
		getSignedString(xv[i], xd[i], str[i]);
		cout << "\"" << str[i] << "\"," << endl;
	}


	string out[inputSize];
	for (int i = 5; i < inputSize; i++) {
		getSignedString(bOutV[i], bOutD[i], out[i]);
		cout << "ONX" << setw(2) << i << " : " << out[i] << endl;
	}


	string mstr;
	cout << "msg diff:" << endl;
	for (int i = 0; i < 16; i++) {
		getSignedString(msgV[i], msgD[i], mstr);
		cout << "\""<<mstr<<"\"," << endl;
	}

	return res;

	/*//check the condition
	for (int i = 5; i < 16; i++) {
		//check onx[i] with xv[i-1],xd[i-1],condition[i-1]
		if (checkONX(str[i - 1], str[i - 2], str[i - 3], out[i]) == false) {
			cout << i << endl;
		}
	}

	//check the condition
	for (int i = 16; i < inputSize; i++) {
		//check onx[i] with xv[i-1],xd[i-1],condition[i-1]
		if (checkIFZ(str[i - 1], str[i - 2], str[i - 3], out[i]) == false) {
			cout << i << endl;
		}
	}*/
}

void Primitive::model_RIPEMD160_FirstRound_Right_NonLinear_36Collision() {
	int s[32] = {
		8,9,9,11,13,15,15,5,7,7,8,11,14,14,12,6,
		9,13,15,7,12,8,9,11,7,7,12,7,6,15,13,11,
	};

	int p[32] = {
		5,14,7,0,9,2,11,4,13,6,15,8,1,10,3,12,
		6,11,3,7,0,13,5,10,14,15,8,12,4,9,1,2
	};

	string msgDiff[16] = {
		"=========n==================n===",
		"================================",
		"================================",
		"================================",
		"================================",
		"================================",
		"===u============u===============",
		"================================",
		"================================",
		"================n============n==",
		"================================",
		"================================",
		"================================",
		"================================",
		"================================",
		"================================",
	};

	const int inputSize = 30;

	string input[inputSize] = {
		"================================",//y0
		"================================",//y1
		"================================",//y2
		"=================n============n=",//y3
		"====n============n==============",//y4
		"=====================n==========",//y5
		
		"=====nuunnnnnnnnnnnnnn=un=======",//y6
		"=u=n=uun==n==nu==nnun=nuuuuuuuuu",//y7
		"n=un=nuuuu===u=un=unnnn=nn=nunuu",//y8

		"======u====n==u==uu===n====n===n",//y9
		"u===u====uu=u==========u========",//10

		"===n=========n============n=====",//y11
		"=====================n==========",//y12
		"========u========n======u===u===",//y13
		"==u======n===n===========n======",//y14
		"======u============uu======n====",//y15
		"===============n================",
		"=======u========================",
		"================================",
		"==========n=============n=======",
		"===nu===========================",
		"================================",
		"====================u===========",
		"n===========nu==================",
		"=======u==================u=====",
		"================================",//y25
		"================================",//y26
		"================================",//y27
		"================================",//y28
		"================================",//y29
	};

	GRBEnv env;
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);

	//model.set(GRB_IntParam_Aggregate, 0);

	//condition
	vector<vector<GRBVar> > condition;//x0-x15
	condition.resize(inputSize);
	initializeVar(model, condition, 32);

	//zero
	vector<GRBVar> zerov, zerod;
	zerov.resize(32);
	zerod.resize(32);
	initializeVar(model, zerov, 32);
	initializeVar(model, zerod, 32);

	//state
	vector<vector<GRBVar> > xv, xd;
	xv.resize(inputSize);
	xd.resize(inputSize);
	initializeVar(model, xv, 32);
	initializeVar(model, xd, 32);

	//intermediate state (carry, rotate)
	vector<vector<GRBVar> > tv, td;
	tv.resize(inputSize * 6);
	td.resize(inputSize * 6);
	initializeVar(model, tv, 33);
	initializeVar(model, td, 33);

	//output of boolean functions
	vector<vector<GRBVar> > bOutV, bOutD;
	bOutV.resize(inputSize);
	bOutD.resize(inputSize);
	initializeVar(model, bOutV, 32);
	initializeVar(model, bOutD, 32);

	//diff after adding the msg difference
	vector<vector<GRBVar> > msgAV(inputSize), msgAD(inputSize);
	vector<vector<GRBVar> > msgCV(inputSize), msgCD(inputSize);
	vector<vector<GRBVar> > msgOV(inputSize), msgOD(inputSize);
	initializeVar(model, msgAV, 32);
	initializeVar(model, msgAD, 32);
	initializeVar(model, msgCV, 33);
	initializeVar(model, msgCD, 33);
	initializeVar(model, msgOV, 32);
	initializeVar(model, msgOD, 32);

	//diff of the message difference
	vector<vector<GRBVar> > msgV(16), msgD(16);
	initializeVar(model, msgV, 32);
	initializeVar(model, msgD, 32);

	vector<GRBVar> highXD, highXV;
	highXD.resize(32);
	highXV.resize(32);
	initializeVar(model, highXD, 32);
	initializeVar(model, highXV, 32);

	int size = 12;
	vector<vector<GRBVar> > q(size);
	vector<vector<GRBVar> > qCarry(size);
	vector<vector<GRBVar> > outV(size), outD(size);
	vector<vector<GRBVar> > outVC(size), outDC(size);
	vector<vector<GRBVar> > innerCV(size), innerCD(size);
	vector<vector<GRBVar> > outerCV(size),outerCD(size);
	vector<vector<GRBVar> > value(size);
	initializeVar(model, q, 32);
	initializeVar(model, qCarry, 33);
	initializeVar(model, innerCV, 33);
	initializeVar(model, innerCD, 33);
	initializeVar(model, outerCV, 33);
	initializeVar(model, outerCD, 33);
	initializeVar(model, value, 32);
	initializeVar(model, outV, 32);
	initializeVar(model, outD, 32);
	initializeVar(model, outVC, 33);
	initializeVar(model, outDC, 33);

	//loading constraints
	for (int i = 0; i < inputSize; i++) {
		loadConstraint(model, input[i], xv[i], xd[i]);
	}

	for (int i = 0; i<16; i++) {
		loadConstraint(model, msgDiff[i], msgV[i], msgD[i]);
	}

	vector<GRBVar> zero(32);
	initializeVar(model, zero, 32);
	for (int i = 0; i < 32; i++) {
		model.addConstr(zero[i] == 0);
	}

	//starting from x5
	for (int i = 3; i < inputSize; i++) {
		//x[i] = x[i-4]<<<10 + (F(x[i-1],x[i-2],x[i-3]<<<10) + x[i-5]<<<10 )<<<s[i]
		//compute F(x[i-1],x[i-2],x[i-3])
		if (i < 16) {
			for (int j = 0; j < 32; j++) {
				ONX(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j]);

				//condition
				ONXCut(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j], condition[i - 1][j], condition[i - 2][j], condition[i - 3][(j + 22) % 32]);
			}
		}
		else {
			for (int j = 0; j < 32; j++) {
				IFZ(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j]);

				//condition
				IFZCut(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j], condition[i - 1][j], condition[i - 2][j], condition[i - 3][(j + 22) % 32]);
			}
		}

		//compute bOut + mdiff
		model.addConstr(msgCV[i][0] == 0);
		model.addConstr(msgCD[i][0] == 0);
		for (int j = 0; j < 32; j++) {
			addition(model, bOutV[i][j], bOutD[i][j],
				msgV[p[i]][j], msgD[p[i]][j],
				msgCV[i][j], msgCD[i][j],
				msgOV[i][j], msgOD[i][j],
				msgCV[i][j + 1], msgCD[i][j + 1]);
		}
		//compute msgDiff + x[i-5]<<<10
		if (i >= 5) {
			model.addConstr(tv[i * 6][0] == 0);
			model.addConstr(td[i * 6][0] == 0);
			for (int j = 0; j < 32; j++) {
				addition(model, msgOV[i][j], msgOD[i][j],
					xv[i - 5][(j + 22) % 32], xd[i - 5][(j + 22) % 32],
					tv[i * 6][j], td[i * 6][j],
					tv[i * 6 + 1][j], td[i * 6 + 1][j],
					tv[i * 6][j + 1], td[i * 6][j + 1]);;
			}
		}
		else {
			for (int j = 0; j < 32; j++) {
				model.addConstr(msgOV[i][j] == tv[i * 6 + 1][j]);
				model.addConstr(msgOD[i][j] == td[i * 6 + 1][j]);
			}
		}

		//rotate tv[i*4+1] to tv[i*4+2], carry is tv[i*4+3]
		rotateState(model, s[i], tv[i * 6 + 1], td[i * 6 + 1], tv[i * 6 + 2], td[i * 6 + 2],
			tv[i * 6 + 3][0], td[i * 6 + 3][0], highXV[i], highXD[i]);
		//add tv[i*4+2] + x[i-4]<<<10

		//addition
		if (i >= 4) {
			for (int j = 0; j < 32; j++) {
				addition(model, tv[i * 6 + 2][j], td[i * 6 + 2][j],
					xv[i - 4][(j + 22) % 32], xd[i - 4][(j + 22) % 32],
					tv[i * 6 + 3][j], td[i * 6 + 3][j],
					tv[i * 6 + 4][j], td[i * 6 + 4][j],
					tv[i * 6 + 3][j + 1], td[i * 6 + 3][j + 1]);
			}
		}
		else {
			for (int j = 0; j < 32; j++) {
				addition(model, tv[i * 6 + 2][j], td[i * 6 + 2][j],
					zero[j], zero[j],
					tv[i * 6 + 3][j], td[i * 6 + 3][j],
					tv[i * 6 + 4][j], td[i * 6 + 4][j],
					tv[i * 6 + 3][j + 1], td[i * 6 + 3][j + 1]);
			}
		}
		//expansion (expand tv[i*6+4] to xv[i])
		model.addConstr(tv[i * 6 + 5][0] == 0);
		model.addConstr(td[i * 6 + 5][0] == 0);//no carry for LSB

		if (i <= 10) {
			expandState(model, 0,
				tv[i * 6 + 4], td[i * 6 + 4],
				xv[i], xd[i],
				tv[i * 6 + 5], td[i * 6 + 5]);
		}
		else {
			for (int j = 0; j < 32; j++) {
				zeroAddition(model,
					tv[i * 6 + 4][j], td[i * 6 + 4][j],
					xv[i][j], xd[i][j],
					tv[i * 6 + 5][j], td[i * 6 + 5][j],
					tv[i * 6 + 5][j + 1], td[i * 6 + 5][j + 1]);
			}
		}


		//check the modular difference of y10, y9
		if (i == 10 || i==11 || i==9) {
			//first compute q[i]=xv[i] - xv[i-4]<<<10
			//convert the signed diff to conditions
			for (int j = 0; j < 32; j++) {
				model.addConstr(-xv[i][j] + condition[i][j] >= 0);
				model.addConstr(xv[i][j] - xd[i][j] - condition[i][j] + 1 >= 0);
				model.addConstr(-xv[i-4][j] + condition[i-4][j] >= 0);
				model.addConstr(xv[i-4][j] - xd[i-4][j] - condition[i-4][j] + 1 >= 0);
			}
			model.addConstr(qCarry[i][0] == 0);
			valueAddition(model, 10, condition[i - 4], q[i],condition[i], qCarry[i]);
			//second (partially) compute outVC,outVD = (tv[i * 6 + 2][j], td[i * 6 + 2][j]) + td[i*6+3][0]
			model.addConstr(outVC[i][0] == 0);
			model.addConstr(outDC[i][0] == 0);
			for (int j = 0; j <= s[i]; j++) {
				if (j == 0) {
					addition(model, tv[i * 6 + 2][j], td[i * 6 + 2][j],
						tv[i * 6 + 3][j], td[i * 6 + 3][j],
						outVC[i][j], outDC[i][j],
						outV[i][j], outD[i][j],
						outVC[i][j+1], outDC[i][j+1]);
				}
				else {
					addition(model, tv[i * 6 + 2][j], td[i * 6 + 2][j],
						zero[j], zero[j],
						outVC[i][j], outDC[i][j],
						outV[i][j], outD[i][j],
						outVC[i][j + 1], outDC[i][j + 1]);
				}
			}
			//compute value=q>>>s[i] + (tv[i * 6 + 1], td[i * 6 + 1]), until 32-s[i]
			model.addConstr(innerCV[i][0] == 0);
			model.addConstr(innerCD[i][0] == 0);
			for (int j = 0; j <= 32 - s[i]; j++) {
				model.addConstr(
					2 * (innerCD[i][j + 1] - 2 * innerCV[i][j + 1]) + value[i][j]
					== (innerCD[i][j] - 2 * innerCV[i][j]) + (td[i * 6 + 1][j] - 2 * tv[i * 6 + 1][j])
					+ q[i][(j + s[i]) % 32]
				);
				model.addConstr(innerCD[i][j + 1] >= innerCV[i][j + 1]);
			}
			//compute value<<<s[i] = q + (outVC,outVD), until s[i]
			model.addConstr(outerCV[i][0] == 0);
			model.addConstr(outerCD[i][0] == 0);
			for (int j = 0; j <= s[i]; j++) {
				model.addConstr(
					2 * (outerCD[i][j + 1] - 2 * outerCV[i][j + 1]) + value[i][(j + 32-s[i]) % 32]
					== (outerCD[i][j] - 2 * outerCV[i][j]) + (outD[i][j] - 2 * outV[i][j])
					+ q[i][j]
				);
				model.addConstr(outerCD[i][j + 1] >= outerCV[i][j + 1]);
			}
		}
	}


	GRBLinExpr wei = 0;
	for (int j = 3; j < 6; j++) {
		for (int i = 0; i < 32; i++) {
			wei += xd[j][i];
		}
	}
	//model.setObjective(wei, GRB_MINIMIZE);

	model.optimize();

	//output
	string str[inputSize];
	for (int i = 0; i < inputSize; i++) {
		getSignedString(xv[i], xd[i], str[i]);
		cout << "X" << setw(2) << i << " : " << str[i];
		outputSignedPos(str[i]);
	}

	for (int i = 0; i < inputSize; i++) {
		cout << "\"" << str[i] << "\"," << endl;
	}

	string out[inputSize];
	for (int i = 3; i < inputSize; i++) {
		getSignedString(bOutV[i], bOutD[i], out[i]);
		cout << "ONX" << setw(2) << i << " : " << out[i] << endl;
	}

	for (int i = 5; i < inputSize; i++) {
		cout << "\"" << out[i] << "\"," << endl;
	}

	for (int i = 0; i < inputSize; i++) {
		getSignedString(xv[i], xd[i], condition[i], str[i]);
	}

	//check the condition
	/*for (int i = 3; i < 16; i++) {
		//check onx[i] with xv[i-1],xd[i-1],condition[i-1]
		if (checkONX(str[i - 1], str[i - 2], str[i - 3], out[i]) == false) {
			cout << i << endl;
		}
	}

	//check the condition
	for (int i = 16; i < inputSize; i++) {
		//check onx[i] with xv[i-1],xd[i-1],condition[i-1]
		if (checkIFZ(str[i - 1], str[i - 2], str[i - 3], out[i]) == false) {
			cout << i << endl;
		}
	}*/
}

void Primitive::model_RIPEMD160_FirstRound_Right() {
	int s[32] = {
		8,9,9,11,13,15,15,5,7,7,8,11,14,14,12,6,
		9,13,15,7,12,8,9,11,7,7,12,7,6,15,13,11,
	};

	string m7 = "-n--------------n----------n----";//message difference
	string m7Inv = "-u--------------u----------u----";

	const int inputSize = 18;

	string input[inputSize] = {
		"================================",//x0
		"================================",//x1
		"=======n==========n=====n=======",//x2
		"=======n=====n==============n===",//x3
		"n==============n================",//x4
		//"????????????????????????????????",//x5
		//"????????????????????????????????",//x6
		//"????????????????????????????????",//x7
		//"????????????????????????????????",//x8
		//"????????????????????????????????",//x9
		"u==========un==nnn======uunnnnnn",
		"n=nu=nuuu=uuuu=nu===unuuuuuuu===",
		"=nuu=unuu=nun=n===n=unun==uu==n=",
		"u=uu=u==uuunn==u=nn==u==n=n===n=",
		"u===n==unnu======uu=========u==n",
		"===u=n=n=======un==n=n========n=",//x10
		"==========n=====u==============n",//x11
		"=u===n=u==============u=========",//x12
		"=========n==========n=====n=====",//x13
		"=====u=====u==============u=====",//x14
		"================================",//x15
		"================================",//x16
		"================================",//x17
	};

	string empty = "================================";

	GRBEnv env;
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);

	//condition
	vector<vector<GRBVar> > condition;//x0-x15
	condition.resize(inputSize);
	initializeVar(model, condition, 32);

	//zero
	vector<GRBVar> zerov,zerod;
	zerov.resize(32);
	zerod.resize(32);
	initializeVar(model,zerov,32);
	initializeVar(model, zerod, 32);
	loadConstraint(model, empty, zerov, zerod);

	//state
	vector<vector<GRBVar> > xv,xd;
	xv.resize(inputSize);
	xd.resize(inputSize);
	initializeVar(model, xv, 32);
	initializeVar(model, xd, 32);

	//intermediate state (carry, rotate)
	vector<vector<GRBVar> > tv, td;
	tv.resize(inputSize * 6);
	td.resize(inputSize * 6);
	initializeVar(model, tv, 33);
	initializeVar(model, td, 33);

	//output of boolean functions
	vector<vector<GRBVar> > bOutV,bOutD;
	bOutV.resize(inputSize);
	bOutD.resize(inputSize);
	initializeVar(model, bOutV, 32);
	initializeVar(model, bOutD, 32);

	vector<GRBVar> highXD, highXV;
	highXD.resize(32);
	highXV.resize(32);
	initializeVar(model, highXD, 32);
	initializeVar(model, highXV, 32);

	//loading constraints
	for (int i = 0; i < inputSize; i++) {
		loadConstraint(model, input[i], xv[i], xd[i]);
	}

	//the first N1 steps (N1=3) [---sparse]
	//load the constraints

	//the middle N2 steps (N2=5) [---dense]
	//starting from x5
	for (int i = 5; i < inputSize; i++) {
		//x[i] = x[i-4]<<<10 + (F(x[i-1],x[i-2],x[i-3]<<<10) + x[i-5]<<<10 )<<<s[i]
		//compute F(x[i-1],x[i-2],x[i-3])
		
		if (i < 16) {
			for (int j = 0; j < 32; j++) {
				ONX(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j]);

				//condition
				ONXCut(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j], condition[i - 1][j], condition[i - 2][j], condition[i - 3][(j + 22) % 32]);
			}
		}
		else if (i >= 16 && i < 32) {
			for (int j = 0; j < 32; j++) {
				//ONX(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					//bOutV[i][j], bOutD[i][j]);

				//condition
				IFZFull(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j], condition[i - 1][j], condition[i - 2][j], condition[i - 3][(j + 22) % 32]);
			}
		}



		//compute bOut + x[i-5]<<<10
		model.addConstr(tv[i*6][0] == 0);
		model.addConstr(td[i*6][0] == 0);
		for (int j = 0; j < 32; j++) {
			addition(model, bOutV[i][j], bOutD[i][j],
				xv[i-5][(j + 22) % 32], xd[i-5][(j + 22) % 32],
				tv[i*6][j], td[i*6][j],
				tv[i*6+1][j], td[i*6+1][j],
				tv[i * 6][j + 1], td[i * 6][j + 1]);;
		}
		//rotate tv[i*4+1] to tv[i*4+2], carry is tv[i*4+3]
		rotateState(model, s[i], tv[i*6+1], td[i*6+1], tv[i*6+2], td[i*6+2],
			tv[i*6+3][0], td[i*6+3][0], highXV[i], highXD[i]);
		//add tv[i*4+2] + x[i-4]<<<10
		//addition
		for (int j = 0; j < 32; j++) {
			addition(model, tv[i * 6 + 2][j], td[i * 6 + 2][j],
				xv[i - 4][(j + 22) % 32], xd[i - 4][(j + 22) % 32],
				tv[i * 6 + 3][j], td[i * 6 + 3][j],
				tv[i * 6 + 4][j], td[i * 6 + 4][j],
				tv[i * 6 + 3][j + 1], td[i * 6 + 3][j+1]);
		}
		//expansion (expand tv[i*6+4] to xv[i])
		model.addConstr(tv[i * 6 + 5][0] == 0);
		model.addConstr(td[i * 6 + 5][0] == 0);//no carry for LSB

		if (i < 10) {
			expandState(model, 0,
				tv[i * 6 + 4], td[i * 6 + 4],
				xv[i], xd[i],
				tv[i * 6 + 5], td[i * 6 + 5]);
		}
		else {
			for (int j = 0; j < 32; j++) {
				zeroAddition(model,
					tv[i * 6 + 4][j], td[i * 6 + 4][j],
					xv[i][j], xd[i][j],
					tv[i * 6 + 5][j], td[i * 6 + 5][j],
					tv[i * 6 + 5][j+1], td[i * 6 + 5][j+1]);
			}
		}
	}

	GRBLinExpr wei = 0;
	for (int j = 5; j < 10; j++) {
		for (int i = 0; i < 32; i++) {
			wei += xd[j][i];
		}
	}
	model.setObjective(wei, GRB_MINIMIZE);

	model.optimize();

	//the last N3 steps (N3=5) [---sparse]


	//output
	string str[inputSize];
	for (int i = 0; i < inputSize; i++) {
		getSignedString(xv[i], xd[i], condition[i],str[i]);
		cout << "X" << setw(2) << i << " : " << str[i] << endl;
	}

	string out[inputSize];
	for (int i = 5; i < inputSize; i++) {
		getSignedString(bOutV[i], bOutD[i], out[i]);
		cout << "ONX" <<setw(2)<<i<<" : "<< out[i] << endl;
	}

	//check the condition
	for (int i = 5; i < 16; i++) {
		//check onx[i] with xv[i-1],xd[i-1],condition[i-1]
		if (checkONX(str[i - 1], str[i - 2], str[i - 3], out[i])==false) {
			cout << i << endl;
		}
	}

	//check the condition
	for (int i = 16; i < inputSize; i++) {
		//check onx[i] with xv[i-1],xd[i-1],condition[i-1]
		if (checkIFZ(str[i - 1], str[i - 2], str[i - 3], out[i]) == false) {
			cout << i << endl;
		}
	}
}

void Primitive::model_RIPEMD160_FirstRound_Right_48Steps() {
	int s[32] = {
		9,13,15,7,12,8,9,11,7,7,12,7,6,15,13,11,
		9,7,15,11,8,6,6,14,12,13,5,14,13,13,7,5
	};

	const int inputSize = 20;

	string input[inputSize] = {
		"================================",//x0
		"================================",//x1
		"================================",//x2
		"=========================u======",//x3
		"=============u==================",//x4
		//"????????????????????????????????",//x5
		//"????????????????????????????????",//x6
		//"????????????????????????????????",//x7
		//"????????????????????????????????",//x8
		//"????????????????????????????????",//x9

		//"11011u10001unnnnnn01011000001000",
		//"nnnnu1nnn001010010010unnnnnnn001",
		//"0uunuu011nu1u1uuuuuuu0u1uun1un00",
		//"nnnuununn11u0un0uuun1u0000101nun",
		//"11unu111110100100n0uuuu00unnunuu",
		//"=====u===u=====n=====n=========u",//x10
		//"====u=====n===============n===n=",//x11
		//"===========u==1==u==1u=========0",//x12
		//"====n=====n=====1===1===========",//x13
		//"==========1==========u====1=====",//x14
		//"==========u==========1====0=====",//x15
		"00000u1unnnnnnnnnn01101010111101",
		"uun1n11nu1101010000000010111unuu",
		"01n0nu0110110unnu0nnnu1unnun0unu",
		"uuun10nuuun1u0u10100n1100u1nu0nu",
		"0110001010n1uuu010u00nu01nu0uu0n",
		"=====n===n=====n===============u",//x10
		"====u=====?===============u===u=",//x11
		"===========n==1==n==1n=========0",//x12
		"====u=====u=========0===========",//x13
		"==========1==========n==========",//x14
		"==========u==========1====-=====",//x15
		"================================",//x16
		"================================",//x17
		"================================",//x18
		"================================",//x19
	};

	GRBEnv env;
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);

	//condition
	vector<vector<GRBVar> > condition;//x0-x15
	condition.resize(inputSize);
	initializeVar(model, condition, 32);

	//zero
	vector<GRBVar> zerov, zerod;
	zerov.resize(32);
	zerod.resize(32);
	initializeVar(model, zerov, 32);
	initializeVar(model, zerod, 32);
	//loadConstraint(model, empty, zerov, zerod);

	//state
	vector<vector<GRBVar> > xv, xd;
	xv.resize(inputSize);
	xd.resize(inputSize);
	initializeVar(model, xv, 32);
	initializeVar(model, xd, 32);

	//intermediate state (carry, rotate)
	vector<vector<GRBVar> > tv, td;
	tv.resize(inputSize * 6);
	td.resize(inputSize * 6);
	initializeVar(model, tv, 33);
	initializeVar(model, td, 33);

	//output of boolean functions
	vector<vector<GRBVar> > bOutV, bOutD;
	bOutV.resize(inputSize);
	bOutD.resize(inputSize);
	initializeVar(model, bOutV, 32);
	initializeVar(model, bOutD, 32);

	vector<GRBVar> highXD, highXV;
	highXD.resize(32);
	highXV.resize(32);
	initializeVar(model, highXD, 32);
	initializeVar(model, highXV, 32);

	//loading constraints
	for (int i = 0; i < inputSize; i++) {
		loadConstraint(model, input[i], xv[i], xd[i]);
		loadCondition(model, input[i], condition[i]);
	}

	//starting from x5
	for (int i = 5; i < inputSize; i++) {
		//x[i] = x[i-4]<<<10 + (F(x[i-1],x[i-2],x[i-3]<<<10) + x[i-5]<<<10 )<<<s[i]
		//compute F(x[i-1],x[i-2],x[i-3])
		if (i < 16) {
			for (int j = 0; j < 32; j++) {
				IFZ(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j]);

				//condition and difference transition
				IFZCut(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j], condition[i - 1][j], condition[i - 2][j], condition[i - 3][(j + 22) % 32]);
			}
		}
		else if (i < 32) {
			for (int j = 0; j < 32; j++) {
				ONZ(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j]);

				//condition and difference transition
				ONZCut(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j], condition[i - 1][j], condition[i - 2][j], condition[i - 3][(j + 22) % 32]);
			}
		}

		//compute bOut + x[i-5]<<<10
		model.addConstr(tv[i * 6][0] == 0);
		model.addConstr(td[i * 6][0] == 0);
		for (int j = 0; j < 32; j++) {
			addition(model, bOutV[i][j], bOutD[i][j],
				xv[i - 5][(j + 22) % 32], xd[i - 5][(j + 22) % 32],
				tv[i * 6][j], td[i * 6][j],
				tv[i * 6 + 1][j], td[i * 6 + 1][j],
				tv[i * 6][j + 1], td[i * 6][j + 1]);;
		}
		//rotate tv[i*4+1] to tv[i*4+2], carry is tv[i*4+3]
		rotateState(model, s[i], tv[i * 6 + 1], td[i * 6 + 1], tv[i * 6 + 2], td[i * 6 + 2],
			tv[i * 6 + 3][0], td[i * 6 + 3][0], highXV[i], highXD[i]);
		//add tv[i*4+2] + x[i-4]<<<10
		//addition
		for (int j = 0; j < 32; j++) {
			addition(model, tv[i * 6 + 2][j], td[i * 6 + 2][j],
				xv[i - 4][(j + 22) % 32], xd[i - 4][(j + 22) % 32],
				tv[i * 6 + 3][j], td[i * 6 + 3][j],
				tv[i * 6 + 4][j], td[i * 6 + 4][j],
				tv[i * 6 + 3][j + 1], td[i * 6 + 3][j + 1]);
		}
		//expansion (expand tv[i*6+4] to xv[i])
		model.addConstr(tv[i * 6 + 5][0] == 0);
		model.addConstr(td[i * 6 + 5][0] == 0);//no carry for LSB

		if (i < 10) {
			expandState(model, 0,
				tv[i * 6 + 4], td[i * 6 + 4],
				xv[i], xd[i],
				tv[i * 6 + 5], td[i * 6 + 5]);
		}
		else {
			for (int j = 0; j < 32; j++) {
				zeroAddition(model,
					tv[i * 6 + 4][j], td[i * 6 + 4][j],
					xv[i][j], xd[i][j],
					tv[i * 6 + 5][j], td[i * 6 + 5][j],
					tv[i * 6 + 5][j + 1], td[i * 6 + 5][j + 1]);
			}
		}
	}

	GRBLinExpr wei = 0;
	for (int j = 11; j < inputSize; j++) {
		for (int i = 0; i < 32; i++) {
			wei += xd[j][i];
			//wei += condition[j][i];
		}
	}
	GRBLinExpr wei10 = 0;
	for (int j = 10; j < 11; j++) {
		for (int i = 0; i < 32; i++) {
			wei10 += xd[j][i];
			//wei += condition[j][i];
		}
	}
	//model.addConstr(wei10 <= 4);
	//model.addConstr(wei <= 11);
	//model.setObjective(wei10, GRB_MINIMIZE);

	model.optimize();


	//output
	string str[inputSize];
	for (int i = 0; i < inputSize; i++) {
		//getSignedString(xv[i], xd[i], condition[i], str[i]);
		getSignedString(xv[i], xd[i], str[i]);
		cout << "Y" << setw(2)   <<i+16+1<<" : " << str[i] << endl;
	}

	string out[inputSize];
	for (int i = 5; i < inputSize; i++) {
		getSignedString(bOutV[i], bOutD[i], out[i]);
		cout << "IFZ" << setw(2) << i << " : " << out[i] << endl;
	}

	//check the condition
	/*for (int i = 5; i < 16; i++) {
		if (checkIFZ(str[i - 1], str[i - 2], str[i - 3], out[i]) == false) {
			cout << i << endl;
		}
	}*/

	//get rotation constraints
	u32 inner, outer;
	string inStr, outStr;
	for (int i = 5; i < inputSize; i++) {
		outer = getValFromSignedDiff(str[i],0) - getValFromSignedDiff(str[i - 4], 10);
		inner = getValFromSignedDiff(out[i],0) + getValFromSignedDiff(str[i - 5], 10);
		getBinaryString(inner, inStr);
		getBinaryString(outer, outStr);
		cout << "step : " << i+16+1 << endl;
		cout << "shift: " << s[i] << endl;
		cout << "inner: " << inStr << endl;
		cout << "outer: " << outStr << endl;
		string st;
		getSignedString(tv[i * 6 + 1], td[i * 6 + 1], st);
		cout << st << endl;
		getSignedString(tv[i * 6 + 2], td[i * 6 + 2], st);
		cout << st << endl;
		cout << endl;
	}
}

void Primitive::autoMSGModification() {
	int s[32] = {
		8,9,9,11,13,15,15,5,7,7,8,11,14,14,12,6,
		9,13,15,7,12,8,9,11,7,7,12,7,6,15,13,11,
	};

	int p[32] = {
		5,14,7,0,9,2,11,4,13,6,15,8,1,10,3,12,
		6,11,3,7,0,13,5,10,14,15,8,12,4,9,1,2
	};

	string msgDiff[16] = {
		/*"===========n====================",//m0[+20]
		"================================",//m1
		"================================",//m2
		"================================",//m3
		"================================",//m4
		"================================",//m5
		"==================u=============",//m6[-13]
		"================================",//m7
		"================================",//m8
		"===============================n",//m9[+0]
		"================================",//m10
		"================================",//m11
		"================================",//m12
		"================================",//m13
		"================================",//m14
		"================================",//m15*/

		"=========n==================n===",//+3,+22
		"================================",
		"================================",
		"================================",
		"================================",
		"================================",
		"===u============u===============",//-28,-15
		"================================",
		"================================",
		"================n============n==",//+2,+15
		"================================",
		"================================",
		"================================",
		"================================",
		"================================",
		"================================",
	};

	const int inputSize = 30;

	string IV[5] = {
		"================================",//v0
		"================================",//v1
		"================================",//v2
		"================================",//v3
		"================================",//v4
	};

	string input[inputSize] = {

"================================",
"================================",
"================================",
"=================n============n=",
"====n============n==============",
"=====================n==========",


"=====nuunnnnnnnnnnnnnn=un=======",
"=u=n=uun==n==nu==nnun=nuuuuuuuuu",
"n=un=nuuuu===u=un=unnnn=nn=nunuu",
"======u====n==u==uu===n====n===n",
"u===u====uu=u==========u========",

"===n=========n============n=====",
"=====================n==========",
"========u========n======u===u===",
"==u======n===n===========n======",
"======u============uu======n====",
"===============n================",
"=======u========================",
"================================",
"==========n=============n=======",
"===nu===========================",
"================================",
"====================u===========",
"n===========nu==================",
"=======u==================u=====",
"================================",
"================================",
"================================",
"================================",
"================================",
	};

	string boolOut[30]={
		"================================",
		"================================",
		"================================",
		"================================",
		"=================u============u=",
"====n===========================",
"====u==u=============n==========",
"=====nn=nununuuunuuunu=uu=======",
"=n=n==n===u==unnu=un==u=nnnnnnnu",
"=n=u===n======uuu====u==uuu=nu==",
"====u==nu==uu===u=nnunu========u",
"n==nnnn=n=====nn=u=n=n=nu====uuu",
"===u====u=n=un===n===n=n==u=n===",
"=============n=======un===n====u",
"===u====u========u======n===n===",
"==u======u===u===========u======",
"====================u====n======",
"===================u=======n===n",

"================================",
"================================",
"==========n=====================",
"================================",
"================================",
"=========================n======",
"=============u==================",
"==========================u=====",
"======================n=========",
"=============================u==",
"================================",
"================================",
	};


	GRBEnv env;
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);

	vector<vector<GRBVar> > xv(inputSize), xd(inputSize);
	vector<vector<GRBVar> > con(inputSize);
	vector<vector<GRBVar> > bOutV(inputSize), bOutD(inputSize);
	vector<vector<GRBVar> > msgV(inputSize), msgD(inputSize);
	initializeVar(model, xv, 32);
	initializeVar(model, xd, 32);
	initializeVar(model, con, 32);
	initializeVar(model, bOutV, 32);
	initializeVar(model, bOutD, 32);
	initializeVar(model, msgV, 32);
	initializeVar(model, msgD, 32);

	for (int i = 0; i < inputSize; i++) {
		loadConstraint(model, input[i], xv[i], xd[i]);
		loadConstraint(model, boolOut[i], bOutV[i], bOutD[i]);
	}

	for (int i = 0; i < 16; i++) {
		loadConstraint(model, msgDiff[i], msgV[i], msgD[i]);
	}

	//convert the signed diff to conditions
	for (int i = 0; i < inputSize; i++) {
		for (int j = 0; j < 32; j++) {
			model.addConstr(-xv[i][j] + con[i][j] >= 0);
			model.addConstr(xv[i][j] - xd[i][j] - con[i][j] + 1 >= 0);
		}
	}

	for (int i = 3; i < inputSize; i++) {
		if (i < 16) {
			for (int j = 0; j < 32; j++) {
				//condition and difference transition
				ONXCut(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j], con[i - 1][j], con[i - 2][j], con[i - 3][(j + 22) % 32]);
			}
		}
		else if (i < 32) {
			for (int j = 0; j < 32; j++) {
				//condition and difference transition
				IFZCut(model, xv[i - 1][j], xd[i - 1][j], xv[i - 2][j], xd[i - 2][j], xv[i - 3][(j + 22) % 32], xd[i - 3][(j + 22) % 32],
					bOutV[i][j], bOutD[i][j], con[i - 1][j], con[i - 2][j], con[i - 3][(j + 22) % 32]);
			}
		}
	}

	u32 inner[inputSize], outer[inputSize];
	string inStr, outStr;
	for (int i = 3; i < inputSize; i++) {
		if (i >= 4) {
			outer[i] = getValFromSignedDiff(input[i], 0) - getValFromSignedDiff(input[i - 4], 10);
		}
		else {
			outer[i] = getValFromSignedDiff(input[i], 0);
		}

		if (i >= 5) {
			inner[i] = getValFromSignedDiff(boolOut[i], 0) + getValFromSignedDiff(input[i - 5], 10);
		}
		else {
			inner[i] = getValFromSignedDiff(boolOut[i], 0);
		}
		inner[i] = inner[i] + getValFromSignedDiff(msgDiff[p[i]], 0);
	}

	//satisfy the modular difference

	vector<vector<GRBVar> > tmp(inputSize * 2);
	initializeVar(model, tmp, 32);
	vector<vector<GRBVar> > q(inputSize);
	initializeVar(model, q, 32);
	vector<vector<GRBVar> > c(inputSize * 3);
	initializeVar(model, c, 33);
	vector<vector<GRBVar> > inVar(inputSize);
	initializeVar(model, inVar, 32);
	vector<vector<GRBVar>> outVar(inputSize);
	initializeVar(model, outVar, 32);

	for (int i = 5; i < inputSize; i++) {
		loadConditionFromValue(model, inner[i], inVar[i]);
		loadConditionFromValue(model, outer[i], outVar[i]);
	}

	for (int i = 0; i < inputSize; i++) {
		toConditionFromSignedDiff(model, input[i], con[i]);
	}

	int testLen = inputSize-1;
	for (int i = 5; i < testLen; i++) {
		model.addConstr(c[3 * i][0] == 0);
		model.addConstr(c[3 * i + 1][0] == 0);
		model.addConstr(c[3 * i + 2][0] == 0);
		valueAddition(model, 32 - s[i], q[i], inVar[i],
			tmp[2 * i], c[3 * i]);
		valueAddition(model, 0, q[i], outVar[i],
			tmp[2 * i + 1], c[3 * i + 1]);
		valueAddition(model, 10, con[i - 4], q[i],
			con[i], c[3 * i + 2]);

		for (int j = 0; j < 32; j++) {
			model.addConstr(tmp[2 * i][(32 - s[i] + j) % 32] == tmp[2 * i + 1][j]);
		}
	}


	model.optimize();

	for (int i = 5; i < testLen; i++) {
		getBinaryString(inner[i], inStr);
		getBinaryString(outer[i], outStr);
		cout << "step : " << i << endl;
		cout << "shift: " << s[i] << endl;
		cout << "inner: " << inStr << endl;
		cout << "outer: " << outStr << endl << endl;
	}

	if (model.get(GRB_IntAttr_Status) == 3) {
		cout << "The model is infeasible" << endl;
	}
	else {
		cout << "feasible" << endl;
		cout << "y10:";
		for (int i = 31; i >= 0; i--) {
			if (con[10][i].get(GRB_DoubleAttr_X) == 0) {
				cout << "0";
			}
			else
				cout << "1";
		}
		cout << endl;
		cout << "y6:";
		for (int i = 31; i >= 0; i--) {
			if (con[6][i].get(GRB_DoubleAttr_X) == 0) {
				cout << "0";
			}
			else
				cout << "1";
		}
		cout << endl;
		cout << "q10:";
		for (int i = 31; i >= 0; i--) {
			if (q[10][i].get(GRB_DoubleAttr_X) == 0) {
				cout << "0";
			}
			else
				cout << "1";
		}
		cout << endl;

		cout << "tmp20:";
		for (int i = 31; i >= 0; i--) {
			if (tmp[20][i].get(GRB_DoubleAttr_X) == 0) {
				cout << "0";
			}
			else
				cout << "1";
		}
		cout << endl;

		cout << "tmp21:";
		for (int i = 31; i >= 0; i--) {
			if (tmp[21][i].get(GRB_DoubleAttr_X) == 0) {
				cout << "0";
			}
			else
				cout << "1";
		}
		cout << endl;
		system("pause");
	}


	string cx[inputSize];
	for (int i = 0; i < inputSize; i++) {
		cx[i] = input[i];
	}
	for (int i = 3; i < 16; i++) {
		getConditonONX(input[i - 1], input[i - 2], input[i - 3], boolOut[i],
			cx[i - 1], cx[i - 2], cx[i - 3]);
	}
	for (int i = 16; i < inputSize; i++) {
		getConditonIFZ(input[i - 1], input[i - 2], input[i - 3], boolOut[i],
			cx[i - 1], cx[i - 2], cx[i - 3]);
	}
	for (int i = 0; i < inputSize; i++) {
		cout << setw(2) << i << ": " << cx[i] << " m"<<setw(2)<< p[i] << endl;
		//cout << "\"" << cx[i] << "\","<< endl;
	}

	cout << "Expressions:" << endl;
	for (int i = 3; i < inputSize; i++) {
		cout << "Bool:" << i << endl;
		cout << cx[i - 1] << endl;
		cout << cx[i - 2] << endl;
		//cout << cx[i - 3] << endl;
		for (int j = 0; j < 32; j++) {
			cout << cx[i - 3][(j + 10) % 32];
		}
		cout << endl;
		cout << boolOut[i] << endl << endl;
	}

	//information of inner and outer
	cout << "inner:" << endl;
	for (int i = 0; i < inputSize; i++) {
		if (i < 3) {
			cout << "0,";
		}
		else {
			cout << "0x"<<hex<<inner[i] << ", ";
		}
	}
	cout << endl;

	cout << endl << "outer:" << endl;
	for (int i = 0; i < inputSize; i++) {
		if (i < 3) {
			cout << "0,";
		}
		else {
			cout << "0x" << hex << outer[i] << ", ";
		}
	}
	cout << endl;
}

void Primitive::initializeVar(GRBModel &model, vector<vector <GRBVar> >& var, int size) {
	for (int i = 0; i < var.size(); i++) {
		var[i].resize(size);
		for (int j = 0; j < size; j++) {
			var[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
}

void Primitive::initializeVar(GRBModel& model, vector<GRBVar>& var, int size) {
	for (int i = 0; i < size; i++) {
		var[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	}
}



void Primitive::initializeVar(GRBModel& model, vector<vector <Var> >& var, int size) {
	for (int i = 0; i < var.size(); i++) {
		var[i].resize(size);
		for (int j = 0; j < size; j++) {
			var[i][j].v = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			var[i][j].d = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
}


void Primitive::initializeVar(GRBModel& model, vector<Var>& var, int size) {
	for (int i = 0; i < size; i++) {
		var[i].v = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		var[i].d = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	}
}


void Primitive::initializeVar(GRBModel& model, vector<Aux>& var) {
	for (int i = 0; i < var.size(); i++) {
		var[i].midV = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		var[i].midD = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		var[i].highV = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		var[i].highD = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	}
}

//source -> des
void Primitive::rotateState(GRBModel& model, int rot,
	vector<GRBVar>& sv, vector<GRBVar>& sd,
	vector<GRBVar>& dv, vector<GRBVar>& dd,
	GRBVar& midV, GRBVar& midD,
	GRBVar& highV, GRBVar& highD){
	//v[31:30]
	//v[31-rot,30-rot]
	for (int i = 0; i < 30 - rot; i++) {
		model.addConstr(sv[i] == dv[(i + rot) % 32]);
		model.addConstr(sd[i] == dd[(i + rot) % 32]);
	}

	rotation(model, sv[31 - rot], sd[31 - rot], sv[30 - rot], sd[30 - rot],
		dv[31], dd[31], dv[30], dd[30],
		midV, midD);

	for (int i = 32 - rot; i < 30; i++) {
		model.addConstr(sv[i] == dv[(i + rot) % 32]);
		model.addConstr(sd[i] == dd[(i + rot) % 32]);
	}
	//we do not need the carry for the MSB
	rotation(model, sv[31], sd[31], sv[30], sd[30],
		dv[(31+rot)%32], dd[(31+rot)%32], dv[(30+rot)%32], dd[(30+rot)%32], highV, highD);
}

void Primitive::expandState(GRBModel& model, int rot,
	vector<GRBVar>& sv, vector<GRBVar>& sd,
	vector<GRBVar>& ev, vector<GRBVar>& ed,
	vector<GRBVar>& cv, vector<GRBVar>& cd) {

	for (int i = 0; i < 31; i++) {
		expansion(model, sv[i], sd[i],
			cv[i], cd[i],
			ev[(i + rot) % 32], ed[(i + rot) % 32],
			cv[i + 1], cd[i + 1]);
	}
	model.addConstr(ed[(31 + rot) % 32] >= ev[(31 + rot) % 32]);
	model.addConstr(ed[(31 + rot) % 32] + sd[31] - cd[31] >= 0);
	model.addConstr(ed[(31 + rot) % 32] - sd[31] + cd[31] >= 0);
	model.addConstr(-ed[(31 + rot) % 32] + sd[31] + cd[31] >= 0);
	model.addConstr(-ed[(31 + rot) % 32] - sd[31] - cd[31] + 2 >= 0);
}

void Primitive::addState(GRBModel& model, int s0, int s1,
	vector<GRBVar>& sv0, vector<GRBVar>& sd0,
	vector<GRBVar>& sv1, vector<GRBVar>& sd1,
	vector<GRBVar>& ev, vector<GRBVar>& ed,
	vector<GRBVar>& cv, vector<GRBVar>& cd) {
	for (int j = 0; j < 32; j++) {
		addition(model, sv0[(j+32-s0)%32], sd0[(j + 32 - s0) % 32],
			sv1[(j + 32 - s1) % 32], sd1[(j + 32 - s1) % 32],
			cv[j], cd[j],
			ev[j], ed[j],
			cv[j+1], cd[j+1]);
	}
}

void Primitive::addExpandState(GRBModel& model, int s0, int s1,
	vector<GRBVar>& sv0, vector<GRBVar>& sd0,
	vector<GRBVar>& sv1, vector<GRBVar>& sd1,
	vector<GRBVar>& ev, vector<GRBVar>& ed,
	vector<GRBVar>& cv, vector<GRBVar>& cd) {
	for (int j = 0; j < 32; j++) {
		additionExpansion(model, sv0[(j + 32 - s0) % 32], sd0[(j + 32 - s0) % 32],
			sv1[(j + 32 - s1) % 32], sd1[(j + 32 - s1) % 32],
			cv[j], cd[j],
			ev[j], ed[j],
			cv[j + 1], cd[j + 1]);
	}
}

//a5 = a1<<s0 + b4<<s1
void Primitive::addExpandState(GRBModel& model, int s0, int s1,
	vector<Var>& a1, vector<Var>& b4, vector<Var>& a5, vector<Var>& c) {
	for (int j = 0; j < 32; j++) {
		additionExpansion(model, a1[(j + 32 - s0) % 32].v, a1[(j + 32 - s0) % 32].d,
			b4[(j + 32 - s1) % 32].v, b4[(j + 32 - s1) % 32].d,
			c[j].v, c[j].d,
			a5[j].v, a5[j].d,
			c[j + 1].v, c[j + 1].d);
	}
}

void Primitive::zeroState(GRBModel& model ,
	vector<GRBVar>& sv0, vector<GRBVar>& sd0,
	vector<GRBVar>& sv1, vector<GRBVar>& sd1,
	vector<GRBVar>& cv, vector<GRBVar>& cd) {
	for (int i = 0; i < 32; i++) {
		zeroAddition(model,
			sv0[i], sd0[i], sv1[i], sd1[i],
			cv[i], cd[i], cv[i + 1], cd[i + 1]);
	}
}


//0 = x - y
void Primitive::zeroState(GRBModel& model,
	vector<Var>& x, vector<Var>& y, vector<Var>& c) {
	for (int i = 0; i < 32; i++) {
		zeroAddition(model,
			x[i].v, x[i].d, y[i].v, y[i].d,
			c[i].v, c[i].d, c[i + 1].v, c[i + 1].d);
	}
}


void Primitive::rotateStateShift10(GRBModel& model, int rot,
	vector<GRBVar>& sv, vector<GRBVar>& sd,
	vector<GRBVar>& dv, vector<GRBVar>& dd,
	GRBVar& midV, GRBVar& midD,
	GRBVar& highV, GRBVar& highD) {

	for (int i = 0; i < 30 - rot; i++) {
		model.addConstr(sv[(i+22)%32] == dv[(i + rot) % 32]);
		model.addConstr(sd[(i + 22) % 32] == dd[(i + rot) % 32]);
	}
	rotation(model, sv[(31 - rot+22)%32], sd[(31 - rot+22)%32], sv[(30 - rot+22)%32], sd[(30 - rot+22)%32],
		dv[31], dd[31], dv[30], dd[30],
		midV, midD);
	for (int i = 32 - rot; i < 30; i++) {
		model.addConstr(sv[(i + 22) % 32] == dv[(i + rot) % 32]);
		model.addConstr(sd[(i + 22) % 32] == dd[(i + rot) % 32]);
	}
	//we do not need the carry for the MSB
	rotation(model, sv[(31+22)%32], sd[(31+22)%32], sv[(30+22)%32], sd[(30+22)%32],
		dv[(31 + rot) % 32], dd[(31 + rot) % 32], dv[(30 + rot) % 32], dd[(30 + rot) % 32], highV, highD);
}

void Primitive::getSignedString(vector<GRBVar>& v, vector<GRBVar>& d, string &str) {
	str.clear();
	for (int i = 31; i >= 0; i--) {
		if (d[i].get(GRB_DoubleAttr_X) == 0) {
			str = str + "=";
		}
		else {
			if (v[i].get(GRB_DoubleAttr_X) == 1) {
				str = str + "u";
			}
			else {
				str = str + "n";
			}
		}
	}
}


void Primitive::getSignedString(vector<Var>& x, string& str) {
	str.clear();
	for (int i = 31; i >= 0; i--) {
		if (x[i].d.get(GRB_DoubleAttr_X) == 0) {
			str = str + "=";
		}
		else {
			if (x[i].v.get(GRB_DoubleAttr_X) == 1) {
				str = str + "u";
			}
			else {
				str = str + "n";
			}
		}
	}
}

void Primitive::getSignedString(vector<GRBVar>& v, vector<GRBVar>& d, vector<GRBVar>& c, string& str) {
	str.clear();
	for (int i = 31; i >= 0; i--) {
		if (d[i].get(GRB_DoubleAttr_X) == 0) {
			if (c[i].get(GRB_DoubleAttr_X) == 0)
				str = str + "0";
			else
				str = str + "1";
			//str = str + "=";
		}
		else {
			if (v[i].get(GRB_DoubleAttr_X) == 1) {
				str = str + "u";
			}
			else {
				str = str + "n";
			}
		}
	}
}

void Primitive::outputSignedPos(string str) {
	for (int i = 0; i < 32; i++) {
		if (str[i] == 'n') {
			cout << "+" << 31 - i << ", ";
		}
		else if (str[i] == 'u') {
			cout << "-" << 31 - i << ", ";
		}
	}
	cout << endl;
}

////////////////////////////////
int Primitive::deduceConditions(u32 v3, u32 v2, u32 v1,GRBEnv &env) {
/*
X7: "n===unnnnnnnnnnnnnnnnnuuuuuuuuuu",
X8: "=u==n==nn===u=unu==uuu==uu=n=uu=",
X9: "==nn=un==nnn=u===nn====u===u=nn=",
X10: "==================nuu===========",
X11: "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn",
"n===nnnunuunuunnnuunuununnnunnuu",
"nn===un==unn=u===nn===uu==u=u==n"
"n=====================nuuuuuuuuu",
"===un========nuu==========uu====",
"==u=unu==nnuu=unu==uunnu==u==uu=",
"nnunnunu===nnuuuuuuuuuuunuununuu",
"uunuunnnunuuuunuuuuununnnunuuunn",
*/
	string state[11] = {
		"================================",//x0
		"================================",//x1
		/*"uuuuuuuuuuuuuuuuuunuuuuuuuuuuuu=",
		"=====u==nn===ununnuu=========n=u",
		"uuuuu=un==nn=ununn==u=uunnuu=n==",
		"===nuu=======nuuuuu=============",
		"nnnnnnnn==unnnnnnnnnnnnnnnnunnnn",*/
		"nuuuuuuuuuuuuuuuuu=nuuuuuuuuuuu=",
		"n===u=u==n=un====uuu=u=nn===u=uu",
		//"=nunxu=n==n==nnxxu==uun==nnu=un=",//has 10 free bits
		//"=xxxnux=xxxx=xx==nu==xx===x=xxxx",//has more than 10 free bits
		"=nun=u=n==n==nn==u==uun==nnu=un=",//has 10 free bits
		"====nu===========nu=============",//has more than 10 free bits
		"nnnnnnnn===unnnnnnnnnnnnnnnunnnn",//has 3 free bits
		"================================",//x12
		"================================",//x13
		"================================",//x14
		"================================",//x15
	};

	string xord[9] = {
		/*"nunuuuuuuuuuununnnnnnuununnnnun=",
		"nunnu=nn==nnu=======nnnnuunun=nn",
		"============nuuuuu==========nuu=",
		"nnnnn===unnn======unnnnnnnn==unn",
		"==nuun=====n=uuunn==un=====u==nn",
		"nnn=====n=uunuununnuuunnn===uunu",
		"nunuuuuunuununuununnnuunnnnnuu==",*/
		"nuunnununnnunununn=nnnuuunnunnn=",
		"=nnn=u=nu=n==unuu=u=u=n==nnu=u=u",
		"===============nuu============n=",
		"===un===unnn=============u=unnnn",
		"=nu=n=u===n==nu==u===nu===n=n=nu",
		"nnnuunu=n==nnnnuunnnnunnnn==uunu",
		"=nnnuuuunuunuunnunuuuunnnunuun==",
		"================================",
		"================================"
	};

	string m0 = "=========n==================n===";
	string m6 = "===u============u===============";

	char a[4];
	bool eq[4];
	bool iscondition = true;


	//GRBEnv env = GRBEnv();
	//env.set(GRB_IntParam_OutputFlag, 0);
	//env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);
	vector<vector<GRBVar> > condition;
	condition.resize(11);
	for (int i = 0; i < 11; i++) {
		condition[i].resize(32);
		for (int j = 0; j < 32; j++) {
			condition[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
	int pos[2];
	int ind[2];
	for (int i = 0; i < 7; i++) {
		//state[i+2], state[i+1], state[i]<<<10
		cout << endl<<"xor:" << i << endl;
		for (int j = 0; j < 32 ; j++) {
			a[0] = state[i + 2][31-j];
			a[1] = state[i + 1][31-j];
			a[2] = state[i][(31 - (j + 22) % 32)];
			a[3] = xord[i][31 - j];
			getConditionXOR(a,eq);
			iscondition = eq[0] | eq[1] | eq[2];
			int cnt = 0;
			if (eq[0]) {
				cout << "[" << i + 2 << "][" << j << "] ";
				ind[cnt] = i + 2;
				pos[cnt] = j;
				cnt++;
			}
			if (eq[1]) {
				cout << "[" << i + 1 << "][" << j << "] ";
				ind[cnt] = i + 1;
				pos[cnt] = j;
				cnt++;
			}
			if (eq[2]) {
				cout << "[" << i << "][" << (j+22)%32 << "] ";
				ind[cnt] = i;
				pos[cnt] = (j + 22) % 32;
				cnt++;
			}
			if (iscondition) {
				cout << "= " << eq[3] << endl;
				if(eq[3]){
					model.addConstr(condition[ind[0]][pos[0]] + condition[ind[1]][pos[1]] == 1);
				}
				else {
					model.addConstr(condition[ind[0]][pos[0]] == condition[ind[1]][pos[1]]);
				}
			}
		}
	}

	//satisfy the modular difference
	//compute q for i=0,1,2,3,4,5,6
	u32 in[9], out[9];
	in[0] = getValFromSignedDiff(m0, 0);
	out[0] = getValFromSignedDiff(state[2], 0);
	string biInString[9], biOutString[9];
	for (int i = 1; i < 4; i++) {
		in[i] = getValFromSignedDiff(xord[i - 1], 0);
		out[i] = getValFromSignedDiff(state[i + 2], 0);
	}
	for (int i = 4; i < 9; i++) {
		in[i] = getValFromSignedDiff(xord[i - 1], 0) + getValFromSignedDiff(state[i - 3], 10);
		out[i] = getValFromSignedDiff(state[i + 2], 0) - getValFromSignedDiff(state[i - 2], 10);
		if (i == 6) {
			in[i] = in[i] + getValFromSignedDiff(m6, 0);
		}
	}

	cout << "inner:" << endl;
	for (int i = 0; i < 9; i++) {
		cout << hex << "0x" << in[i] << ",";
	}
	cout << endl;
	cout << "outer:" << endl;
	for (int i = 0; i < 9; i++) {
		cout << hex << "0x" << out[i] << ",";
	}
	cout << endl;
	//system("pause");

	int s[32] = {
		11,14,15,12,5,8,7,9,11,13,14,15,6,7,9,8,
		7,6,8,13,11,9,7,15,7,12,15,9,11,7,13,12
	};

	//finding conforming internal state words
	vector<vector<GRBVar> > tmp(7 * 2);
	initializeVar(model, tmp, 32);
	vector<vector<GRBVar> > q(7);
	initializeVar(model, q, 32);
	vector<vector<GRBVar> > c(7 * 3);
	initializeVar(model, c, 33);
	vector<vector<GRBVar> > inVar(7);
	initializeVar(model, inVar, 32);
	vector<vector<GRBVar>> outVar(7);
	initializeVar(model, outVar, 32);

	for (int i = 0; i<7; i++) {
		loadConditionFromValue(model, in[i], inVar[i]);
		loadConditionFromValue(model, out[i], outVar[i]);
	}

	for (int i = 0; i < 11; i++) {
		toConditionFromSignedDiff(model, state[i], condition[i]);
	}

	//add conditions on condition[0]
	loadConditionFromValue(model, v3, condition[0]);

	//tmp[2*i] = q>>>s +in[i], tmp[2*i+1]=q + out[i], 
	//condition[i+2] = q + condition[i-2]<<<10
	//tmp[2*i]<<<s = tmp[2*i+1]
	for (int i = 2; i < 7; i++) {
		valueAddition(model, 32 - s[i], q[i], inVar[i], 
			tmp[2 * i], c[3 * i]);
		valueAddition(model, 0, q[i], outVar[i], 
			tmp[2 * i + 1], c[3 * i + 1]);
		valueAddition(model, 10, condition[i - 2], q[i],
				condition[i + 2], c[3 * i + 2]);

		for (int j = 0; j < 32; j++) {
			model.addConstr(tmp[2 * i][(32 - s[i] + j) % 32] == tmp[2 * i + 1][j]);
		}
	}

	vector<vector<GRBVar> > v(5);//v0,v1,v2,v3,v4
	initializeVar(model, v, 32);
	//compute v1,v2
	loadConditionFromValue(model, v2, v[2]);
	loadConditionFromValue(model, v1, v[1]);
	for (int i = 0; i < 7; i++) {
		valueAddition(model, 32 - s[i], q[i], inVar[i],
			tmp[2 * i], c[3 * i]);
		valueAddition(model, 0, q[i], outVar[i],
			tmp[2 * i + 1], c[3 * i + 1]);

		if (i == 0) {
			valueAddition(model, 10, v[1], q[i],
				condition[i + 2], c[3 * i + 2]);
		}
		else if (i == 1) {
			valueAddition(model, 10, v[2], q[i],
				condition[i + 2], c[3 * i + 2]);
		}

		else {
			valueAddition(model, 10, condition[i - 2], q[i],
				condition[i + 2], c[3 * i + 2]);
		}
		

		for (int j = 0; j < 32; j++) {
			model.addConstr(tmp[2 * i][(32 - s[i] + j) % 32] == tmp[2 * i + 1][j]);
		}
	}

	model.optimize();

	if (model.get(GRB_IntAttr_Status) == 3) {
		//cout << "The model is infeasible" << endl;
		return 0;
	}
	return 1;

	string output[11];
	for (int i = 0; i <9; i++) {
		output[i].clear();
		cout << "x" << i << ": ";
		for (int j = 31; j >= 0; j--) {
			if (state[i][31 - j] == '=') {
				//cout << condition[i][j].get(GRB_DoubleAttr_X);
				if (condition[i][j].get(GRB_DoubleAttr_X) == 1)
					output[i] = output[i] + "1";
				else
					output[i] = output[i] + "0";
			}
			else {
				//cout << state[i][31 - j];
				output[i] = output[i] + state[i][31 - j];
			}
		}
		cout << output[i]<<endl;
	}

	//processing output
	cout << endl;
	for (int i = 0; i < 7; i++) {
		cout << "x" << i + 2 << " : " << output[i + 2] << endl;
		cout << "x" << i + 1 << " : " << output[i + 1] << endl;
		cout << "x'" << i << ": ";
		for (int j = 31; j >= 0; j--) {
			cout << output[i][31 - (j + 22) % 32];
		}
		cout << endl;
		cout << "xr : " << xord[i] << endl << endl;
	}

	//manual check
	for (int i = 0; i < 7; i++) {
		getBinaryString(in[i], biInString[i]);
		getBinaryString(out[i], biOutString[i]);
		cout << "shift:" << s[i] << endl;
		cout << "inner:" << biInString[i] << endl;
		cout << "outer:" << biOutString[i] << endl << endl;
	}

	return 1;
}

void Primitive::getConditionXOR(char a[], bool eq[]) {
	if (a[0] == '=' && a[1] == '=' && a[2] == '=' && a[3] == '=') {
		eq[0] = 0;
		eq[1] = 0;
		eq[2] = 0;
		eq[3] = 0;
	}
	//==n ==u
	else if (a[0] == '=' && a[1] == '=' && a[2] == 'n' && a[3] == 'n') {
		eq[0] = 1;
		eq[1] = 1;
		eq[2] = 0;
		eq[3] = 0;
	}
	else if (a[0] == '=' && a[1] == '=' && a[2] == 'n' && a[3] == 'u') {
		eq[0] = 1;
		eq[1] = 1;
		eq[2] = 0;
		eq[3] = 1;
	}
	else if (a[0] == '=' && a[1] == '=' && a[2] == 'u' && a[3] == 'n') {
		eq[0] = 1;
		eq[1] = 1;
		eq[2] = 0;
		eq[3] = 1;
	}
	else if (a[0] == '=' && a[1] == '=' && a[2] == 'u' && a[3] == 'u') {
		eq[0] = 1;
		eq[1] = 1;
		eq[2] = 0;
		eq[3] = 0;
	}

	//=n= =u=
	else if (a[0] == '=' && a[1] == 'n' && a[2] == '=' && a[3] == 'n') {
		eq[0] = 1;
		eq[1] = 0;
		eq[2] = 1;
		eq[3] = 0;
	}
	else if (a[0] == '=' && a[1] == 'n' && a[2] == '=' && a[3] == 'u') {
		eq[0] = 1;
		eq[1] = 0;
		eq[2] = 1;
		eq[3] = 1;
	}
	else if (a[0] == '=' && a[1] == 'u' && a[2] == '=' && a[3] == 'n') {
		eq[0] = 1;
		eq[1] = 0;
		eq[2] = 1;
		eq[3] = 1;
	}
	else if (a[0] == '=' && a[1] == 'u' && a[2] == '=' && a[3] == 'u') {
		eq[0] = 1;
		eq[1] = 0;
		eq[2] = 1;
		eq[3] = 0;
	}


	//n== u==
	else if (a[0] == 'n' && a[1] == '=' && a[2] == '=' && a[3] == 'n') {
		eq[0] = 0;
		eq[1] = 1;
		eq[2] = 1;
		eq[3] = 0;
	}
	else if (a[0] == 'n' && a[1] == '=' && a[2] == '=' && a[3] == 'u') {
		eq[0] = 0;
		eq[1] = 1;
		eq[2] = 1;
		eq[3] = 1;
	}
	else if (a[0] == 'u' && a[1] == '=' && a[2] == '=' && a[3] == 'n') {
		eq[0] = 0;
		eq[1] = 1;
		eq[2] = 1;
		eq[3] = 1;
	}
	else if (a[0] == 'u' && a[1] == '=' && a[2] == '=' && a[3] == 'u') {
		eq[0] = 0;
		eq[1] = 1;
		eq[2] = 1;
		eq[3] = 0;
	}

	else {
		eq[0] = 0;
		eq[1] = 0;
		eq[2] = 0;
		eq[3] = 0;
	}

	//if there are fewer than 1 '='
	//there is no extra condition
}

bool Primitive::checkONX(string& x, string& y, string& z, string& w) {
	bool x0, x1, y0, y1, z0, z1, w0, w1;
	for (int i = 0; i < 32; i++) {
		getVal(x[31 - i], x0, x1);
		getVal(y[31 - i], y0, y1);
		getVal(z[31 - (i + 22) % 32], z0, z1);
		w0 = x0 ^ (y0 | (z0 ^ 1));
		w1 = x1 ^ (y1 | (z1 ^ 1));

		if (w0 == w1 && w[31 - i] == '=') {

		}
		else if (w0 == 0 && w1 == 1 && w[31 - i] == 'n') {

		}
		else if (w0 == 1 && w1 == 0 && w[31 - i] == 'u') {

		}
		else {
			cout << x[31 - i] << y[31 - i] << z[31 - (i + 22) % 32] << w[31 - i] << endl;
			cout << x0 << y0 << z0 << w0 << endl;
			cout << x1 << y1 << z1 << w1 << endl;

			cout << "ONX wrong" << endl;
			return 0;
		}
	}
	return 1;
}

bool Primitive::checkIFZ(string& x, string& y, string& z, string& w) {
	bool x0, x1, y0, y1, z0, z1, w0, w1;
	for (int i = 0; i < 32; i++) {
		getVal(x[31 - i], x0, x1);
		getVal(y[31 - i], y0, y1);
		getVal(z[31 - (i + 22) % 32], z0, z1);
		w0 = (x0 & z0) ^ (y0 & (z0 ^ 1));
		w1 = (x1 & z1) ^ (y1 & (z1 ^ 1));

		if (w0 == w1 && w[31 - i] == '=') {

		}
		else if (w0 == 0 && w1 == 1 && w[31 - i] == 'n') {

		}
		else if (w0 == 1 && w1 == 0 && w[31 - i] == 'u') {

		}
		else {
			cout << x[31 - i] << y[31 - i] << z[31 - (i + 22) % 32] << w[31 - i] << endl;
			cout << x0 << y0 << z0 << w0 << endl;
			cout << x1 << y1 << z1 << w1 << endl;

			cout << "IFZ wrong" << endl;
			return 0;
		}
	}
	return 1;
}

void Primitive::getVal(char& x, bool& val0, bool& val1) {
	if (x == 'n') {
		val0 = 0;
		val1 = 1;
	}
	else if (x == 'u') {
		val0 = 1;
		val1 = 0;
	}
	else if (x == '0') {
		val0 = val1 = 0;
	}
	else if (x == '1') {
		val0 = val1 = 1;
	}
}

void Primitive::getConditonIFX(string& x, string& y, string& z, string& w,
	string& cx, string& cy, string& cz) {
	//IFX(x,y,z<<<10) = IFZ(y,z<<<10,x)
	string tx, ty, tz;
	string tcx, tcy, tcz;
	tcx = "????????????????????????????????";
	tcy = "????????????????????????????????";
	tcz = "????????????????????????????????";
	//tx = y
	tx = y;
	//ty = z<<<10
	ty = z;
	for (int i = 0; i < 32; i++) {
		int j = (i+10) % 32;
		ty[i] = z[j];
	}
	//tz = x>>>10
	tz = x;
	for (int i = 0; i < 32; i++) {
		int j = (i - 10 + 32) % 32;
		tz[i] = x[j];
	}
	tcx=cy;
	//tcy = cz<<<10
	for (int i = 0; i < 32; i++) {
		int j = (i + 10 ) % 32;
		tcy[i] = cz[j];
	}
	//tcz = cx>>>10
	for (int i = 0; i < 32; i++) {
		int j = (i - 10 + 32) % 32;
		tcz[i] = cx[j];
	}

	getConditonIFZ(tx, ty, tz, w, tcx, tcy, tcz);

	cy = tcx;
	//cz = tcy>>>10
	for (int i = 0; i < 32; i++) {
		int j = (i - 10 + 32) % 32;
		cz[i] = tcy[j];
	}
	//cx = tcz<<<10
	for (int i = 0; i < 32; i++) {
		int j = (i + 10) % 32;
		cx[i] = tcz[j];
	}
}

void Primitive::getConditonIFZ(string& x, string& y, string& z,string& w,
	string& cx, string& cy, string& cz) {
	//w = xz + y(z+1)
	int iz = 0;
	for (int i = 0; i < 32; i++) {
		iz = (i + 10) % 32;
		if ((x[i] == 'n'|| x[i] == 'u') && y[i] == '=' && z[iz] == '=' && w[i] == '=') {
			cz[iz] = '0';
		}
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == '=' && w[i] == 'u') {
			cz[iz] = '1';
		}
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == '=' && w[i] == 'n') {
			cz[iz] = '1';
		}

		else if (x[i] == '=' && (y[i] == 'n'|| y[i] == 'u') && z[iz] == '=' && w[i] == '=') {
			cz[iz] = '1';
		}
		else if (x[i] == '=' && y[i] == 'n' && z[iz] == '=' && w[i] == 'n') {
			cz[iz] = '0';
		}
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == '=' && w[i] == 'u') {
			cz[iz] = '0';
		}
		else if (x[i] == '=' && y[i] == '=' && z[iz] == 'n' && w[i] == 'n') {
			cx[i] = '1', cy[i] = '0';
		}
		else if (x[i] == '=' && y[i] == '=' && z[iz] == 'u' && w[i] == 'u') {
			cx[i] = '0', cy[i] = '1';
		}
		//"=nn=", "=nnn",//15,16 (x=0) (x=1)
		else if (x[i] == '=' && y[i] == 'n' && z[iz] == 'n' && w[i] == '=') {
			cx[i] = '0';
		}
		else if (x[i] == '=' && y[i] == 'n' && z[iz] == 'n' && w[i] == 'n') {
			cx[i] = '1';
		}
		//"=nu=", "=nun",//17,18 (x=1) (x=0)
		else if (x[i] == '=' && y[i] == 'n' && z[iz] == 'u' && w[i] == '=') {
			cx[i] = '1';
		}
		else if (x[i] == '=' && y[i] == 'n' && z[iz] == 'u' && w[i] == 'n') {
			cx[i] = '0';
		}
		//"=un=", "=unu",//19,20 (x=1) (x=0)
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == 'n' && w[i] == '=') {
			cx[i] = '1';
		}
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == 'n' && w[i] == 'u') {
			cx[i] = '0';
		}
		//"=uu=", "=uuu",//21,22 (x=0) (x=1)
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == 'u' && w[i] == '=') {
			cx[i] = '0';
		}
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == 'u' && w[i] == 'u') {
			cx[i] = '1';
		}
		//"n=n=", "n=nn",//23,24 (y=1) (y=0)
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == 'n' && w[i] == '=') {
			cy[i] = '1';
		}
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == 'n' && w[i] == 'n') {
			cy[i] = '0';
		}
		//"n=u=", "n=un",//25,26 (y=0) (y=1)
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == 'u' && w[i] == '=') {
			cy[i] = '0';
		}
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == 'u' && w[i] == 'n') {
			cy[i] = '1';
		}
		//"u=n=", "u=nu",//27,28 (y=0) (y=1)
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == 'n' && w[i] == '=') {
			cy[i] = '0';
		}
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == 'n' && w[i] == 'u') {
			cy[i] = '1';
		}
		//"u=u=", "u=uu",//29,30 (y=1) (y=0)
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == 'u' && w[i] == '=') {
			cy[i] = '1';
		}
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == 'u' && w[i] == 'u') {
			cy[i] = '0';
		}
		//"nn=n", "uu=u",//31,32
		//"nu=n", "nu=u",//33,34 (z=1) (z=0)
		else if (x[i] == 'n' && y[i] == 'u' && z[iz] == '=' && w[i] == 'n') {
			cz[iz] = '1';
		}
		else if (x[i] == 'n' && y[i] == 'u' && z[iz] == '=' && w[i] == 'u') {
			cz[iz] = '0';
		}
		//"un=u", "un=n",//35,36 (z=1) (z=0)
		else if (x[i] == 'u' && y[i] == 'n' && z[iz] == '=' && w[i] == 'u') {
			cz[iz] = '1';
		}
		else if (x[i] == 'u' && y[i] == 'n' && z[iz] == '=' && w[i] == 'n') {
			cz[iz] = '0';
		}
	}
}

void Primitive::getConditonONX(string& x, string& y, string& z,string& w,
	string& cx, string& cy, string& cz) {
	int iz = 0;
	for (int i = 0; i < 32; i++) {
		iz = (i + 10) % 32;
		//"=nun", "=nuu",//x=0 x=1 (19,20)
		if (x[i] == '=' && y[i] == 'n' && z[iz] == 'u' && w[i] == 'n') {
			cx[i] = '0';
		}
		else if (x[i] == '=' && y[i] == 'n' && z[iz] == 'u' && w[i] == 'u') {
			cx[i] = '1';
		}
		//"=unu", "=unn",//x=0 x=1 (21,22)
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == 'n' && w[i] == 'u') {
			cx[i] = '0';
		}
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == 'n' && w[i] == 'n') {
			cx[i] = '1';
		}
		//"nn=u", "nn==",//z=0 z=1 (23,24)
		else if (x[i] == 'n' && y[i] == 'n' && z[iz] == '=' && w[i] == 'u') {
			cz[iz] = '0';
		}
		else if (x[i] == 'n' && y[i] == 'n' && z[iz] == '=' && w[i] == '=') {
			cz[iz] = '1';
		}
		//"nu=u", "nu==",//z=0 z=1 (25,26)
		else if (x[i] == 'n' && y[i] == 'u' && z[iz] == '=' && w[i] == 'u') {
			cz[iz] = '0';
		}
		else if (x[i] == 'n' && y[i] == 'u' && z[iz] == '=' && w[i] == '=') {
			cz[iz] = '1';
		}
		//"uu=n", "uu==",//z=0 z=1 (27,28)
		else if (x[i] == 'u' && y[i] == 'u' && z[iz] == '=' && w[i] == 'n') {
			cz[iz] = '0';
		}
		else if (x[i] == 'u' && y[i] == 'u' && z[iz] == '=' && w[i] == '=') {
			cz[iz] = '1';
		}
		//"un=n", "un==",//z=0 z=1 (29,30)
		else if (x[i] == 'u' && y[i] == 'n' && z[iz] == '=' && w[i] == 'n') {
			cz[iz] = '0';
		}
		else if (x[i] == 'u' && y[i] == 'n' && z[iz] == '=' && w[i] == '=') {
			cz[iz] = '1';
		}
		//"n=nu", "n=n=",//y=1 y=0 (31,32)
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == 'n' && w[i] == 'u') {
			cy[i] = '1';
		}
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == 'n' && w[i] == '=') {
			cy[i] = '0';
		}
		//"n=uu", "n=u=",//y=1 y=0 (33,34)
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == 'u' && w[i] == 'u') {
			cy[i] = '1';
		}
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == 'u' && w[i] == '=') {
			cy[i] = '0';
		}
		//"u=nn", "u=n=",//y=1 y=0 (35,36)
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == 'n' && w[i] == 'n') {
			cy[i] = '1';
		}
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == 'n' && w[i] == '=') {
			cy[i] = '0';
		}
		//"u=un", "u=u=",//y=1 y=0 (37,38)
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == 'u' && w[i] == 'n') {
			cy[i] = '1';
		}
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == 'u' && w[i] == '=') {
			cy[i] = '0';
		}

		//"==u=", (y=1)
		else if (x[i] == '=' && y[i] == '=' && z[iz] == 'u' && w[i] == '=') {
			cy[i] = '1';
		}
		//"==uu", (x=1,y=0)
		else if (x[i] == '=' && y[i] == '=' && z[iz] == 'u' && w[i] == 'u') {
			cx[i] = '1', cy[i] = '0';
		}
		//"==un", (x=0,y=0)
		else if (x[i] == '=' && y[i] == '=' && z[iz] == 'u' && w[i] == 'n') {
			cx[i] = '0', cy[i] = '0';
		}
		//"==n=", (y=1)
		else if (x[i] == '=' && y[i] == '=' && z[iz] == 'n' && w[i] == '=') {
			cy[i] = '1';
		}
		//"==nn", (x=1,y=0) 
		else if (x[i] == '=' && y[i] == '=' && z[iz] == 'n' && w[i] == 'n') {
			cx[i] = '1', cy[i] = '0';
		}
		//"==nu", (x=0,y=0)
		else if (x[i] == '=' && y[i] == '=' && z[iz] == 'n' && w[i] == 'u') {
			cx[i] = '0', cy[i] = '0';
		}
		//"=n==", (z=0)
		else if (x[i] == '=' && y[i] == 'n' && z[iz] == '=' && w[i] == '=') {
			cz[iz] = '0';
		}
		//"=n=n", (x=0,z=1)
		else if (x[i] == '=' && y[i] == 'n' && z[iz] == '=' && w[i] == 'n') {
			cx[i] = '0', cz[iz] = '1';
		}
		//"=n=u", (x=1,z=1)
		else if (x[i] == '=' && y[i] == 'n' && z[iz] == '=' && w[i] == 'u') {
		cx[i] = '1', cz[iz] = '1';
		}
		//"=u==", (z=0)
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == '=' && w[i] == '=') {
			cz[iz] = '0';
		}
		//"=u=u", (x=0,z=1)
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == '=' && w[i] == 'u') {
			cx[i] = '0', cz[iz] = '1';
		}
		//"=u=n", (x=1,z=1)
		else if (x[i] == '=' && y[i] == 'u' && z[iz] == '=' && w[i] == 'n') {
			cx[i] = '1', cz[iz] = '1';
		}
		//"n==u", 
		//"n==n", (y=0,z=1)
		else if (x[i] == 'n' && y[i] == '=' && z[iz] == '=' && w[i] == 'n') {
			cy[i] = '0', cz[iz] = '1';
		}
		//"u==n", 
		//"u==u", (y=0,z=1)
		else if (x[i] == 'u' && y[i] == '=' && z[iz] == '=' && w[i] == 'u') {
			cy[i] = '0', cz[iz] = '1';
		}
	}
}

u32 Primitive::getValFromSignedDiff(string& str,int shift) {
	u32 val = 0;
	u32 one = 1;
	for (int i = 0; i < 32; i++) {
		if (str[31 - (i+(32-shift))%32] == 'n') {
			val = val + (one << i);
		}
		else if (str[31 - (i+(32-shift))%32] == 'u') {
			val = val - (one << i);
		}
	}
	return val;
}

void Primitive::getBinaryString(u32 val, string& str) {
	str.clear();
	for (int i = 31; i >= 0; i--) {
		if (((val >> i) & 0x1) == 1) {
			str = str + "1";
		}
		else {
			str = str + "0";
		}
		if (i % 4 == 0) {
			str = str+" ";
		}
	}
}

void Primitive::loadConditionFromValue(GRBModel& model, u32 val, vector<GRBVar>& var) {
	for (int i = 0; i < 32; i++) {
		u32 v = (val >> i) & 0x1;
		if (v == 0) {
			model.addConstr(var[i] == 0);
		}
		else {
			model.addConstr(var[i] == 1);
		}
	}
}


void Primitive::valueAddition(GRBModel& model, int s0,
	vector<GRBVar>& x, vector<GRBVar>& y,
	vector<GRBVar>& z, vector<GRBVar>& c) {
	//z=x<<<s0 + y
	for (int i = 0; i < 32; i++) {
		model.addConstr(
			2 * c[i + 1] + z[i] == x[(i + 32 - s0) % 32] + y[i] + c[i]
		);
	}
}

void Primitive::toConditionFromSignedDiff(GRBModel& model, string diff, vector<GRBVar>& var) {
	for (int i = 0; i < 32; i++) {
		if (diff[31 - i] == 'n') {
			model.addConstr(var[i] == 0);
		}
		else if (diff[31 - i] == 'u') {
			model.addConstr(var[i] == 1);
		}
	}
}





	