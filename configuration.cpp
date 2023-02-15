#include "configuration.h"
#include "gurobi_c++.h"
//#include "/opt/gurobi910/linux64/include/gurobi_c++.h"
#include <string>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

Primitive::Primitive() {
	threadNum = 4;
	outputInfo = 1;
}

void Primitive::setThreadNum(int number) {
	threadNum = number;
}

void Primitive::setOutputInfo(int info) {
	outputInfo = info;
}

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
	else {//the second way to do expansion ... 
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
		q.resize(32);
		carryQ.resize(33);
		initializeVar(model, q);
		initializeVar(model, carryQ);

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
	vector<Var>& carryAdd, vector<Var>& carryExp,
	vector<Var>& q, Par& par) {
	//b4 = expand(b3)
	setZero(model, carryAdd[0]);
	zeroState(model, b3, b4, carryAdd);
	//a5 = a1<<<10 + b4<<<shift
	setZero(model, carryExp[0]);
	addExpandState(model, 10, shift, a1, b4, a5, carryExp);

	if (par.isV) {
		q.resize(32);
		initializeVar(model, q);
		for (int i = 0; i < shift + 1; i++) {
			model.addConstr(q[i].v == b4[(i - shift + 32) % 32].v);
			model.addConstr(q[i].d == b4[(i - shift + 32) % 32].d);
		}
	}
}

//yv = x<<<10 + qv
void Primitive::detectContraditions(GRBModel& model, int shift,
	vector<Var>& a1, vector<Var>& a5,
	vector<Var>& b3, vector<Var>& q,
	vector<GRBVar>& a1v, vector<GRBVar>& a5v, 
	vector<GRBVar>& qv,vector<GRBVar>& outCarry,
	vector<GRBVar>& value, vector<Var>& innerCarry, vector<Var>& outerCarry) {

	qv.resize(32);
	value.resize(32);
	outCarry.resize(33);
	innerCarry.resize(32 - shift + 2);
	outerCarry.resize(shift + 2);

	initializeVar(model, qv);
	initializeVar(model, value);
	initializeVar(model, outCarry);
	initializeVar(model, innerCarry);
	initializeVar(model, outerCarry);

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

	setZero(model, c1[0]);
	setZero(model, c2[0]);

	modAddition(model, m, b1, c1, b2, 0, 0);//b2 = b1+m
	modAddition(model, b2, a0, c2, b3, 0, 10);//b3 = a0<<<10 + b2

	
	//b4 = b3<<< sr + a1 [shift]
	if (par.isF == 1) {//use the first strategy
		rotateFirstStrategy(model, shift, b3, a1, b4, b5, a5, q, c3, c4, c5, aux, par);
	}
	else {//use the second strategy
		rotateSecondStrategy(model, shift, a1, b3, b4, a5, c3, c4, q, par);
	}

	if (par.isV) {//detect more contraditions
		detectContraditions(model, shift, a1, a5, b3, q, v1, v5, qv, outCarry, value, innerCarry, outerCarry);
	}
	
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

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, outputInfo);
	env.set(GRB_IntParam_Threads, threadNum);
	GRBModel model = GRBModel(env);

	initializeVar(model, a, 32);
	initializeVar(model, av, 32);
	initializeVar(model, m, 32);
	initializeVar(model, b, 32);
	//initializeVar(model, q, 32);

	initializeVar(model, cMod1, 33);
	initializeVar(model, cMod2, 33);
	initializeVar(model, cMod3, 33);
	initializeVar(model, cExp, 33);
	//initializeVar(model, cQ, 33);

	initializeVar(model, aux);

	for (int i = 0; i < 16; i++) {
		loadConstraint(model, mPattern[i], m[i]);
	}

	for (int i = 0; i < varLen; i++) {
		loadConstraint(model, diff[i], a[i]);
	}

	for (int i = 0; i < varLen; i++) {
		loadCondition(model, diff[i], av[i]);
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
	if (opSteps > 0) {
		GRBLinExpr hwSum = 0;
		for (int i = 0; i < opSteps; i++) {
			int k = i + (opStart - start) + 5;
			for (int j = 0; j < 32; j++) {
				hwSum = hwSum + a[k][j].d;
			}
		}
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
	}
	for (int i = 0; i < size; i++) {
		getSignedString(b[5 * i], bfDiff[i]);
	}
	return 1;
}


bool Primitive::autoCheck(int start, int isRight, 
	string msgDiff[], vector<string> &input, vector<string>& boolOut) {
	GRBEnv env;
	env.set(GRB_IntParam_OutputFlag, outputInfo);
	env.set(GRB_IntParam_Threads, threadNum);
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
			getConditonONZ(input[i - 1], input[i - 2], input[i - 3], boolOut[u],
				cx[i - 1], cx[i - 2], cx[i - 3]);
		}
		else if (name=="IFX") {
			getConditonIFX(input[i - 1], input[i - 2], input[i - 3], boolOut[u],
				cx[i - 1], cx[i - 2], cx[i - 3]);
		}
	}

	cout << "after adding bit conditions:" << endl;
	for (int i = 0; i < inputSize; i++) {
		int k = i - 5 + start;
		cout << setw(2) <<dec<< k << " " << cx[i] << endl;
	}

	//output the conditions for XOR
	for (int i = 5; i < inputSize; i++) {
		int k = i - 5 + start;
		int u = i - 5;
		string name = fNameL[k / 16];
		if (isRight == 1)
			name = fNameR[k / 16];

		if (name == "XOR") {
			outputXORCondition(k, input[i-1], input[i - 2], input[i - 3], boolOut[u]);
		}
	}

	return 1;
}














void Primitive::initializeVar(GRBModel &model, vector<vector <GRBVar> >& var, int size) {
	for (int i = 0; i < var.size(); i++) {
		var[i].resize(size);
		for (int j = 0; j < size; j++) {
			var[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
}

void Primitive::initializeVar(GRBModel& model, vector<GRBVar>& var) {
	for (int i = 0; i < var.size(); i++) {
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


void Primitive::initializeVar(GRBModel& model, vector<Var>& var) {
	for (int i = 0; i < var.size(); i++) {
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


//0 = x - y
void Primitive::zeroState(GRBModel& model,
	vector<Var>& x, vector<Var>& y, vector<Var>& c) {
	for (int i = 0; i < 32; i++) {
		zeroAddition(model,
			x[i].v, x[i].d, y[i].v, y[i].d,
			c[i].v, c[i].d, c[i + 1].v, c[i + 1].d);
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


void Primitive::outputXORCondition(int step, string x, string y, string z, string w) {
	int pos[2];
	int ind[2];
	char a[4];
	bool eq[4];
	bool isCondition;
		
	for (int j = 0; j < 32; j++) {
		a[0] = x[31 - j];
		a[1] = y[31 - j];
		a[2] = z[(31 - ((j + 22) % 32))];
		a[3] = w[31 - j];
		getConditionXOR(a, eq);
		isCondition = eq[0] | eq[1] | eq[2];
		int cnt = 0;
		if (eq[0]) {
			//cout << "[" << step-1 << "][" << j << "] ";
			ind[cnt] = step - 1;
			pos[cnt] = j;
			cnt++;
		}
		if (eq[1]) {
			//cout << "[" << step - 2 << "][" << j << "] ";
			ind[cnt] = step - 2;
			pos[cnt] = j;
			cnt++;
		}
		if (eq[2]) {
			//cout << "[" << step - 3 << "][" << (j + 22) % 32 << "] ";
			ind[cnt] = step - 3;
			pos[cnt] = (j + 22) % 32;
			cnt++;
		}
		if (isCondition) {
			//cout << "= " << eq[3] << endl;
			cout << "[" << ind[cnt - 2] << "][" << pos[cnt - 2] << "] + ";
			cout << "[" << ind[cnt - 1] << "][" << pos[cnt - 1] << "] = ";
			cout << eq[3] << endl;
		}
	}
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

void Primitive::getConditonONZ(string& x, string& y, string& z, string& w,
	string& cx, string& cy, string& cz) {
	//ONZ(x,y,z<<<10) = ONX(z<<<10, x, y)
	string tx, ty, tz;
	string tcx, tcy, tcz;
	tcx = "????????????????????????????????";
	tcy = "????????????????????????????????";
	tcz = "????????????????????????????????";
	
	//tx = z<<<10
	tx = z;
	for (int i = 0; i < 32; i++) {
		int j = (i + 10) % 32;
		tx[i] = z[j];
	}

	//ty = x
	ty = x;

	//tz = y>>>10
	tz = y;
	for (int i = 0; i < 32; i++) {
		int j = (i - 10 + 32) % 32;
		tz[i] = y[j];
	}

	//tcx = cz<<<10
	for (int i = 0; i < 32; i++) {
		int j = (i + 10) % 32;
		tcx[i] = cz[j];
	}

	tcy = cx;

	//tcz = cy>>>10
	for (int i = 0; i < 32; i++) {
		int j = (i - 10 + 32) % 32;
		tcz[i] = cy[j];
	}

	getConditonONX(tx, ty, tz, w, tcx, tcy, tcz);

	cx = tcy;
	//cy = tcz<<<10
	for (int i = 0; i < 32; i++) {
		int j = (i + 10) % 32;
		cy[i] = tcz[j];
	}
	//cz = tcx>>>10
	for (int i = 0; i < 32; i++) {
		int j = (i - 10 + 32) % 32;
		cz[i] = tcx[j];
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





	