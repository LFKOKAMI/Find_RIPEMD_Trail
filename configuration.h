#ifndef _CONFIFURATION_H_
#define _CONFIFURATION_H_

#include "gurobi_c++.h"
//#include "/opt/gurobi910/linux64/include/gurobi_c++.h"
#include <string>
#include <vector>
using namespace std;

typedef unsigned int u32;

struct Par {
	string fNa;
	int isC;
	int isF;
	int isV;
	int isK;
};

struct Var {
	GRBVar v;
	GRBVar d;
};

struct Aux {
	GRBVar midV;
	GRBVar midD;
	GRBVar highV;
	GRBVar highD;
};

const int sl[80] = { 11, 14, 15, 12, 5, 8, 7, 9, 11, 13, 14, 15, 6, 7, 9, 8,
		7, 6, 8, 13, 11, 9, 7, 15, 7, 12, 15, 9, 11, 7, 13, 12,
		11, 13, 6, 7, 14, 9, 13, 15, 14, 8, 13, 6, 5, 12, 7, 5,
		11, 12, 14, 15, 14, 15, 9, 8, 9, 14, 5, 6, 8, 6, 5, 12,
		9, 15, 5, 11, 6, 8, 13, 12, 5, 12, 13, 14, 11, 8, 5, 6 };

const int sr[80] = { 8, 9, 9, 11, 13, 15, 15, 5, 7, 7, 8, 11, 14, 14, 12, 6,
	9, 13, 15, 7, 12, 8, 9, 11, 7, 7, 12, 7, 6, 15, 13, 11,
	9, 7, 15, 11, 8, 6, 6, 14, 12, 13, 5, 14, 13, 13, 7, 5,
	15, 5, 8, 11, 14, 14, 6, 14, 6, 9, 12, 9, 12, 5, 15, 8,
	8, 5, 12, 9, 12, 5, 14, 6, 8, 13, 6, 5, 15, 13, 11, 11 };

const int pl[80] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
	7, 4, 13, 1, 10, 6, 15, 3, 12, 0, 9, 5, 2, 14, 11, 8,
	3, 10, 14, 4, 9, 15, 8, 1, 2, 7, 0, 6, 13, 11, 5, 12,
	1, 9, 11, 10, 0, 8, 12, 4, 13, 3, 7, 15, 14, 5, 6, 2,
	4, 0, 5, 9, 7, 12, 2, 10, 14, 1, 3, 8, 11, 6, 15, 13 };

const int pr[80] = { 5, 14, 7, 0, 9, 2, 11, 4, 13, 6, 15, 8, 1, 10, 3, 12,
	6, 11, 3, 7, 0, 13, 5, 10, 14, 15, 8, 12, 4, 9, 1, 2,
	15, 5, 1, 3, 7, 14, 6, 9, 11, 8, 12, 2, 10, 0, 4, 13,
	8, 6, 4, 1, 3, 11, 15, 0, 5, 12, 2, 13, 9, 7, 10, 14,
	12, 15, 10, 4, 1, 5, 8, 7, 6, 2, 13, 14, 0, 3, 9, 11 };


class Primitive {
private:
	const string fNameR[5] = { "ONX","IFZ","ONZ","IFX","XOR" };
	const string fNameL[5] = { "XOR","IFX","ONZ","IFZ","ONX" };
	int threadNum = 4;
	int outputInfo = 1;
public:
	Primitive();
	void setThreadNum(int number);
	void setOutputInfo(int info);

	void addition(
		GRBModel& model, 
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3, 
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7, 
		GRBVar& a8, GRBVar& a9
		);
	void addZero(
		GRBModel& model,
		GRBVar& a0, GRBVar& a1, int a2, int a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
		GRBVar& a8, GRBVar& a9
	);

	void zeroAddition(
		GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7
	);
	void expansion(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7
		);
	void additionExpansion(
		GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
		GRBVar& a8, GRBVar& a9
	);
	void rotation(
		GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
		GRBVar& a8, GRBVar& a9
	);
	void XOR(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7);
	void XORCut(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
		GRBVar& x, GRBVar& y, GRBVar& z);

	void ONX(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7);
	void ONXCut(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7, 
		GRBVar& x, GRBVar& y, GRBVar& z);

	void ONZ(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7);
	void ONZCut(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
		GRBVar& x, GRBVar& y, GRBVar& z);

	void IFZ(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7);
	void IFZCut(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
		GRBVar& x, GRBVar& y, GRBVar& z);
	void IFZFull(GRBModel& model,
		GRBVar& a0, GRBVar& a1, GRBVar& a2, GRBVar& a3,
		GRBVar& a4, GRBVar& a5, GRBVar& a6, GRBVar& a7,
		GRBVar& x, GRBVar& y, GRBVar& z);


	void boolFastModel(GRBModel& model, string &name, vector<Var>& x, vector<Var>& y, vector<Var>& z, vector<Var>& w);
	void boolConditionModel(GRBModel& model, string &name, vector<Var>& x, vector<Var>& y, vector<Var>& z, vector<Var>& w,
		vector<GRBVar>& xv, vector<GRBVar>& yv, vector<GRBVar>& zv);
	void modAddition(GRBModel& model, vector<Var>& x, vector<Var>& y, vector<Var>& carry, vector<Var>& z, int sx, int sy);
	void rotateFirstStrategy(GRBModel& model, int shift,
		vector<Var>& b3, vector<Var>& a1, vector<Var>& b4,
		vector<Var>& b5, vector<Var>& a5, vector<Var>& q,
		vector<Var>& carryAdd, vector<Var>& carryExp, vector<Var>& carryQ,
		Aux& aux, Par& par);
	void rotateSecondStrategy(GRBModel& model, int shift,
		vector<Var>& a1, vector<Var>& b3, vector<Var>& b4, vector<Var>& a5,
		vector<Var>& carryAdd, vector<Var>& carryExp,
		vector<Var>& q, Par& par);
	void detectContraditions(GRBModel& model, int shift,
		vector<Var>& a1, vector<Var>& a5,
		vector<Var>& b3, vector<Var>& q,
		vector<GRBVar>& a1v, vector<GRBVar>& a5v,
		vector<GRBVar>& qv, vector<GRBVar>& outCarry,
		vector<GRBVar>& value, vector<Var>& innerCarry, vector<Var>& outerCarry);

	void rotateState(GRBModel& model, int shift,
		vector<Var>& x, vector<Var>& y, Aux& aux);
	void expandState(GRBModel& model, int shift, int isK,
		vector<Var>& x, vector<Var>& y, vector<Var>& carry);
	void zeroState(GRBModel& model,
		vector<Var>& x, vector<Var>& y, vector<Var>& c);
	void addExpandState(GRBModel& model, int s0, int s1,
		vector<Var>& a1, vector<Var>& b4, vector<Var>& a5, vector<Var>& c);

	//initialize variables
	void initializeVar(GRBModel &model, vector<vector <GRBVar> >& var, int size);
	void initializeVar(GRBModel& model, vector<vector <Var> >& var, int size);
	void initializeVar(GRBModel& model, vector<GRBVar>& var);
	void initializeVar(GRBModel& model, vector<Var>& var);
	void initializeVar(GRBModel& model, vector<Aux>& var);

	void setZero(GRBModel& model, Var& var);
	void loadConstraint(GRBModel& model, string& diff,vector<Var>& x);//load constraints
	void loadCondition(GRBModel& model, string diff,vector<GRBVar>& v);//load constraints

	void modelRoundFunction(
		GRBModel& model, int step,
		Par& par, int shift, 
		vector<Var>& m, vector<Var>& a0, vector<Var>& a1, vector<Var>& a2, vector<Var>& a3, vector<Var>& a4, vector<Var>& a5,
		vector<Var>& b1, vector<Var>& b2, vector<Var>& b3, vector<Var>& b4, vector<Var>& b5,
		vector<Var>& c1, vector<Var>& c2, vector<Var>& c3, vector<Var>& c4, vector<Var>& c5,
		vector<GRBVar>& v1, vector<GRBVar>& v2, vector<GRBVar>& v3, vector<GRBVar>& v4, vector<GRBVar>& v5,
		vector<GRBVar>& qv, vector<GRBVar>& outCarry, vector<GRBVar>& value, vector<Var>& innerCarry, vector<Var>& outerCarry,
		vector<Var>& q, Aux& aux
		);

	bool buildModel(int start, int end, int opStart, int opSteps, int isRight,
		vector<int>& isC, vector<int>& isF, vector<int>& isV,
		vector<string>& diff,
		string mPattern[],
		vector<string>& bfDiff);

	//check the result of the solver by testing all equations on q
	bool autoCheck(int start, int isRight, 
		string msgDiff[], vector<string>& input, vector<string>& boolOut);
	void loadConditionFromValue(GRBModel& model, u32 val, vector<GRBVar>& var);
	void valueAddition(GRBModel& model, int s0,
		vector<GRBVar>& x, vector<GRBVar>& y,
		vector<GRBVar>& z, vector<GRBVar>& c);
	void toConditionFromSignedDiff(GRBModel& model, string diff, vector<GRBVar>& var);

	/*
	//some useful functions
	*/

	//deduce conditions
	void outputXORCondition(int step, string x, string y, string z, string w);
	void getConditionXOR(char a[], bool eq[]);
	void getConditonIFZ(string& x, string& y, string& z, string& w,
		string& cx, string& cy, string& cz);
	void getConditonIFX(string& x, string& y, string& z, string& w,
		string& cx, string& cy, string& cz);
	void getConditonONX(string& x, string& y, string& z, string& w,
		string& cx, string& cy, string& cz);
	void getConditonONZ(string& x, string& y, string& z, string& w,
		string& cx, string& cy, string& cz);

	//compute (i-shift+32)%32
	int minus(int shift);

	void getSignedString(vector<Var>& x, string& str);

	//get value
	u32 getValFromSignedDiff(string &str,int shift=0);
	void getBinaryString(u32 val,string &str);
};

#endif
