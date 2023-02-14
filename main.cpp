#include "configuration.h"
#include "gurobi_c++.h"
#include <random>
#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
using namespace std;
//#include "/opt/gurobi910/linux64/include/gurobi_c++.h"

void testProbability_36Collision_Left() {
	Primitive ripemd;

	//left(36 steps)
	//ripemd.model_RIPEMD160_FirstRound_Left();
	int sh[2];
	int length = 2;
	while (0) {
		cout << "length:";
		cin >> length;
		cout << "bit position:";
		for (int i = 0; i < length; i++) {
			cin >> sh[i];
		}
		ripemd.model_RIPEMD160_FirstRound_Left_36Collision(sh, length);
	}


	std::mt19937 mt;            // メルセンヌ・ツイスタの32ビット版
	std::random_device rnd;     // 非決定的な乱数生成器
	mt.seed(rnd());

	int count[2] = { 0,0 };
	int testTimes = 1 << 0;

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);

	for (int i = 0; i < testTimes; i++) {
		u32 v3 = mt();
		u32 v2 = mt();
		u32 v1 = mt();
		if (ripemd.deduceConditions(v3, v2, v1, env) == 0) {
			count[0]++;
		}
		else {
			count[1]++;
		}
		if (i % (1 << 12) == 0) {
			cout << hex << (i / (1 << 12)) << endl;
			cout << "count[0]:" << hex << count[0] << endl;
			cout << "count[1]:" << hex << count[1] << endl << endl;
		}
	}
	cout << hex << count[0] << endl;
	cout << hex << count[1] << endl;
}

void outputCoefficientMatrices() {
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

	int row = 19;
	int col = 9;
	int index = 1;
	cout << "\\[" << endl;
	cout << "\\mathcal{H}_"<<index<<" | C_"<<index<<"^T =" << endl;
	cout << "\\left[{\\begin{array}{";
	for (int i = 0; i < col - 1; i++)
		cout << "c";
	cout << "|c";
	cout << "}" << endl;

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (j != col - 1) {
				cout << ma[i][j] << "& ";
			}
			else
				cout << 1 - ma[i][j] << "\\\\" << endl;
		}
	}
	cout << "\\end{array}} \\right]," << endl;
	cout << "\\]" << endl;
}

int main0() {
	outputCoefficientMatrices();
	system("pause");

	Primitive ripemd;

	/*GRBEnv env;
	//env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);
	int length = 2;
	int sh[2];
	int minIndex[2] = {0}, min = 10000000;
	int store[32 * 31] = {0};
	int cnt = 0;
	for (int i = 3; i < 4; i++) {
		for (int j = 22; j < 23; j++) {
			sh[0] = i;
			sh[1] = j;
			//sh[0] = 20;
			int res = ripemd.model_RIPEMD160_FirstRound_Right_36Collision(sh, length, env);
			if (res < min) {
				min = res;
				minIndex[0] = i;
				minIndex[1] = j;
			}
			cout << "current min:" << minIndex[0] << " " << minIndex[1] << " : " << min << endl;
			store[cnt] = res;
			cnt++;
			cout << "current:";
			for (int k = 0; k < cnt; k++) {
				cout << store[k] << ",";
			}
			cout << endl << endl;
		}
	}*/
	//ripemd.model_RIPEMD160_FirstRound_Right_NonLinear_36Collision();
	ripemd.autoMSGModification();
	//testProbability_36Collision_Left();

	//right(36 steps)
	//ripemd.model_RIPEMD160_FirstRound_Right();

	//right(48 steps)
	//ripemd.model_RIPEMD160_FirstRound_Right_48Steps();
	return 0;
}


int rightBranch() {
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

	const int varLength = 22;
	const int size = varLength - 5;

	int start = 6, end = 10;
	int isRight = 1;
	int isC[size] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
	int isF[size] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
	int isV[size] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };

	string diff[varLength] = {
		"================================",//y1
		"================================",//y2
		"=================n============n=",//y3
		"====n============n==============",//y4
		"=====================n==========",//y5

		//"????????????????????????????????",//y6
		//"????????????????????????????????",//y7
		//"????????????????????????????????",//y8
		//"????????????????????????????????",//y9
		//"????????????????????????????????",//y10
		"10001nuunnnnnnnnnnnnnn=un1101110",//y6
		"0u0n1uun00n10nu01nnun=nuuuuuuuuu",//y7
		"n1un0nuuuu1=0u0un0unnnn1nn0nunuu",//y8
		"=1=010u1000n00u01uu010n101=n100n",//y9
		"u1=0u0110uu=u011=0=1=0=u1=1=0111",//10

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
	};

	vector<string> input(varLength);
	for (int i = 0; i < varLength; i++) {
		input[i] = diff[i];
	}
	vector<string> bfDiff(size);//signed diff of the output of the boolean function
	

	for (int i = 0; i < 16; i++)
		cout << msgDiff[i] << endl;

	cout << endl;

	for (int i = 0; i < varLength; i++)
		cout << input[i] << endl;

	cout << endl;
	system("pause");

	bool find = 1;
	Primitive primitive;
	find=primitive.buildModel(start, end, isRight,
		isC, isF, isV,
		input,msgDiff,bfDiff);
	
	if (!find) {
		cout << "wrong diff, please adjust inDiff/outDiff" << endl;
		return 1;
	}

	cout << "signed difference of states:" << endl;
	for (int i = 0; i < input.size(); i++)
		cout << "Y" << start - 5 + i << ": " << input[i] << endl;
	cout << endl;

	cout << "signed difference of boolean functions:" << endl;
	for (int i = 0; i < size; i++) {
		cout << "Diff of BoolFun Step " << i+start << ": " << bfDiff[i] << endl;
	}
	cout << endl;


	//output other information
	find = primitive.autoCheck(start, isRight, msgDiff, input, bfDiff);
	if (!find)
		cout << "wrong diff, please adjust isF/isV/isC or inDiff/outDiff" << endl;

	return 1;
}




int leftBranch() {
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


	const int varLength = 15;
	const int size = varLength - 5;

	int start = 0, end = 4;
	int isRight = 0;//left branch
	int isC[10] = { 0,0,0,0,0,0,0,0,0,0 };
	int isF[10] = { 1,1,1,1,0,1,1,1,1,1 };
	int isV[10] = { 0,0,0,0,0,0,0,0,0,0 };

	string diff[varLength] = {
		"================================",//x-5
		"================================",//x-4
		"================================",//x-3
		"================================",//x-2
		"================================",//x-1

		//"????????????????????????????????",//x0
		//"????????????????????????????????",//x1
		//"????????????????????????????????",//x2
		//"????????????????????????????????",//x3
		//"????????????????????????????????",//x4
		"nuuuuuuuuuuuuuuuuu=nuuuuuuuuuuu=",//x0
		"n===u=u==n=un====uuu=u=nn===u=uu",//x1
		"=nun=u=n==n==nn==u==uun==nnu=un=",//x2
		"====nu===========nu=============",//x3
		"nnnnnnnn===unnnnnnnnnnnnnnnunnnn",//x4

		"================================",//x5
		"================================",//x6
		"================================",//x7
		"================================",//x8
		"================================",//x9
	};

	vector<string> input(varLength);
	for (int i = 0; i < varLength; i++) {
		input[i] = diff[i];
	}
	vector<string> bfDiff(size);//signed diff of the output of the boolean function


	bool find = 1;
	Primitive primitive;
	find = primitive.buildModel(start, end, isRight,
		isC, isF, isV,
		input, msgDiff, bfDiff);

	if (!find) {
		cout << "wrong diff, please adjust inDiff/outDiff" << endl;
		return 1;
	}

	cout << "signed difference of states:" << endl;
	for (int i = 0; i < input.size(); i++)
		cout << "Y" << start - 5 + i << ": " << input[i] << endl;
	cout << endl;

	cout << "signed difference of boolean functions:" << endl;
	for (int i = 0; i < size; i++) {
		cout << "Diff of BoolFun Step " << i + start << ": " << bfDiff[i] << endl;
	}
	cout << endl;


	//output other information
	find = primitive.autoCheck(start, isRight, msgDiff, input, bfDiff);
	if (!find)
		cout << "wrong diff, please adjust isF/isV/isC or inDiff/outDiff" << endl;

	//....................

	return 1;
}


//search the sparse part
int linearPropagation() {

	return 1;
}

int search(int start, int end, int opStart, int opSteps, int isRight, 
	vector<int>& isC, vector<int>& isF, vector<int>& isV,
	vector<string>& diff, string msgDiff[], vector<string>& bfDiff) {

	bool find = 1;
	Primitive primitive;
	find = primitive.buildModel(start, end, opStart, opSteps, isRight,
		isC, isF, isV,
		diff, msgDiff, bfDiff);

	if (!find) {
		cout << "wrong diff, please adjust inDiff/outDiff" << endl;
		return 1;
	}

	cout << "signed difference of states:" << endl;
	for (int i = 0; i < diff.size(); i++)
		cout << "Y" << start - 5 + i << ": " << diff[i] << endl;
	cout << endl;

	cout << "signed difference of boolean functions:" << endl;
	for (int i = 0; i < diff.size()-5; i++) {
		cout << "Diff of BoolFun Step " << i + start << ": " << bfDiff[i] << endl;
	}
	cout << endl;


	//output other information
	find = primitive.autoCheck(start, isRight, msgDiff, diff, bfDiff);
	if (!find)
		cout << "wrong diff, please adjust isF/isV/isC or inDiff/outDiff" << endl;
}

void outputParameters(int varLength, int start, int end, int opStart, int opSteps, int isRight, vector<int>& isC, vector<int>&  isF, vector<int>& isV) {
	int size = varLength - 5;
	cout << "length of steps:" << varLength << endl;
	cout << "start = " << start << " , end=" << end << endl;
	if (opSteps > 0) {
		cout << "optimize step ";
		for (int i = 0; i < opSteps; i++) {
			cout << i + opStart << " ";
		}
		cout << endl;
	}

	cout << "isRight:" << isRight << endl;
	cout << "pos: ";
	for (int i = 0; i < size; i++) {
		cout << setw(3) << i + start << " ";
	}
	cout << endl;
	cout << "isC: ";
	for (int i = 0; i < size; i++) {
		cout << setw(3) << isC[i] << " ";
	}
	cout << endl;
	cout << "isF: ";
	for (int i = 0; i < size; i++) {
		cout << setw(3) << isF[i] << " ";
	}
	cout << endl;
	cout << "isV: ";
	for (int i = 0; i < size; i++) {
		cout << setw(3) << isV[i] << " ";
	}
	cout << endl;
}

int searchTrail_readParameters(string msgFile, string diffFile, string parFile,
	int& varLength, int& start, int& end, int& opStart, int& opSteps, int& isRight,
	vector<int>& isC, vector<int>& isF, vector<int>& isV,
	string msgDiff[], vector<string> &diff, vector<string>& bfDiff) {
	
	ifstream file(msgFile);
	string str;
	if (file.is_open()) {
		for (int i = 0; i < 16; i++) {
			getline(file, str, ' ');
			getline(file, msgDiff[i]);
			//cout << str << ": " << msgDiff[i] << endl;
		}
	}
	else
		cout << "opening "<<msgFile <<" failed!" << endl;
	file.close();

	int size = 0;
	file.open(parFile);
	if(file.is_open()) {
		file >> varLength;
		file >> start;
		file >> end;
		file >> opStart;
		file >> opSteps;
		file >> isRight;
		size = varLength - 5;
		isC.resize(size);
		isF.resize(size);
		isV.resize(size);
		for (int i = 0; i < size; i++)
			file >> isC[i];
		for (int i = 0; i < size; i++)
			file >> isF[i];
		for (int i = 0; i < size; i++)
			file >> isV[i];
	}
	else
		cout << "opening " << parFile << " failed!" << endl;
	file.close();

	file.open(diffFile);
	diff.resize(varLength);
	if (file.is_open()) {
		for (int i = 0; i < varLength; i++) {
			getline(file, str, ' ');
			getline(file, diff[i]);
			//cout << str << ": " << diff[i] << endl;
		}
	}
	else
		cout << "opening " << diffFile << " failed!" << endl;

	bfDiff.resize(size);//signed diff of the output of the boolean function

	return 1;
}

int searchTrail(string msgFile, string diffFile, string parFile) {
	string msgDiff[16];

	int varLength = 0, start = 0, end = 0, opStart=0, opSteps =0, isRight = 0;
	vector<int> isC, isF, isV;
	isC.clear();
	isF.clear();
	isV.clear();

	vector<string> diff,bfDiff;
	diff.clear();
	bfDiff.clear();

	searchTrail_readParameters(msgFile, diffFile, parFile,
		varLength, start, end, opStart, opSteps, isRight,
		isC, isF, isV,
		msgDiff, diff, bfDiff);

	outputParameters(varLength, start, end, opStart, opSteps, isRight, isC, isF, isV);

	search(start, end, opStart, opSteps, isRight, isC, isF, isV, diff, msgDiff, bfDiff);

	isC.clear();
	isF.clear();
	isV.clear();
	diff.clear();
	bfDiff.clear();
}

int main() {
	//searchTrail("data/left/test_msgFile.txt", "data/left/test_diffFile.txt", "data/left/test_parFile.txt");
	//searchTrail("data/right/test_msgFile.txt", "data/right/test_diffFile.txt", "data/right/test_parFile.txt");
	
	//searchTrail("data/right/msgFile.txt", "data/right/diffFile.txt", "data/right/parFile.txt");
	searchTrail("data/left/msgFile.txt", "data/left/diffFile.txt", "data/left/parFile.txt");

	//searchTrail("data/right/sparse_msgFile.txt", "data/right/sparse_diffFile.txt", "data/right/sparse_parFile.txt");
	//searchTrail("data/left/sparse_msgFile.txt", "data/left/sparse_diffFile.txt", "data/left/sparse_parFile.txt");
	return 1;
}