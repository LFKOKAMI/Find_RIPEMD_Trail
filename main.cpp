#include "configuration.h"
#include "gurobi_c++.h"
#include <random>
#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
using namespace std;
//#include "/opt/gurobi910/linux64/include/gurobi_c++.h"

int search(int threadNum, int outputInfo, int start, int end, int opStart, int opSteps, int isRight, 
	vector<int>& isC, vector<int>& isF, vector<int>& isV,
	vector<string>& diff, string msgDiff[], vector<string>& bfDiff) {

	bool find = 1;
	Primitive primitive;

	primitive.setThreadNum(threadNum);
	primitive.setOutputInfo(outputInfo);

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

int searchTrail(int threadNum, int outputInfo, string msgFile, string diffFile, string parFile) {
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

	search(threadNum, outputInfo, start, end, opStart, opSteps, isRight, isC, isF, isV, diff, msgDiff, bfDiff);

	isC.clear();
	isF.clear();
	isV.clear();
	diff.clear();
	bfDiff.clear();
}

int main() {
	int threadNum = 4;
	int outputInfo = 1;

	cout << "Please input the number of used threads:";
	cin >> threadNum;

	cout << "Do you allow the solver to output the solving process? (0 / 1) - (allowed / not allowed):";
	cin >> outputInfo;
	if (outputInfo != 0)
		outputInfo = 1;

	//check the correctness of the trail
	//searchTrail(threadNum, outputInfo, "data/left/test_msgFile.txt", "data/left/test_diffFile.txt", "data/left/test_parFile.txt");
	//searchTrail(threadNum, outputInfo, "data/right/test_msgFile.txt", "data/right/test_diffFile.txt", "data/right/test_parFile.txt");
	
	//search the solution of unknown parts
	//searchTrail(threadNum, outputInfo, "data/right/msgFile.txt", "data/right/diffFile.txt", "data/right/parFile.txt");
	//searchTrail(threadNum, outputInfo, "data/left/msgFile.txt", "data/left/diffFile.txt", "data/left/parFile.txt");

	//search sparse parts
	searchTrail(threadNum, outputInfo, "data/right/sparse_msgFile.txt", "data/right/sparse_diffFile.txt", "data/right/sparse_parFile.txt");
	//searchTrail(threadNum, outputInfo, "data/left/sparse_msgFile.txt", "data/left/sparse_diffFile.txt", "data/left/sparse_parFile.txt");
	return 1;
}