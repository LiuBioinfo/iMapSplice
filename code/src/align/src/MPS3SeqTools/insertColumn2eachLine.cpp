// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>

#include "../general/read_block_test.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputFile outputFile" << endl;
		exit(1);
	}
	string inputFileStr = argv[1];
	ifstream input_ifs(inputFileStr.c_str());
	string outputFileStr = argv[2];
	ofstream output_ofs(outputFileStr.c_str());
	int tmpNum = 0;
	while(1)
	{
		if(input_ifs.eof())
			break;
		string tmpLine;
		getline(input_ifs, tmpLine);
		if((input_ifs.eof())||(tmpLine == ""))
			break;
		//if(tmpLine.at(0) == '@')||()
		//	continue;
		tmpNum ++;
		vector<string> tmpLineFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 3; tmp++)
		{
			int tabLoc = tmpLine.find("\t", startLoc);
			string tmpLineField = tmpLine.substr(startLoc, tabLoc-startLoc);
			tmpLineFieldVec.push_back(tmpLineField);
			startLoc = tabLoc + 1;
		}
		string otherStr = tmpLine.substr(startLoc);
		tmpLineFieldVec.push_back(otherStr);
		output_ofs << tmpLineFieldVec[0] << "\t" << tmpLineFieldVec[1] << "\t"
			<< tmpLineFieldVec[2] << "\t" << "JUNC_" << int_to_str(tmpNum) << "\t" << otherStr << endl;	
	}

	return 0;
}