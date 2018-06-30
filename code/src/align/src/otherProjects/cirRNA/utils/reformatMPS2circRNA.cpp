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
#include <hash_map>
#include <map>
#include <set>

#include "../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputMPS2circRNAresults outputReformatedCircRNAfile" << endl;
		exit(1);
	}	

	string inputFilePath = argv[1];
	string outputFilePath = argv[2];
	ifstream MPS2circRNA_ifs(inputFilePath.c_str());
	ofstream reformatedCircRNA_ofs(outputFilePath.c_str());
	while(!MPS2circRNA_ifs.eof())
	{
		string tmpStr;
		getline(MPS2circRNA_ifs, tmpStr);
		if((MPS2circRNA_ifs.eof())||(tmpStr == ""))
			break;		
		vector<string> backSpliceFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 5; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			string tmpFieldStr = tmpStr.substr(startLoc, tabLoc-startLoc);
			backSpliceFieldVec.push_back(tmpFieldStr);
			startLoc = tabLoc + 1;
		}
		string tmpChrInfoStr = backSpliceFieldVec[0];
		string tmpStartPosStr = backSpliceFieldVec[1];
		string tmpEndPosStr = backSpliceFieldVec[2];
		string tmpJuncIdStr = backSpliceFieldVec[3];
		string tmpSupNumStr = backSpliceFieldVec[4];
		int tlideLoc = tmpChrInfoStr.find("~");
		string tmpChrName = tmpChrInfoStr.substr(0, tlideLoc);
		//int tmpStartPos = atoi(tmpStartPosStr.c_str());
		//int tmpEndPos = atoi(tmpEndPosStr.c_str());
		//int tmpSupNum = atoi(tmpSupNumStr.c_str());
		reformatedCircRNA_ofs << tmpChrName << "\t" 
			<< tmpEndPosStr << "\t" << tmpStartPosStr << "\t" 
			<< tmpJuncIdStr << "\t" << tmpSupNumStr << endl;
	}
	MPS2circRNA_ifs.close();
	reformatedCircRNA_ofs.close();
	return 0;
}
