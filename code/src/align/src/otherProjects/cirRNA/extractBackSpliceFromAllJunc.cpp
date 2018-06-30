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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/alignInferJunctionHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputJuncFile outputFolder" << endl;
		exit(1);
	}
	string inputJuncFilePath = argv[1];
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());
	string outputBackSpliceStr = outputFolderStr + "backSpliceJunc.txt";
	ofstream backSJ_ofs(outputBackSpliceStr.c_str());
	string outputNormalSpliceStr = outputFolderStr + "normalSpliceJunc.txt";
	ofstream normalSJ_ofs(outputNormalSpliceStr.c_str());

	ifstream junc_ifs(inputJuncFilePath.c_str());
	while(!junc_ifs.eof())
	{
		string juncStr;
		getline(junc_ifs, juncStr);
		 if(junc_ifs.eof()||(juncStr == ""))
		 	break;
		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		juncFieldVec.push_back(juncStr.substr(startLoc));
		string tmpSJstartPosStr = juncFieldVec[1];
		string tmpSJendPosStr = juncFieldVec[2];
		int tmpSJstartPosInt = atoi(tmpSJstartPosStr.c_str());
		int tmpSJendPosInt = atoi(tmpSJendPosStr.c_str());		
		if(tmpSJstartPosInt > tmpSJendPosInt)
		{
			backSJ_ofs << juncStr << endl;
		}
		else
		{
			normalSJ_ofs << juncStr << endl;
		}
	}
	junc_ifs.close();
	backSJ_ofs.close();
	normalSJ_ofs.close();
	log_ofs.close();
	return 0;
}