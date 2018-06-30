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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputJuncReadCountFile outputSplicedAlignmentCount" << endl;
		exit(1);
	}
	string inputSJreadCountPath = argv[1];
	string outputSplicedAlignmentCountPath = argv[2];
	ifstream SJreadCount_ifs(inputSJreadCountPath.c_str());
	ofstream splicedAlignmentCount_ofs(outputSplicedAlignmentCountPath.c_str());
	int totalSplicedAlignmentNum = 0;
	while(!SJreadCount_ifs.eof())
	{
		string juncStr;
		getline(SJreadCount_ifs, juncStr);
		 if(SJreadCount_ifs.eof()||(juncStr == ""))
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
		string tmpSJsupportNumStr = juncFieldVec[4];
		int tmpSJsupportNumInt = atoi(tmpSJsupportNumStr.c_str());
		totalSplicedAlignmentNum += tmpSJsupportNumInt;
	}	
	splicedAlignmentCount_ofs << "inputSJreadCountPath: " << inputSJreadCountPath << endl << endl;
	splicedAlignmentCount_ofs << "totalSplicedAlignmentNum: " << totalSplicedAlignmentNum << endl;

	SJreadCount_ifs.close();
	splicedAlignmentCount_ofs.close();
	return 0;
}	