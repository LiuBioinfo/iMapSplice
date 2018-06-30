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
	if(argc != 4)
	{
		cout << "Executable inputJuncListFile lowSupNumThreshold outputJuncFilePrefix" << endl;
		exit(1);
	}

	string inputJuncListFile = argv[1];
	string lowSupNumThresholdStr = argv[2];
	int lowSupNumThreshold = atoi(lowSupNumThresholdStr.c_str());
	string outputFilePrefix = argv[3];
	string output_lowSupJunc = outputFilePrefix + "lowSup.junc";
	string output_highSupJunc = outputFilePrefix + "highSup.junc";
	ofstream lowSup_ofs(output_lowSupJunc.c_str());
	ofstream highSup_ofs(output_highSupJunc.c_str());
	ifstream junc_ifs(inputJuncListFile.c_str());
	while(!junc_ifs.eof())
	{
		string tmpJuncStr;
		getline(junc_ifs, tmpJuncStr);
		if(tmpJuncStr == "")
			break;
		int tabLoc_1 = tmpJuncStr.find("\t");
		int tabLoc_2 = tmpJuncStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpJuncStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpJuncStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpJuncStr.find("\t", tabLoc_4 + 1);
		string tmpJuncSupNumStr = tmpJuncStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		int tmpJuncSupNum = atoi(tmpJuncSupNumStr.c_str());
		if(tmpJuncSupNum <= lowSupNumThreshold)
			lowSup_ofs << tmpJuncStr << endl;
		else
			highSup_ofs << tmpJuncStr << endl;
	}
	junc_ifs.close();
	lowSup_ofs.close();
	highSup_ofs.close();
	return 0;
}