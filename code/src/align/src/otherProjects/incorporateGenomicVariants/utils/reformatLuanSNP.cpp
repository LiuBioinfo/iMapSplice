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
#include <sstream>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputLuanSNPfile outputReformattedSNPfile" << endl;
		exit(1);
	}
	string inputLuanSNPfile = argv[1];
	string outputReformattedSNPfile = argv[2];
	ifstream LuanSNP_ifs(inputLuanSNPfile.c_str());
	ofstream reformatSNP_ofs(outputReformattedSNPfile.c_str());
	while(!LuanSNP_ifs.eof())
	{
		string tmpStr;
		getline(LuanSNP_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		string tmpChr = tmpStr.substr(0, tabLoc_1);
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		string tmpPos = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpBaseChange;
		if(tabLoc_3 != string::npos)
			tmpBaseChange = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		else
			tmpBaseChange = tmpStr.substr(tabLoc_2 + 1);
		string tmpRefBase = tmpBaseChange.substr(0,1);
		string tmpAlterBase = tmpBaseChange.substr(3,1);
		reformatSNP_ofs << tmpChr << "\t" << tmpPos << "\t" << tmpRefBase << "\t" << tmpAlterBase << endl;
	}
	LuanSNP_ifs.close();
	reformatSNP_ofs.close();
	return 0;
}