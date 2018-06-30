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
#include "../../general/read_block_test.h"
#include "../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputFeatureCountsResults outputFile" << endl;
		exit(1);
	}
	string inputFeatureCountsResults = argv[1];
	string outputFile = argv[2];

	ifstream fc_ifs(inputFeatureCountsResults.c_str());
	ofstream geneCount_ofs(outputFile.c_str());
	string str_1, str_2;
	getline(fc_ifs, str_1); // 1st line
	getline(fc_ifs, str_2); // 2nd line
	while(!fc_ifs.eof())
	{
		string tmpFcStr;
		getline(fc_ifs, tmpFcStr);
		if(tmpFcStr == "")
			break;
		int firstTabLoc = tmpFcStr.find("\t");
		int lastTabLoc = tmpFcStr.rfind("\t");
		//cout << "lastTabLoc: " << lastTabLoc << endl;
		string tmpGeneName = tmpFcStr.substr(0, firstTabLoc);
		string tmpCount = tmpFcStr.substr(lastTabLoc + 1);
		//cout << "tmpCount: " << tmpCount << endl;
		geneCount_ofs << tmpGeneName << "\t" << tmpCount << endl;
	}
	fc_ifs.close();
	geneCount_ofs.close();
	return 0;
}