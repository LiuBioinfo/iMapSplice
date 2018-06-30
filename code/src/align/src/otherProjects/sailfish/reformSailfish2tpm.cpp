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
		cout << "Executable inputSailfishResults outputFile" << endl;
		exit(1);
	}
	string inputSailfishResults = argv[1];
	string outputFile = argv[2];

	ifstream sf_ifs(inputSailfishResults.c_str());
	ofstream tpm_ofs(outputFile.c_str());
	string str_1;
	getline(sf_ifs, str_1); // 1st line
	while(!sf_ifs.eof())
	{
		string tmpSfStr;
		getline(sf_ifs, tmpSfStr);
		if(tmpSfStr == "")
			break;
		int tabLoc_1 = tmpSfStr.find("\t");
		int tabLoc_2 = tmpSfStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpSfStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpSfStr.find("\t", tabLoc_3 + 1);
		string tmpTranscriptFullName = tmpSfStr.substr(0, tabLoc_1);
		int straightLineLoc = tmpTranscriptFullName.find("|");
		string tmpTranscriptName = tmpTranscriptFullName.substr(0, straightLineLoc);
		string tmpTpmStr = tmpSfStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		tpm_ofs << tmpTranscriptName << "\t" << tmpTpmStr << endl;
	}
	sf_ifs.close();
	tpm_ofs.close();
	return 0;
}