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
		cout << "Executable inputIndividualVCF outputFolderPath" << endl;
		exit(1);
	}
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());

	string inputIndividualVCF = argv[1];
	ifstream individualVCF_ifs();
	while(!individualVCF_ifs.eof())
	{
		string tmpVCFstr;
		getline(individualVCF_ifs, tmpVCFstr);
		if(tmpVCFstr == "")
			break;
		if(tmpVCFstr.substr(0,1) == "#")
			continue;
		vector<string> tmpVCFfieldVec;

		string tmpVar_chrName = "chr" + tmpVCFfieldVec[0];
		string tmpVar_chrPosStr = tmpVCFfieldVec[1];
		int tmpVar_chrPos = atoi(tmpVar_chrPosStr.c_str());
		string tmpVar_ID = tmpVCFfieldVec[2];
		string tmpVar_refBase = tmpVCFfieldVec[3];
		string tmpVar_altBase = tmpVCFfieldVec[4];
		string tmpVar_qualStr = tmpVCFfieldVec[5];
		double tmpVar_qual = atof(tmpVar_qualStr.c_str());
		string tmpVar_filter = tmpVCFfieldVec[6];
		string tmpVar_format = tmpVCFfieldVec[8];
		string tmpVar_baseLabelStr = tmpVCFfieldVec[9];
	}
	individualVCF_ifs.close();
	log_ofs.close();
	return 0;
}