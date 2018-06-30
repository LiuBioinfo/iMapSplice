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
		cout << "Executable specificAluFamilyName inputMergedAluFa outputSpecificAluFamilySeqFile" << endl;
		exit(1);
	}
	string specificAluFamilyName = argv[1];
	string inputMergedAluFa = argv[2];
	string outputSpecificAluFamilySeqFile = argv[3];
	ifstream mergedAlu_ifs(inputMergedAluFa.c_str());
	ofstream specificAlu_ofs(outputSpecificAluFamilySeqFile.c_str());
	while(!mergedAlu_ifs.eof())
	{
		string tmpAluName, tmpAluSeq;
		getline(mergedAlu_ifs, tmpAluName);
		if(tmpAluName == "")
			break;
		getline(mergedAlu_ifs, tmpAluSeq);
		int underlineLoc = tmpAluName.find("_");
		string tmpAluFamilyName = tmpAluName.substr(1, underlineLoc-1);
		cout << "tmpAluFamilyName: " << tmpAluFamilyName << endl;
		if(tmpAluFamilyName == specificAluFamilyName)
		{
			specificAlu_ofs << tmpAluName << endl;
			specificAlu_ofs << tmpAluSeq << endl;
		}
	}
	specificAlu_ofs.close();
	mergedAlu_ifs.close();
	return 0;
}