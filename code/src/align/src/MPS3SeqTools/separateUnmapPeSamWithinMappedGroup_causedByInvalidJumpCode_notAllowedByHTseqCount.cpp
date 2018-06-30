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

bool mappedOrNot(int tmpFlag)
{
	if(tmpFlag & 0x4)
		return false;
	else
		return true;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputPEsam outputSeparatedFilePrefix" << endl;
		exit(1);
	}
	string inputOriPeSAM = argv[1];
	string outputSeparatedFilePrefix = argv[2];
	string outputSeparatedFile_valid = outputSeparatedFilePrefix + "_valid.sam";
	string outputSeparatedFile_invalid = outputSeparatedFilePrefix + "_invalid.sam";
	ifstream oriPeSAM_ifs(inputOriPeSAM.c_str());
	ofstream validPeSAM_ofs(outputSeparatedFile_valid.c_str());
	ofstream invalidPeSAM_ofs(outputSeparatedFile_invalid.c_str());
	while(!oriPeSAM_ifs.eof())
	{
		string tmpSAM_1;
		getline(oriPeSAM_ifs, tmpSAM_1);
		if(oriPeSAM_ifs.eof())
			break;
		if(tmpSAM_1.substr(0,1) == "@")
		{
			validPeSAM_ofs << tmpSAM_1 << endl;
			continue;
		}
		string tmpSAM_2;
		getline(oriPeSAM_ifs, tmpSAM_2);

		int firstTabLoc_1 = tmpSAM_1.find("\t");
		int firstTabLoc_2 = tmpSAM_2.find("\t");
		// string tmpOriReadName_1 = tmpSAM_1.substr(0, firstTabLoc_1);
		// string tmpOriReadName_2 = tmpSAM_2.substr(0, firstTabLoc_2);
		// string tmpOriSamOther_1 = tmpSAM_1.substr(firstTabLoc_1 + 1);
		// string tmpOriSamOther_2 = tmpSAM_2.substr(firstTabLoc_2 + 1);
		// int blankLoc_1 = tmpOriReadName_1.find(" ");
		// int blankLoc_2 = tmpOriReadName_2.find(" ");
		// string tmpModifiedReadName_1 = tmpSAM_1.substr(0, blankLoc_1);
		// string tmpModifiedReadName_2 = tmpSAM_2.substr(0, blankLoc_2);
		// modifiedPeSAM_ofs << tmpModifiedReadName_1 << "\t" << tmpOriSamOther_1 << endl;
		// modifiedPeSAM_ofs << tmpModifiedReadName_2 << "\t" << tmpOriSamOther_2 << endl;
		int secondTabLoc_1 = tmpSAM_1.find("\t", firstTabLoc_1 + 1);
		int secondTabLoc_2 = tmpSAM_2.find("\t", firstTabLoc_2 + 1);
		string flagStr_1 = tmpSAM_1.substr(firstTabLoc_1 + 1, secondTabLoc_1 - 1 - firstTabLoc_1 - 1 + 1);
		string flagStr_2 = tmpSAM_2.substr(firstTabLoc_2 + 1, secondTabLoc_2 - 1 - firstTabLoc_2 - 1 + 1);
		int flag_1 = atoi(flagStr_1.c_str());
		int flag_2 = atoi(flagStr_2.c_str());
		bool mappedOrNotBool_1 = mappedOrNot(flag_1);
		bool mappedOrNotBool_2 = mappedOrNot(flag_2);
		if((mappedOrNotBool_1 && (!mappedOrNotBool_2))
			||((!mappedOrNotBool_1) && mappedOrNotBool_2))
		{
			invalidPeSAM_ofs << tmpSAM_1 << endl;
			invalidPeSAM_ofs << tmpSAM_2 << endl;
		}
		else
		{
			validPeSAM_ofs << tmpSAM_1 << endl;
			validPeSAM_ofs << tmpSAM_2 << endl;
		}

	}
	oriPeSAM_ifs.close();
	validPeSAM_ofs.close();
	invalidPeSAM_ofs.close();
	return 0;
}