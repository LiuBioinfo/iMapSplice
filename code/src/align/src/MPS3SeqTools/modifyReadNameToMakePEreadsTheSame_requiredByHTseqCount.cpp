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
		cout << "Executable inputPEsam outputModifiedPEsam" << endl;
		exit(1);
	}
	string inputOriPeSAM = argv[1];
	string outputModifiedPeSAM = argv[2];
	ifstream oriPeSAM_ifs(inputOriPeSAM.c_str());
	ofstream modifiedPeSAM_ofs(outputModifiedPeSAM.c_str());
	while(!oriPeSAM_ifs.eof())
	{
		string tmpSAM_1;
		getline(oriPeSAM_ifs, tmpSAM_1);
		if(oriPeSAM_ifs.eof())
			break;
		if(tmpSAM_1.substr(0,1) == "@")
		{
			modifiedPeSAM_ofs << tmpSAM_1 << endl;
			continue;
		}
		string tmpSAM_2;
		getline(oriPeSAM_ifs, tmpSAM_2);

		int firstTabLoc_1 = tmpSAM_1.find("\t");
		int firstTabLoc_2 = tmpSAM_2.find("\t");
		string tmpOriReadName_1 = tmpSAM_1.substr(0, firstTabLoc_1);
		string tmpOriReadName_2 = tmpSAM_2.substr(0, firstTabLoc_2);
		string tmpOriSamOther_1 = tmpSAM_1.substr(firstTabLoc_1 + 1);
		string tmpOriSamOther_2 = tmpSAM_2.substr(firstTabLoc_2 + 1);
		int blankLoc_1 = tmpOriReadName_1.find(" ");
		int blankLoc_2 = tmpOriReadName_2.find(" ");
		string tmpModifiedReadName_1 = tmpSAM_1.substr(0, blankLoc_1);
		string tmpModifiedReadName_2 = tmpSAM_2.substr(0, blankLoc_2);
		modifiedPeSAM_ofs << tmpModifiedReadName_1 << "\t" << tmpOriSamOther_1 << endl;
		modifiedPeSAM_ofs << tmpModifiedReadName_2 << "\t" << tmpOriSamOther_2 << endl;
	}
	oriPeSAM_ifs.close();
	modifiedPeSAM_ofs.close();
	return 0;
}