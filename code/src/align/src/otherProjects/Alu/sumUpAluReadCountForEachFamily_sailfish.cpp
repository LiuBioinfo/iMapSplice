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

void extractAluRegionInfoFromStr(string& tmpSailfishStr, 
	string& tmpAluFamilyName, double& tmpReadCount)
{
	int tabLoc_1 = tmpSailfishStr.find("\t");
	int tabLoc_2 = tmpSailfishStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpSailfishStr.find("\t", tabLoc_2 + 1);	
	int tabLoc_4 = tmpSailfishStr.find("\t", tabLoc_3 + 1);
	int tabLoc_5 = tmpSailfishStr.find("\t", tabLoc_4 + 1);	
	int tabLoc_6 = tmpSailfishStr.find("\t", tabLoc_5 + 1);	
	string tmpAluName = tmpSailfishStr.substr(0, tabLoc_1);
	int lineLoc = tmpAluName.find("_");
	tmpAluFamilyName = tmpAluName.substr(0, lineLoc);
	string tmpReadCountStr = tmpSailfishStr.substr(tabLoc_6+1);
	tmpReadCount = atof(tmpReadCountStr.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{	
		cout << "Executable inputSailfishResults outputSumUpResults" << endl;
		exit(1);
	}
	string inputSailfishResults = argv[1];
	string outputSumUpResults = argv[2];
	ifstream sailfish_ifs(inputSailfishResults.c_str());
	ofstream sumUp_ofs(outputSumUpResults.c_str());
	cout << "start to read sailfish results and do sum up " << endl;
	vector<string> aluNameVec;
	vector<double> aluReadCountVec;
	while(!sailfish_ifs.eof())
	{
		string tmpStr;
		getline(sailfish_ifs, tmpStr);
		if(tmpStr.substr(0,1) == "#")
			continue;
		if(tmpStr == "")
			break;
		cout << "tmpStr: " << tmpStr << endl;
		string tmpAluFamilyName;
		double tmpReadCount;
		extractAluRegionInfoFromStr(tmpStr, 
			tmpAluFamilyName, tmpReadCount);
		int matchedIndex = -1;
		for(int tmp = 0; tmp < aluNameVec.size(); tmp++)
		{
			string tmpNameInVec = aluNameVec[tmp];
			if(tmpNameInVec == tmpAluFamilyName)
			{
				matchedIndex = tmp;
				break;
			}
		}
		if(matchedIndex < 0)
		{
			aluNameVec.push_back(tmpAluFamilyName);
			aluReadCountVec.push_back(0);
		}
		else
		{
			double oriCount = aluReadCountVec[matchedIndex];
			double newCount = oriCount + tmpReadCount;
			aluReadCountVec[matchedIndex] = newCount;
		}
	}
	cout << "start to output sum up results " << endl;
	for(int tmp = 0; tmp < aluNameVec.size(); tmp++)
	{
		double tmpReadCount = aluReadCountVec[tmp];
		string tmpAluFamilyName = aluNameVec[tmp];
		sumUp_ofs << tmpAluFamilyName << "\t" << tmpReadCount << endl;
	}

	sailfish_ifs.close();
	sumUp_ofs.close();
	return 0;
}