// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

#include "../../../general/read_block_test.h"

using namespace std;

void extractRealFusionDetectedFusionInfoFromFusionBreakPointStr(
	string& tmpFusionName, 
	string& realChrName_1, string& realChrName_2, 
	int& realBreakPoint_1, int& realBreakPoint_2, 
	string& detectedChrName_1, string& detectedChrName_2,
	int& detectedBreakPoint_1, int& detectedBreakPoint_2, string& tmpFusionStr)
{
	int midTabLoc = tmpFusionStr.find("\t", 0);
	string tmpFusionStr_1stPart = tmpFusionStr.substr(0,midTabLoc);

	int secondTabLoc = (tmpFusionStr.substr(midTabLoc+1)).find("\t", 0);
	string tmpFusionStr_2ndPart;
	if(secondTabLoc == string::npos)
		tmpFusionStr_2ndPart = tmpFusionStr.substr(midTabLoc+1);
	else
		tmpFusionStr_2ndPart = tmpFusionStr.substr(midTabLoc + 1, secondTabLoc);
	vector<string> tmpRealFusionFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 6; tmp++)
	{
		int tabLoc = tmpFusionStr_1stPart.find("_", startLoc);
		string tmpFusionField = tmpFusionStr_1stPart.substr(startLoc, tabLoc - startLoc);
		tmpRealFusionFieldVec.push_back(tmpFusionField);
		startLoc = tabLoc + 1;
	}
	realChrName_1 = tmpRealFusionFieldVec[1];
	realBreakPoint_1 = atoi(tmpRealFusionFieldVec[2].c_str());
	realChrName_2 = tmpRealFusionFieldVec[4];
	realBreakPoint_2 = atoi(tmpRealFusionFieldVec[5].c_str());
	tmpFusionName = tmpRealFusionFieldVec[0] + "_" + tmpRealFusionFieldVec[1]
		+ "_" + tmpRealFusionFieldVec[2] + "_" + tmpRealFusionFieldVec[3]
		+ "_" + tmpRealFusionFieldVec[4] + "_" + tmpRealFusionFieldVec[5];

	int midCharCloc = tmpFusionStr_2ndPart.find("c",2);
	string tmpDetectedFusionInfo_1 = tmpFusionStr_2ndPart.substr(0, midCharCloc);
	string tmpDetectedFusionInfo_2 = tmpFusionStr_2ndPart.substr(midCharCloc);
	int commaLoc_1 = tmpDetectedFusionInfo_1.find(":");
	int commaLoc_2 = tmpDetectedFusionInfo_2.find(":");
	detectedChrName_1 = tmpDetectedFusionInfo_1.substr(0,commaLoc_1);
	detectedChrName_2 = tmpDetectedFusionInfo_2.substr(0,commaLoc_2);
	detectedBreakPoint_1 = atoi((tmpDetectedFusionInfo_1.substr(commaLoc_1 + 1)).c_str());
	detectedBreakPoint_2 = atoi((tmpDetectedFusionInfo_2.substr(commaLoc_2 + 1)).c_str());
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputFusionNameFromFusionCancer inputFusionBreakPointResults outputFolder" << endl;
		exit(1);
	}
	string inputFusionNameFromFusionCancer = argv[1];
	string inputFusionBreakPointResults = argv[2];
	ifstream fusion_ifs(inputFusionNameFromFusionCancer.c_str());
	ifstream fusionBreakPointResults_ifs(inputFusionBreakPointResults.c_str());

	string outputFolder = argv[3];
	outputFolder += "/";
	string mkdir= "mkdir -p " + outputFolder;
	system(mkdir.c_str());	
	string output_detectedFusion = outputFolder + "/detectedFusion.txt";
	string output_missedFusion = outputFolder + "/missedFusion.txt";
	string output_correctFusionBreakPoint = outputFolder + "/correctFusionBreakPoint.txt";		
	string output_incorrectFusionBreakPoint = outputFolder + "/incorrectFusionBreakPoint.txt";
	
	ofstream detectedFusion_ofs(output_detectedFusion.c_str());
	ofstream missedFusion_ofs(output_missedFusion.c_str());
	ofstream correctFusionBreakPoint_ofs(output_correctFusionBreakPoint.c_str());
	ofstream incorrectFusionBreakPoint_ofs(output_incorrectFusionBreakPoint.c_str());

	set<string> fusionNameSet;
	while(!fusion_ifs.eof())
	{
		string tmpFusionStr;
		getline(fusion_ifs, tmpFusionStr);
		if((fusion_ifs.eof())||(tmpFusionStr == ""))
			break;
		string tmpFusionName = tmpFusionStr.substr(1);
		fusionNameSet.insert(tmpFusionName);
	}

	set<string> detectedFusionNameSet;
	while(!fusionBreakPointResults_ifs.eof())
	{
		string tmpFusionStr;
		getline(fusionBreakPointResults_ifs, tmpFusionStr);
		if((fusionBreakPointResults_ifs.eof())||(tmpFusionStr == ""))
			break;		
		string tmpFusionName;
		string realChrName_1;
		string realChrName_2;
		int realBreakPoint_1;
		int realBreakPoint_2;
		string detectedChrName_1;
		string detectedChrName_2;
		int detectedBreakPoint_1;
		int detectedBreakPoint_2;
		extractRealFusionDetectedFusionInfoFromFusionBreakPointStr(
			tmpFusionName, realChrName_1, realChrName_2,
			realBreakPoint_1, realBreakPoint_2,
			detectedChrName_1, detectedChrName_2,
			detectedBreakPoint_1, detectedBreakPoint_2,
			tmpFusionStr);
		// cout << "realChrName_1: " << realChrName_1 << endl;
		// cout << "realChrName_2: " << realChrName_2 << endl;
		// cout << "detectedChrName_1: " << detectedChrName_1 << endl;
		// cout << "detectedChrName_2: " << detectedChrName_2 << endl;
		// cout << "realBreakPoint_1: " << realBreakPoint_1 << endl;
		// cout << "realBreakPoint_2: " << realBreakPoint_2 << endl;
		// cout << "detectedBreakPoint_1: " << detectedBreakPoint_1 << endl;
		// cout << "detectedBreakPoint_2: " << detectedBreakPoint_2 << endl;
		// if(((realChrName_1 == detectedChrName_1)&&((realBreakPoint_1 - 300000 <= detectedBreakPoint_1)&&(detectedBreakPoint_1 <= realBreakPoint_1 + 300000)))
		// 	||((realChrName_2 == detectedChrName_2)&&((realBreakPoint_2 - 300000 <= detectedBreakPoint_2)&&(detectedBreakPoint_2 <= realBreakPoint_2 + 300000)))  
		// 	||((realChrName_1 == detectedChrName_2)&&((realBreakPoint_1 - 300000 <= detectedBreakPoint_2)&&(detectedBreakPoint_2 <= realBreakPoint_1 + 300000)))
		// 	||((realChrName_2 == detectedChrName_1)&&((realBreakPoint_2 - 300000 <= detectedBreakPoint_1)&&(detectedBreakPoint_1 <= realBreakPoint_2 + 300000)))
		// 		)

		if((((realChrName_1 == detectedChrName_1)&&((realBreakPoint_1 - 300000 <= detectedBreakPoint_1)&&(detectedBreakPoint_1 <= realBreakPoint_1 + 300000)))
			&&((realChrName_2 == detectedChrName_2)&&((realBreakPoint_2 - 300000 <= detectedBreakPoint_2)&&(detectedBreakPoint_2 <= realBreakPoint_2 + 300000)))  
			)
			||
			(((realChrName_1 == detectedChrName_2)&&((realBreakPoint_1 - 300000 <= detectedBreakPoint_2)&&(detectedBreakPoint_2 <= realBreakPoint_1 + 300000)))
			&&((realChrName_2 == detectedChrName_1)&&((realBreakPoint_2 - 300000 <= detectedBreakPoint_1)&&(detectedBreakPoint_1 <= realBreakPoint_2 + 300000)))
				))
		{
			detectedFusionNameSet.insert(tmpFusionName);
			correctFusionBreakPoint_ofs << tmpFusionStr << endl;
		}
		else
		{
			incorrectFusionBreakPoint_ofs << tmpFusionStr << endl;
		}
	}

	for(set<string>::iterator strSetIter = detectedFusionNameSet.begin();
		strSetIter != detectedFusionNameSet.end(); strSetIter++)
	{
		string tmpFusionName = (*strSetIter);
		detectedFusion_ofs << tmpFusionName << endl;
		//detectedFusionNameSet.insert(tmpFusionName);
	}

	for(set<string>::iterator strSetIter = fusionNameSet.begin();
		strSetIter != fusionNameSet.end(); strSetIter++)
	{
		string tmpFusionName = (*strSetIter);
		if(detectedFusionNameSet.find(tmpFusionName) == detectedFusionNameSet.end())
			missedFusion_ofs << tmpFusionName << endl;
	}

	return 0;
}