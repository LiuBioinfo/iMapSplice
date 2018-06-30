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

#include "../../../general/read_block_test.h"
//#include "../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputValidatedBackSpliceJuncList inputSam2alignInferJuncHashFile outputFolderPath" << endl;
		exit(1);
	}
	string inputValidatedBackSpliceJuncList = argv[1];
	string inputSam2alignInferJuncHashFile = argv[2];
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());
	string outputValidatedSpliceDetectionResultsFilePath = outputFolderStr + "detectionResults.txt";
	ofstream detectionResults_ofs(outputValidatedSpliceDetectionResultsFilePath.c_str());
	string outputFoundValidatedBackSpliceFilePath = outputFolderStr + "foundBackSplice.txt";
	ofstream foundBackSplice_ofs(outputFoundValidatedBackSpliceFilePath.c_str());
	string outputUnfoundValidatedBackSpliceFilePath = outputFolderStr + "unfoundBackSplice.txt";
	ofstream unfoundBackSplice_ofs(outputUnfoundValidatedBackSpliceFilePath.c_str());
	string outputOtherSpliceFilePath = outputFolderStr + "otherSplice.txt";
	ofstream otherSplice_ofs(outputOtherSpliceFilePath.c_str());

	int validatedBackSpliceNum = 0;
	vector<int> supportNumInDetectedSpliceResultsVec;

	vector<string> juncNameVec_validatedBackSplice;
	vector<string> chrNameStrVec_validatedBackSplice;
	vector<int> startPosVec_validatedBackSplice;
	vector<int> endPosVec_validatedBackSplice;

	ifstream validatedBackSpliceList_ifs(inputValidatedBackSpliceJuncList.c_str());
	while(!validatedBackSpliceList_ifs.eof())
	{
		string backSpliceStr;
		getline(validatedBackSpliceList_ifs, backSpliceStr);
		if(backSpliceStr == "")
			break;
		vector<string> backSpliceFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = backSpliceStr.find("\t", startLoc);
			string tmpBackSpliceField = backSpliceStr.substr(startLoc, tabLoc-startLoc);
			backSpliceFieldVec.push_back(tmpBackSpliceField);
			startLoc = tabLoc + 1;
		}
		string tmpBackSpliceNameStr = backSpliceFieldVec[0];
		string tmpBackSpliceChrNameStr = backSpliceFieldVec[1];
		string tmpBackSpliceEndPosStr = backSpliceFieldVec[2];
		string tmpBackSpliceStartPosStr = backSpliceFieldVec[3];
		int tmpBackSpliceEndPosInt = atoi(tmpBackSpliceEndPosStr.c_str());
		int tmpBackSpliceStartPosInt = atoi(tmpBackSpliceStartPosStr.c_str());
		supportNumInDetectedSpliceResultsVec.push_back(0);
		juncNameVec_validatedBackSplice.push_back(tmpBackSpliceNameStr);
		chrNameStrVec_validatedBackSplice.push_back(tmpBackSpliceChrNameStr);
		startPosVec_validatedBackSplice.push_back(tmpBackSpliceStartPosInt);
		endPosVec_validatedBackSplice.push_back(tmpBackSpliceEndPosInt);
		validatedBackSpliceNum ++;
	}
	validatedBackSpliceList_ifs.close();

	cout << "validatedBackSpliceNum: " << validatedBackSpliceNum << endl;
	log_ofs << "validatedBackSpliceNum: " << validatedBackSpliceNum << endl;
	
	ifstream sam2alignInferJuncHash_ifs(inputSam2alignInferJuncHashFile.c_str());
	while(!sam2alignInferJuncHash_ifs.eof())
	{
		string tmpSpliceStr;
		getline(sam2alignInferJuncHash_ifs, tmpSpliceStr);
		if(tmpSpliceStr == "")
			break;
		vector<string> tmpSpliceFieldVec;
		int startLoc2 = 0;
		for(int tmp2 = 0; tmp2 < 4; tmp2++)	
		{
			int tabLoc2 = tmpSpliceStr.find("\t", startLoc2);
			string tmpSpliceField = tmpSpliceStr.substr(startLoc2, tabLoc2 - startLoc2);
			tmpSpliceFieldVec.push_back(tmpSpliceField);
			startLoc2 = tabLoc2 + 1;
		}
		tmpSpliceFieldVec.push_back(tmpSpliceStr.substr(startLoc2));
		string tmpSpliceChrNameStr = tmpSpliceFieldVec[0];
		string tmpSpliceStartPosStr = tmpSpliceFieldVec[1];
		string tmpSpliceEndPosStr = tmpSpliceFieldVec[2];
		string tmpSpliceSupportNumStr = tmpSpliceFieldVec[4];
		int tmpSpliceStartPosInt = atoi(tmpSpliceStartPosStr.c_str());
		int tmpSpliceEndPosInt = atoi(tmpSpliceEndPosStr.c_str());
		int tmpSpliceSupportNumInt = atoi(tmpSpliceSupportNumStr.c_str());
		bool tmpFoundInValidatedSpliceListBool = false;
		int tmpFoundValidatedBackSpliceIndex = -1;
		for(int tmpSplice = 0; tmpSplice < validatedBackSpliceNum; tmpSplice++)
		{
			if((tmpSpliceChrNameStr == chrNameStrVec_validatedBackSplice[tmpSplice])
				&& (tmpSpliceStartPosInt == startPosVec_validatedBackSplice[tmpSplice])
				&& (tmpSpliceEndPosInt == endPosVec_validatedBackSplice[tmpSplice]))
			{
				tmpFoundValidatedBackSpliceIndex = tmpSplice;
				tmpFoundInValidatedSpliceListBool = true;
			}
		}
		if(tmpFoundInValidatedSpliceListBool)
		{
			foundBackSplice_ofs << tmpSpliceStr << endl;
			supportNumInDetectedSpliceResultsVec[tmpFoundValidatedBackSpliceIndex] 
				= tmpSpliceSupportNumInt;
		}
		else
		{
			if(tmpSpliceStartPosInt > tmpSpliceEndPosInt)
				unfoundBackSplice_ofs << tmpSpliceStr << endl;
			else
				otherSplice_ofs << tmpSpliceStr << endl;
		}
	}
	sam2alignInferJuncHash_ifs.close();

	for(int tmp = 0; tmp < validatedBackSpliceNum; tmp++)
	{
		detectionResults_ofs << juncNameVec_validatedBackSplice[tmp] << "\t"
			<< chrNameStrVec_validatedBackSplice[tmp] << "\t"
			<< endPosVec_validatedBackSplice[tmp] << "\t"
			<< startPosVec_validatedBackSplice[tmp] << "\t"
			<< supportNumInDetectedSpliceResultsVec[tmp] << endl; 
	}

	log_ofs.close();
	detectionResults_ofs.close();
	foundBackSplice_ofs.close();
	unfoundBackSplice_ofs.close();
	otherSplice_ofs.close();
	return 0;
}