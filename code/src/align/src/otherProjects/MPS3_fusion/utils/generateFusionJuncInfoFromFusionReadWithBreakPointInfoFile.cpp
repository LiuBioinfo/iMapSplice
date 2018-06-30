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
#include "../../../general/index_info.h"

using namespace std;

typedef map<int, int> FusionJuncEndPos2supNumMap;
typedef map<int, FusionJuncEndPos2supNumMap > FusionJuncPosPairSupNumMap;

// bool searchForCanonicalFusionBreakPointWithOffset(
// 	bool leftGenePoint_NorOrRcm_bool, bool rightGene_NorOrRcm_bool,// here, left & right indicates the location in read seq.
// 	int leftChrNameInt, int rightChrNameInt,
// 	int leftBreakPointPos, int rightBreakPointPos,
// 	int& adjustedCanonicalLeftFusionBreakPoint,
// 	int& adjustedCanonicalRightFusionBreakPoint, Index_Info* indexInfo, int offset_max)
// {
// 	string rawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
// 		leftChrNameInt, rightChrNameInt,
// 		leftBreakPointPos, rightBreakPointPos,
// 		leftGenePoint_NorOrRcm_bool, rightGene_NorOrRcm_bool);
// 	vector<int> candiOffsetVec;
// 	for(int tmp = 0; tmp < offset_max; tmp++)
// 	{	
// 		candiOffsetVec.push_back(tmp+1);
// 		candiOffsetVec.push_back(-1-tmp);
// 	}
// 	if(leftGenePoint_NorOrRcm_bool && rightGene_NorOrRcm_bool)
// 	{
// 		if((rawFusionJuncFlankString == "GTAG")||(rawFusionJuncFlankString == "CTAC"))
// 			return false;
// 		for(int tmpOffsetIndex = 0; tmpOffsetIndex < offset_max*2; tmpOffsetIndex ++)
// 		{
// 			int tmpOffset = candiOffsetVec[tmpOffsetIndex];
// 			string tmpRawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
// 				leftChrNameInt, rightChrNameInt,
// 				leftBreakPointPos + tmpOffset, rightBreakPointPos + tmpOffset,
// 				true, true);
// 			if((tmpRawFusionJuncFlankString == "GTAG")||(tmpRawFusionJuncFlankString == "CTAC"))
// 			{
// 				adjustedCanonicalLeftFusionBreakPoint = leftBreakPointPos + tmpOffset;
// 				adjustedCanonicalRightFusionBreakPoint = rightBreakPointPos + tmpOffset;
// 				return true;
// 			}
// 		}	
// 	}
// 	else if((!leftGenePoint_NorOrRcm_bool) && (!rightGene_NorOrRcm_bool))
// 	{
// 		if((rawFusionJuncFlankString == "AGGT")||(rawFusionJuncFlankString == "ACCT"))
// 			return false;
// 		for(int tmpOffsetIndex = 0; tmpOffsetIndex < offset_max*2; tmpOffsetIndex ++)
// 		{
// 			int tmpOffset = candiOffsetVec[tmpOffsetIndex];
// 			string tmpRawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
// 				leftChrNameInt, rightChrNameInt,
// 				leftBreakPointPos + tmpOffset, rightBreakPointPos + tmpOffset,
// 				false, false);
// 			if((tmpRawFusionJuncFlankString == "AGGT")||(tmpRawFusionJuncFlankString == "ACCT"))
// 			{
// 				adjustedCanonicalLeftFusionBreakPoint = leftBreakPointPos + tmpOffset;
// 				adjustedCanonicalRightFusionBreakPoint = rightBreakPointPos + tmpOffset;
// 				return true;
// 			}
// 		}
// 	}
// 	else if(leftGenePoint_NorOrRcm_bool && (!rightGene_NorOrRcm_bool)) // For Rev
// 	{
// 		if((rawFusionJuncFlankString == "GTCT")||(rawFusionJuncFlankString == "CTGT"))
// 			return false;
// 		for(int tmpOffsetIndex = 0; tmpOffsetIndex < offset_max*2; tmpOffsetIndex ++)
// 		{
// 			int tmpOffset = candiOffsetVec[tmpOffsetIndex];
// 			string tmpRawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
// 				leftChrNameInt, rightChrNameInt,
// 				leftBreakPointPos + tmpOffset, rightBreakPointPos + tmpOffset,
// 				true, false);
// 			if((tmpRawFusionJuncFlankString == "GTCT")||(tmpRawFusionJuncFlankString == "CTGT"))
// 			{
// 				adjustedCanonicalLeftFusionBreakPoint = leftBreakPointPos + tmpOffset;
// 				adjustedCanonicalRightFusionBreakPoint = rightBreakPointPos - tmpOffset;
// 				return true;
// 			}
// 		}
// 	}
// 	else // REV FOR
// 	{
// 		if((rawFusionJuncFlankString == "ACAG")||(rawFusionJuncFlankString == "AGAC"))
// 			return false;
// 		for(int tmpOffsetIndex = 0; tmpOffsetIndex < offset_max*2; tmpOffsetIndex ++)
// 		{
// 			int tmpOffset = candiOffsetVec[tmpOffsetIndex];
// 			string tmpRawFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
// 				leftChrNameInt, rightChrNameInt,
// 				leftBreakPointPos + tmpOffset, rightBreakPointPos + tmpOffset,
// 				true, false);
// 			if((tmpRawFusionJuncFlankString == "ACAG")||(tmpRawFusionJuncFlankString == "AGAC"))
// 			{
// 				adjustedCanonicalLeftFusionBreakPoint = leftBreakPointPos + tmpOffset;
// 				adjustedCanonicalRightFusionBreakPoint = rightBreakPointPos - tmpOffset;
// 				return true;
// 			}
// 		}
// 	}
// }

void extractFusionInfoFromFusionBreakPointStr(
	int& detectedChrNameInt_1, int& detectedChrNameInt_2,
	int& detectedBreakPoint_1, int& detectedBreakPoint_2, 
	Index_Info* indexInfo, string& tmpFusionStr)
{
	int midTabLoc = tmpFusionStr.find("\t", 0);
	string tmpFusionStr_1stPart = tmpFusionStr.substr(0,midTabLoc);

	int secondTabLoc = (tmpFusionStr.substr(midTabLoc+1)).find("\t", 0);
	string tmpFusionStr_2ndPart;
	if(secondTabLoc == string::npos)
		tmpFusionStr_2ndPart = tmpFusionStr.substr(midTabLoc+1);
	else
		tmpFusionStr_2ndPart = tmpFusionStr.substr(midTabLoc + 1, secondTabLoc);
	//cout << "tmpFusionStr_2ndPart: " << tmpFusionStr_2ndPart << endl;
	int midCharCloc = tmpFusionStr_2ndPart.find("c",2);
	string tmpDetectedFusionInfo_1 = tmpFusionStr_2ndPart.substr(0, midCharCloc-1);
	string tmpDetectedFusionInfo_2 = tmpFusionStr_2ndPart.substr(midCharCloc);
	//cout << "tmpDetectedFusionInfo_1: " << tmpDetectedFusionInfo_1 << endl;
	//cout << "tmpDetectedFusionInfo_2: " << tmpDetectedFusionInfo_2 << endl;
	int commaLoc_1 = tmpDetectedFusionInfo_1.find(":");
	int commaLoc_2 = tmpDetectedFusionInfo_2.find(":");
	string detectedChrName_1 = tmpDetectedFusionInfo_1.substr(0,commaLoc_1);
	string detectedChrName_2 = tmpDetectedFusionInfo_2.substr(0,commaLoc_2);
	//cout << "detectedChrName_1: " << detectedChrName_1 << endl;
	//cout << "detectedChrName_2: " << detectedChrName_2 << endl;
	detectedChrNameInt_1 = indexInfo->convertStringToInt(detectedChrName_1);
	detectedChrNameInt_2 = indexInfo->convertStringToInt(detectedChrName_2);
	//cout << "detectedChrNameInt_1: " << detectedChrNameInt_1 << endl;
	detectedBreakPoint_1 = atoi((tmpDetectedFusionInfo_1.substr(commaLoc_1 + 1)).c_str());
	detectedBreakPoint_2 = atoi((tmpDetectedFusionInfo_2.substr(commaLoc_2 + 1)).c_str());	
}
		
void insertUpdataFusionJuncMapVecVecWithDetectedFusionJunc(
	vector < vector< FusionJuncPosPairSupNumMap > >& fusionJuncMapVecVec, 
	int detectedChrNameInt_1, int detectedChrNameInt_2,
	int detectedBreakPoint_1, int detectedBreakPoint_2)
{
	FusionJuncPosPairSupNumMap::iterator tmp1stMapIter 
		= ((fusionJuncMapVecVec[detectedChrNameInt_1])[detectedChrNameInt_2]).find(detectedBreakPoint_1);
	if(tmp1stMapIter == ((fusionJuncMapVecVec[detectedChrNameInt_1])[detectedChrNameInt_2]).end())
	{
		FusionJuncEndPos2supNumMap tmp2ndMap;
		tmp2ndMap.insert(pair<int,int>(detectedBreakPoint_2, 1));
		((fusionJuncMapVecVec[detectedChrNameInt_1])[detectedChrNameInt_2]).insert(
			pair<int,FusionJuncEndPos2supNumMap>(detectedBreakPoint_1, tmp2ndMap));
	}
	else
	{
		FusionJuncEndPos2supNumMap::iterator tmp2ndMapIter
			= (tmp1stMapIter->second).find(detectedBreakPoint_2);
		if(tmp2ndMapIter == (tmp1stMapIter->second).end())
			(tmp1stMapIter->second).insert(pair<int,int>(detectedBreakPoint_2, 1));
		else
			(tmp2ndMapIter->second)++;
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexInfoPath inputFuionReadWithBreakPointInfoFilePath outputFusionJuncPath " << endl;
		exit(1);
	}

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "end of initiating indexInfo" << endl;
	string inputFuionReadWithBreakPointInfoFilePath = argv[2];
	ifstream fusionBreakPointResults_ifs(inputFuionReadWithBreakPointInfoFilePath.c_str());

	string outputPath = argv[3];
	ofstream fusionJunc_ofs(outputPath.c_str());
	cout << "chromNum: " << chromNum << endl;
	vector < vector< FusionJuncPosPairSupNumMap > > fusionJuncMapVecVec;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++ )
	{
		//cout << "tmpChr: " << tmpChr << endl; 
		vector<FusionJuncPosPairSupNumMap> tmpFusionJuncMapVec;
		for(int tmpChr2 = 0; tmpChr2 < chromNum; tmpChr2 ++)
		{	
			//cout << "tmpChr2: " << tmpChr2 << endl;
			FusionJuncPosPairSupNumMap tmpFusionJuncPosPairSupNumMap;
			tmpFusionJuncMapVec.push_back(tmpFusionJuncPosPairSupNumMap);
		}
		fusionJuncMapVecVec.push_back(tmpFusionJuncMapVec);
	}

	cout << "start to insertUpdataFusionJuncMapVecVec with detected fusion junctions ...." << endl;
	//insertUpdataFusionJuncMapVecVec with detected Fusion junctions
	while(!fusionBreakPointResults_ifs.eof())
	{
		string tmpFusionStr;
		getline(fusionBreakPointResults_ifs, tmpFusionStr);
		//cout << "tmpFusionStr: " << tmpFusionStr << endl;
		if((fusionBreakPointResults_ifs.eof())||(tmpFusionStr == ""))
			break;		
		int detectedChrNameInt_1;
		int detectedChrNameInt_2;
		int detectedBreakPoint_1;
		int detectedBreakPoint_2;
		extractFusionInfoFromFusionBreakPointStr(
			detectedChrNameInt_1, detectedChrNameInt_2,
			detectedBreakPoint_1, detectedBreakPoint_2, indexInfo, tmpFusionStr);
		//cout << "detectedChrNameInt_1: " << detectedChrNameInt_1 << endl;
		//cout << "detectedChrNameInt_2: " << detectedChrNameInt_2 << endl;
		//cout << "detectedBreakPoint_1: " << detectedBreakPoint_1 << endl;
		//cout << "detectedBreakPoint_2: " << detectedBreakPoint_2 << endl;
		insertUpdataFusionJuncMapVecVecWithDetectedFusionJunc(
			fusionJuncMapVecVec, detectedChrNameInt_1, detectedChrNameInt_2,
			detectedBreakPoint_1, detectedBreakPoint_2);
	}
	cout << "start to output Fusion junction info ..." << endl;
	//output Fusion junction info
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		string tmpChrNameStr_1 = indexInfo->returnChrNameStr(tmpChr);
		for(int tmpChr2 = 0; tmpChr2 < chromNum; tmpChr2++)
		{
			string tmpChrNameStr_2 = indexInfo->returnChrNameStr(tmpChr2);
			for(FusionJuncPosPairSupNumMap::iterator tmp1stMapIter 
				= (fusionJuncMapVecVec[tmpChr])[tmpChr2].begin();
				tmp1stMapIter != (fusionJuncMapVecVec[tmpChr])[tmpChr2].end();
				tmp1stMapIter ++)
			{
				int tmpJunc_startPos = tmp1stMapIter->first;
				for(FusionJuncEndPos2supNumMap::iterator tmp2ndMapIter
					= (tmp1stMapIter->second).begin();
					tmp2ndMapIter != (tmp1stMapIter->second).end();
					tmp2ndMapIter ++)
				{
					int tmpJunc_endPos = tmp2ndMapIter->first;
					int tmpJunc_supportNum = tmp2ndMapIter->second;
					fusionJunc_ofs << tmpChrNameStr_1 << "\t" << tmpChrNameStr_2 << "\t"
						<< tmpJunc_startPos << "\t" << tmpJunc_endPos << "\t" 
						<< tmpJunc_supportNum << endl;
				}
			}
		}
	}

	delete indexInfo;
	fusionBreakPointResults_ifs.close();
	parameter_ifs.close();
	fusionJunc_ofs.close();
	return 0;
}