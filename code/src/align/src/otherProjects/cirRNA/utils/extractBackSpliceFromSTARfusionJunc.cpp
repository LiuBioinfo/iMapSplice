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

#include "../../../general/index_info.h"

using namespace std;

typedef map<int, int> SJendPos2supNumMap;
typedef map<int, SJendPos2supNumMap > SJchrPosPairMap;

void checkBackSpliceJunctionWithTargetSJchrPosPairMapVec(
	vector<SJchrPosPairMap>& targetSJchrPosPairMapVec,
	int tmpChrNameInt, int tmpChrPos_start, int tmpChrPos_end)
{
	SJchrPosPairMap::iterator tmpSJchrPosPairMapIter
		= targetSJchrPosPairMapVec[tmpChrNameInt].find(tmpChrPos_start);
	if(tmpSJchrPosPairMapIter 
		== targetSJchrPosPairMapVec[tmpChrNameInt].end()) // SJstartPos notFound, addSJ
	{
		SJendPos2supNumMap tmpSJendPos2supNumMap;
		tmpSJendPos2supNumMap.insert(pair<int,int>(tmpChrPos_end, 1));
		targetSJchrPosPairMapVec[tmpChrNameInt].insert(
			pair<int,SJendPos2supNumMap>(tmpChrPos_start, tmpSJendPos2supNumMap));
	}
	else // SJstartPos found, check SJendPos
	{
		SJendPos2supNumMap::iterator tmpSJendPos2supNumMapIter
			= (tmpSJchrPosPairMapIter->second).find(tmpChrPos_end);
		if(tmpSJendPos2supNumMapIter 
			== (tmpSJchrPosPairMapIter->second).end()) // SJendPos notFound, addSJ
			(tmpSJchrPosPairMapIter->second).insert(pair<int,int>(
				tmpChrPos_end, 1));
		else // SJendPos found, supportNum ++
			(tmpSJendPos2supNumMapIter->second) ++;
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath inputFusionSJsam_STARoutput backSpliceOutputFile" << endl;
		exit(1);
	}
	int backSpliceDistanceMax = 500000;

	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	string indexStr = indexFolderPath;
	indexStr += "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;
	vector<SJchrPosPairMap> spliceJunctionChrPosPairMapVec;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		SJchrPosPairMap tmpSJchrPosPairMap;
		spliceJunctionChrPosPairMapVec.push_back(tmpSJchrPosPairMap);
	}

	string inputFusionSJpath = argv[2];
	cout << "inputFusionSJpath: " << inputFusionSJpath << endl;
	string outputBackSplicePath = argv[3];
	cout << "outputBackSplicePath: " << outputBackSplicePath << endl;
	ifstream fusionSJ_ifs(inputFusionSJpath.c_str());
	ofstream backSplice_ofs(outputBackSplicePath.c_str());

	while(!fusionSJ_ifs.eof())
	{
		string tmpFusionSJstr;
		getline(fusionSJ_ifs, tmpFusionSJstr);
		if((fusionSJ_ifs.eof())||(tmpFusionSJstr == ""))
			break;
		//cout << "tmpFusionSJstr: " << tmpFusionSJstr << endl;
		vector<string> fusionSJfieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 7; tmp++)
		{
			int tabLoc = tmpFusionSJstr.find("\t", startLoc);
			string tmpFusionSJfieldStr = tmpFusionSJstr.substr(startLoc, tabLoc-startLoc);
			fusionSJfieldVec.push_back(tmpFusionSJfieldStr);
			startLoc = tabLoc + 1;
		}
		string tmpFusionSJ_startChrNameStr = fusionSJfieldVec[0];
		string tmpFusionSJ_startChrPosStr = fusionSJfieldVec[1];
		string tmpFusionSJ_1stStrandStr = fusionSJfieldVec[2];
		string tmpFusionSJ_endChrNameStr = fusionSJfieldVec[3];
		string tmpFusionSJ_endChrPosStr = fusionSJfieldVec[4];
		string tmpFusionSJ_2ndStrandStr = fusionSJfieldVec[5];
		string tmpFusionSJ_spanningOrEncompassingStr = fusionSJfieldVec[6];
		// cout << "tmpFusionSJ_startChrNameStr: " << tmpFusionSJ_startChrNameStr << endl;
		// cout << "tmpFusionSJ_startChrPosStr: " << tmpFusionSJ_startChrPosStr << endl;
		// cout << "tmpFusionSJ_endChrNameStr: " << tmpFusionSJ_endChrNameStr << endl;
		// cout << "tmpFusionSJ_endChrPosStr: " << tmpFusionSJ_endChrPosStr << endl;
		int tmpFusionSJ_startChrNameInt = indexInfo->convertStringToInt(
			tmpFusionSJ_startChrNameStr);
		int tmpFusionSJ_startChrPosInt = atoi(tmpFusionSJ_startChrPosStr.c_str()) - 1;
		int tmpFusionSJ_endChrNameInt = indexInfo->convertStringToInt(
			tmpFusionSJ_endChrNameStr);
		int tmpFusionSJ_endChrPosInt = atoi(tmpFusionSJ_endChrPosStr.c_str()) + 1;
		int tmpFusionSJ_spanningOrEncompassingCase = atoi(tmpFusionSJ_spanningOrEncompassingStr.c_str());
		//cout << "tmpFusionCase: " << tmpFusionSJ_spanningOrEncompassingCase << endl;
		if((tmpFusionSJ_startChrNameInt == tmpFusionSJ_endChrNameInt)
			&&(tmpFusionSJ_1stStrandStr == tmpFusionSJ_2ndStrandStr)
			&&(tmpFusionSJ_spanningOrEncompassingCase >= 0))
		{
			if((tmpFusionSJ_1stStrandStr == "+")
				&&(tmpFusionSJ_startChrPosInt > tmpFusionSJ_endChrPosInt)
				&&(tmpFusionSJ_startChrPosInt - tmpFusionSJ_endChrPosInt <= backSpliceDistanceMax))
			{
				checkBackSpliceJunctionWithTargetSJchrPosPairMapVec(spliceJunctionChrPosPairMapVec, 
					tmpFusionSJ_startChrNameInt, tmpFusionSJ_startChrPosInt, tmpFusionSJ_endChrPosInt);
			}
			else if((tmpFusionSJ_1stStrandStr == "-")
				&&(tmpFusionSJ_endChrPosInt > tmpFusionSJ_startChrPosInt)
				&&(tmpFusionSJ_endChrPosInt - tmpFusionSJ_startChrPosInt <= backSpliceDistanceMax))
			{
				checkBackSpliceJunctionWithTargetSJchrPosPairMapVec(spliceJunctionChrPosPairMapVec, 
					tmpFusionSJ_startChrNameInt, tmpFusionSJ_endChrPosInt, tmpFusionSJ_startChrPosInt);				
			}
			else
			{}
		}
	}

	int tmpJuncNum = 0;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++)
	{
		int tmpBackSpliceChrNameInt = tmpChr;
		string tmpBackSpliceChrNameStr = indexInfo->returnChrNameStr(tmpChr);
		for(SJchrPosPairMap::iterator tmpSJchrPosPairMapIter 
			= spliceJunctionChrPosPairMapVec[tmpChr].begin();
			tmpSJchrPosPairMapIter != spliceJunctionChrPosPairMapVec[tmpChr].end();
			tmpSJchrPosPairMapIter ++)
		{
			int tmpBackSpliceChrPos_start = tmpSJchrPosPairMapIter->first;
			for(SJendPos2supNumMap::iterator tmpSJendPos2supNumMapIter
				= (tmpSJchrPosPairMapIter->second).begin();
				tmpSJendPos2supNumMapIter != (tmpSJchrPosPairMapIter->second).end();
				tmpSJendPos2supNumMapIter ++)
			{
				tmpJuncNum ++;
				int tmpBackSpliceChrPos_end = tmpSJendPos2supNumMapIter->first;
				int tmpBackSpliceSupportNum = tmpSJendPos2supNumMapIter->second;
				backSplice_ofs << tmpBackSpliceChrNameStr << "\t"
					<< tmpBackSpliceChrPos_start << "\t"
					<< tmpBackSpliceChrPos_end << "\t"
					<< "JUNC_" << tmpJuncNum << "\t" << tmpBackSpliceSupportNum << endl;
			}
		}
	}

	delete indexInfo;
	fusionSJ_ifs.close();
	backSplice_ofs.close();
	return 0;
}