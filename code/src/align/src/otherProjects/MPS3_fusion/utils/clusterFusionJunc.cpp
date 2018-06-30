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
#include "../../../general/splice_info.h"

using namespace std;

void extractFusionJuncChrNamePosSupNumFromFusionJuncStr(
	int& tmpChrNameInt_1, int& tmpChrNameInt_2, 
	int& tmpChrPosInt_1, int& tmpChrPosInt_2, int& tmpSupNum,
	string& tmpFusionJuncStr, Index_Info* indexInfo)
{
	vector<string> tmpFusionJuncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 4; tmp++)
	{
		int tabLoc = tmpFusionJuncStr.find("\t", startLoc);
		tmpFusionJuncFieldVec.push_back(tmpFusionJuncStr.substr(startLoc, tabLoc - startLoc));
		startLoc = tabLoc + 1;
	}
	tmpFusionJuncFieldVec.push_back(tmpFusionJuncStr.substr(startLoc));

	string tmpChrNameStr_1 = tmpFusionJuncFieldVec[0];
	string tmpChrNameStr_2 = tmpFusionJuncFieldVec[1];
 	string tmpChrPosStr_1 = tmpFusionJuncFieldVec[2];
	string tmpChrPosStr_2 = tmpFusionJuncFieldVec[3];
	string tmpSupportNumStr = tmpFusionJuncFieldVec[4];

	tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1);
	tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
	tmpChrPosInt_1 = atoi(tmpChrPosStr_1.c_str());
	tmpChrPosInt_2 = atoi(tmpChrPosStr_2.c_str());
	tmpSupNum = atoi(tmpSupportNumStr.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolder inputRawFusionJuncPath outputClusteredFusionJuncPath" << endl;
		exit(1);
	}
	int clusterOffset = 100;

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "end of initiating indexInfo" << endl;

	string inputRawFusionJuncPath = argv[2];
	string outputClusteredFusionJuncPath = argv[3];
	ifstream rawFusionJunc_ifs(inputRawFusionJuncPath.c_str());
	ofstream clusteredFusionJunc_ofs(outputClusteredFusionJuncPath.c_str());

	vector < pair<int,int> > fusionJuncChrNameIntPairVec;
	vector < vector< pair<int,int> > > fusionJuncPosPairVecVec;
	vector < vector<int> > fusionJuncSupNumVecVec;
	while(!rawFusionJunc_ifs.eof())
	{
		string tmpRawFusionJuncStr;
		getline(rawFusionJunc_ifs, tmpRawFusionJuncStr);
		if((rawFusionJunc_ifs.eof())||(tmpRawFusionJuncStr == ""))
			break;
		int tmpChrNameInt_1, tmpChrNameInt_2, tmpChrPosInt_1, tmpChrPosInt_2, tmpSupNum;
		extractFusionJuncChrNamePosSupNumFromFusionJuncStr(
			tmpChrNameInt_1, tmpChrNameInt_2, 
			tmpChrPosInt_1, tmpChrPosInt_2, tmpSupNum,
			tmpRawFusionJuncStr, indexInfo);
		bool alreadyWithSomeClusterBool = false;
		for(int tmp = 0; tmp < fusionJuncChrNameIntPairVec.size(); tmp++)
		{
			int tmpExistingFusionJunc_chrNameInt_1 = (fusionJuncChrNameIntPairVec[tmp]).first;
			int tmpExistingFusionJunc_chrNameInt_2 = (fusionJuncChrNameIntPairVec[tmp]).second;
			int tmpExistingFusionJunc_pos_1 = (fusionJuncPosPairVecVec[tmp])[0].first;
			int tmpExistingFusionJunc_pos_2 = (fusionJuncPosPairVecVec[tmp])[0].second;
			if((tmpChrNameInt_1 == tmpExistingFusionJunc_chrNameInt_1)
				&&(tmpChrNameInt_2 == tmpExistingFusionJunc_chrNameInt_2)
				&&((tmpChrPosInt_1 >= tmpExistingFusionJunc_pos_1 - 100)
					&&(tmpChrPosInt_1 <= tmpExistingFusionJunc_pos_1 + 100))
				&&((tmpChrPosInt_2 >= tmpExistingFusionJunc_pos_2 - 100)
					&&(tmpChrPosInt_2 <= tmpExistingFusionJunc_pos_2 + 100)))
			{
				alreadyWithSomeClusterBool = true;
				(fusionJuncPosPairVecVec[tmp]).push_back(pair<int,int>(tmpChrPosInt_1, tmpChrPosInt_2));
				(fusionJuncSupNumVecVec[tmp]).push_back(tmpSupNum);
			}
		}
		if(!alreadyWithSomeClusterBool)
		{
			fusionJuncChrNameIntPairVec.push_back(pair<int,int>(tmpChrNameInt_1, tmpChrNameInt_2));
			vector< pair<int,int> > tmpFusionJuncPosPairVec;
			vector< int > tmpFusionSupNumVec;
			tmpFusionJuncPosPairVec.push_back(pair<int,int>(tmpChrPosInt_1, tmpChrPosInt_2));
			tmpFusionSupNumVec.push_back(tmpSupNum);
			fusionJuncPosPairVecVec.push_back(tmpFusionJuncPosPairVec);
			fusionJuncSupNumVecVec.push_back(tmpFusionSupNumVec);
		}
	}

	for(int tmp = 0; tmp < fusionJuncChrNameIntPairVec.size(); tmp++)
	{
		int candiJuncPointPairNum = fusionJuncPosPairVecVec[tmp].size();
		string tmpFusionChrName_1 = indexInfo->returnChrNameStr(fusionJuncChrNameIntPairVec[tmp].first);
		string tmpFusionChrName_2 = indexInfo->returnChrNameStr(fusionJuncChrNameIntPairVec[tmp].second);
		if(candiJuncPointPairNum == 1)
		{
			clusteredFusionJunc_ofs << tmpFusionChrName_1 << "\t" << tmpFusionChrName_2 << "\t"
				<< (fusionJuncPosPairVecVec[tmp])[0].first << "\t"
				<< (fusionJuncPosPairVecVec[tmp])[0].second << "\t" 
				<< (fusionJuncSupNumVecVec[tmp])[0] << endl;
		}
		else
		{
			// select the one with max support # for now.
			int tmpSupNum_max = 0;
			int tmpSupNum_max_indexInTmpCluster = 0;
			int tmpSupNum_sum = 0;
			int tmpCandiJuncPosPairNum = (fusionJuncSupNumVecVec[tmp]).size();
			for(int tmp2 = 0; tmp2 < tmpCandiJuncPosPairNum; tmp2 ++)
			{
				int tmpSupNum = (fusionJuncSupNumVecVec[tmp])[tmp2];
				tmpSupNum_sum += tmpSupNum;
				if(tmpSupNum > tmpSupNum_max)
				{
					tmpSupNum_max = tmpSupNum;
					tmpSupNum_max_indexInTmpCluster = tmp2;
				}
			}
			clusteredFusionJunc_ofs << tmpFusionChrName_1 << "\t" << tmpFusionChrName_2 << "\t"
				<< (fusionJuncPosPairVecVec[tmp])[tmpSupNum_max_indexInTmpCluster].first << "\t"
				<< (fusionJuncPosPairVecVec[tmp])[tmpSupNum_max_indexInTmpCluster].second << "\t" 
				<< tmpSupNum_sum;
				//<< (fusionJuncSupNumVecVec[tmp])[tmpSupNum_max_indexInTmpCluster];
			for(int tmp2 = 0; tmp2 < tmpCandiJuncPosPairNum; tmp2 ++)
			{
				int tmpPos_1 = (fusionJuncPosPairVecVec[tmp])[tmp2].first;
				int tmpPos_2 = (fusionJuncPosPairVecVec[tmp])[tmp2].second;
				int tmpSupNum = (fusionJuncSupNumVecVec[tmp])[tmp2];
				//int tmpOffset_1 = tmpPos_1 - (fusionJuncPosPairVecVec[tmp])[tmpSupNum_max_indexInTmpCluster].first;
				//int tmpOffset_2 = tmpPos_2 - (fusionJuncPosPairVecVec[tmp])[tmpSupNum_max_indexInTmpCluster].second;
				clusteredFusionJunc_ofs << "\t" << tmpPos_1 << "\t" << tmpPos_2 << "\t" << tmpSupNum;
			}
			clusteredFusionJunc_ofs << endl;
		}
	}

	rawFusionJunc_ifs.close();
	clusteredFusionJunc_ofs.close();
	return 0;
}