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
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

void extractFusionJuncInfoFromStr(
	int& tmpChrNameInt_1, int& tmpChrNameInt_2,
	int& tmpBreakPointPos_1, int& tmpBreakPointPos_2,
	//string& tmpStrand_1, string& tmpStrand_2, 
	//int& tmpAnchorSize_1, int& tmpAnchorSize_2,
	string& tmpFusionJuncStr, Index_Info* indexInfo)
{
	vector<string> tmpFusionJuncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 4; tmp++)
	{
		int tabLoc = tmpFusionJuncStr.find("\t", startLoc);
		string tmpFusionJuncField = tmpFusionJuncStr.substr(startLoc, tabLoc-startLoc);
		tmpFusionJuncFieldVec.push_back(tmpFusionJuncField);
		startLoc = tabLoc + 1;
	}
	tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpFusionJuncFieldVec[0]);
	tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpFusionJuncFieldVec[1]);	
	string tmpChrPosStr_1 = tmpFusionJuncFieldVec[2];
	string tmpChrPosStr_2 = tmpFusionJuncFieldVec[3];
	tmpBreakPointPos_1 = atoi(tmpChrPosStr_1.c_str());
	tmpBreakPointPos_2 = atoi(tmpChrPosStr_2.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolderPath inputGroundTruthSJ inputFusionSJ outputFilePrefix offset" << endl;
		exit(1);
	}
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	string offsetStr = argv[5];
	int offset = atoi(offsetStr.c_str());

	string outputFilePrefix = argv[4];
	string output_filterOut_path = outputFilePrefix + "_filterOut.txt";
	string output_kept_path = outputFilePrefix + "_kept.txt";
	ofstream filterOut_ofs(output_filterOut_path.c_str());
	ofstream kept_ofs(output_kept_path.c_str());

	string inputGroundTruthSJ = argv[2];
	AlignInferJunctionHash_Info* juncHash_groundTruth = new AlignInferJunctionHash_Info();
	juncHash_groundTruth->initiateAlignInferJunctionInfo(chromNum);
	juncHash_groundTruth->insertJuncFromJuncFile_chrNamePosOnly(inputGroundTruthSJ, indexInfo);

	string fusionJunc_path = argv[3];
	ifstream fusionJunc_ifs(fusionJunc_path.c_str());
	while(!fusionJunc_ifs.eof())
	{
		string tmpFusionJuncStr;
		getline(fusionJunc_ifs, tmpFusionJuncStr);
		//cout << "tmpFusionStr: " << tmpFusionStr << endl;
		if((fusionJunc_ifs.eof())||(tmpFusionJuncStr == ""))
			break;
		int detectedChrNameInt_1, detectedChrNameInt_2;
		int detectedBreakPoint_1, detectedBreakPoint_2;
		extractFusionJuncInfoFromStr(
			detectedChrNameInt_1, detectedChrNameInt_2,
			detectedBreakPoint_1, detectedBreakPoint_2,
			tmpFusionJuncStr, indexInfo);	
		if(detectedChrNameInt_1 != detectedChrNameInt_2)
		{	
			kept_ofs << tmpFusionJuncStr << endl;
			continue;
		}
		// int foundIndexInAlignInferJuncHash_1 
		// 	= juncHash_groundTruth->searchAndReturnAlignInferInfoVecIndex(
		// 		detectedChrNameInt_1, detectedBreakPoint_1, detectedBreakPoint_2);
		// int foundIndexInAlignInferJuncHash_2
		// 	= juncHash_groundTruth->searchAndReturnAlignInferInfoVecIndex(
		// 		detectedChrNameInt_1, detectedBreakPoint_2, detectedBreakPoint_1);
		bool tmpJuncFoundInGroundTruthSJhash = false;
		for(int tmp = 0-offset; tmp < offset; tmp++)
		{
			for(int tmp2 = 0-offset; tmp2 < offset; tmp2++)
			{
				int tmpBreakPointPos_1 = detectedBreakPoint_1 + tmp;
				int tmpBreakPointPos_2 = detectedBreakPoint_2 + tmp2;
				int tmpFoundIndex_1 = juncHash_groundTruth->searchAndReturnAlignInferInfoVecIndex(
					detectedChrNameInt_1, tmpBreakPointPos_1, tmpBreakPointPos_2);
				int tmpFoundIndex_2 = juncHash_groundTruth->searchAndReturnAlignInferInfoVecIndex(
					detectedChrNameInt_1, tmpBreakPointPos_2, tmpBreakPointPos_1);
				if((tmpFoundIndex_1 >= 0)||(tmpFoundIndex_2 >= 0))
				{
					tmpJuncFoundInGroundTruthSJhash = true;
					break;
				}
			}
			if(tmpJuncFoundInGroundTruthSJhash)
				break;
		}

		if(tmpJuncFoundInGroundTruthSJhash)
			filterOut_ofs << tmpFusionJuncStr << endl;
		else
			kept_ofs << tmpFusionJuncStr << endl;
	}
	fusionJunc_ifs.close();
	delete juncHash_groundTruth;
	filterOut_ofs.close();
	kept_ofs.close();
	return 0;
}