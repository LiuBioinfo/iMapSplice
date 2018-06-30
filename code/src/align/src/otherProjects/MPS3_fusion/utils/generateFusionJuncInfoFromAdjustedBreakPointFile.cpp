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

#include "../../../general/extractUnmapAlignment2ReadFile.h"
#include "../../../phase1/arrayQueue.h"
#include "../../../stats_info.h"
#include "../../../constantDefinitions.h"
#include "../../../general/option_info.h"
#include "../../../general/read_block_test.h"
#include "../../../general/bwtmap_info.h"
#include "../../../general/DoubleAnchorScore.h"
#include "../../../general/sbndm.h"
#include "../../../general/otherFunc.h"
#include "../../../general/index_info.h"
#include "../../../general/enhanced_suffix_array_info.h"
#include "../../../general/annotation_info.h"
#include "../../../phase1/repeatRegion.h"
#include "../../../general/segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "../../../general/splice_info.h"
#include "../../../general/fixGapRelationParameters.h"
#include "../../../general/read_info.h"
#include "../../../general/seg_info.h"
//#include "general/fixDoubleAnchor_annotation_info.h"
#include "../../../general/fixDoubleAnchorNWDP_info.h"
#include "../../../general/fixDoubleAnchorMatch_info.h"
#include "../../../general/fixDoubleAnchorInsertion_info.h"
#include "../../../general/fixDoubleAnchorDeletion_info.h"
#include "../../../general/fixDoubleAnchorSplice_complicate_info.h"
#include "../../../general/fixDoubleAnchorSplice_info.h"
#include "../../../general/fixDoubleAnchorCirRNA_info.h"
#include "../../../general/path_info.h"
#include "../../../general/gap_info.h"
#include "../../../general/align_info.h"
#include "../../../general/peAlign_info.h"
#include "../../../general/groupSeg_info.h"
#include "../../../general/alignInferJunctionHash_info_vec.h"
#include "../../../phase2/spliceJunctionHash_info.h"
#include "../../../phase2/unmapEnd_info.h"
#include "../../../phase2/unfixedHead.h"
#include "../../../phase2/unfixedTail.h"
#include "../../../phase2/incompleteLongHead.h"
#include "../../../phase2/incompleteLongTail.h"
#include "../../../phase2/sam2junc.h"
#include "../../../fixHeadTail.h"
#include "../../../phase2/fixOneEndUnmapped.h"
#include "../../../fixPhase1.h"
#include "../../../general/readSeqPreProcessing.h"
#include "../../../general/headerSection_info.h"
#include "../../../general/otherFunc2.h"
#include "../../../general/alignmentToJunc.h"

using namespace std;

typedef map<int, int> FusionJuncEndPos2fusionInfoVecIndexMap;
typedef map<int, FusionJuncEndPos2fusionInfoVecIndexMap > FusionJuncPosPair2fusionInfoVecIndexMap;

void extractFusionInfoFromAdjustedFusionBreakPointStr(
	int& detectedChrNameInt_1, int& detectedChrNameInt_2,
	int& detectedBreakPoint_1, int& detectedBreakPoint_2,
	string& detectedFusionJuncStrand_1, string& detectedFusionJuncStrand_2,
	string& detectedFusionJuncFlankString,
	int& detectedAnchorLength_1, int& detectedAnchorLength_2,
	int& XM, //string& nonStrandedFusionJuncCaseStr,
	string& fusionJuncCaseStr, int& mapRangeMax_1, int& mapRangeMax_2,
	string& tmpFusionBreakPointInfoStr, Index_Info* indexInfo)
{
	vector<string> fusionJuncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 13; tmp++)
	{
		int tabLoc = tmpFusionBreakPointInfoStr.find("\t", startLoc);
		string tmpDetectedFusionJuncField = tmpFusionBreakPointInfoStr.substr(startLoc, tabLoc-startLoc);
		fusionJuncFieldVec.push_back(tmpDetectedFusionJuncField);
		startLoc = tabLoc + 1;
	}
	fusionJuncFieldVec.push_back(tmpFusionBreakPointInfoStr.substr(startLoc));
	// int possiLastTabLoc = tmpFusionBreakPointInfoStr.find("\t", startLoc);
	// if(possiLastTabLoc == string::npos) // stranded, nonStrandedFusionJuncCaseStr is not available,
	// 	fusionJuncFieldVec.push_back(tmpFusionBreakPointInfoStr.substr(startLoc)); // push back XMstr
	// else // nonStranded, nonStrandedFusionJuncCaseStr is available,
	// {
	// 	fusionJuncFieldVec.push_back(tmpFusionBreakPointInfoStr.substr(startLoc, possiLastTabLoc-startLoc)); // push back XMstr
	// 	startLoc = possiLastTabLoc + 1;
	// 	fusionJuncFieldVec.push_back(tmpFusionBreakPointInfoStr.substr(startLoc));
	// }

	string tmpChrNameStr_1 = fusionJuncFieldVec[1];
	string tmpChrNameStr_2 = fusionJuncFieldVec[2];
	string tmpChrPosStr_1 = fusionJuncFieldVec[3];
	string tmpChrPosStr_2 = fusionJuncFieldVec[4];
	detectedChrNameInt_1 = indexInfo->convertStringToInt(fusionJuncFieldVec[1]);
	detectedChrNameInt_2 = indexInfo->convertStringToInt(fusionJuncFieldVec[2]);
	detectedBreakPoint_1 = atoi(tmpChrPosStr_1.c_str());
	detectedBreakPoint_2 = atoi(tmpChrPosStr_2.c_str());
	detectedFusionJuncStrand_1 = fusionJuncFieldVec[5];
	detectedFusionJuncStrand_2 = fusionJuncFieldVec[6];
	detectedFusionJuncFlankString = fusionJuncFieldVec[7];
	string tmpAnchorLengthStr_1 = fusionJuncFieldVec[8];
	string tmpAnchorLengthStr_2 = fusionJuncFieldVec[9];
	detectedAnchorLength_1 = atoi(tmpAnchorLengthStr_1.c_str());
	detectedAnchorLength_2 = atoi(tmpAnchorLengthStr_2.c_str());
	string XMstr = fusionJuncFieldVec[10];
	XM = atoi(XMstr.c_str());
	//if(fusionJuncFieldVec.size() == 12)
	//	nonStrandedFusionJuncCaseStr = fusionJuncFieldVec[11];
	fusionJuncCaseStr = fusionJuncFieldVec[11];
	string mapRangeMax_1_str = fusionJuncFieldVec[12];
	string mapRangeMax_2_str = fusionJuncFieldVec[13];
	mapRangeMax_1 = atoi(mapRangeMax_1_str.c_str());
	mapRangeMax_2 = atoi(mapRangeMax_2_str.c_str());
}

void updateFusionJuncMapRange(string& tmpFusionCaseStr, int existingMapRangeMax_1, int existingMapRangeMax_2,
	int tmpMapRangeMax_1, int tmpMapRangeMax_2, int& updatedMapRangeMax_1, int& updatedMapRangeMax_2)
{
	// update mapRangeMax_1
	if((tmpFusionCaseStr == "1,")||(tmpFusionCaseStr == "2,")||(tmpFusionCaseStr == "7,")||(tmpFusionCaseStr == "8,")
		||(tmpFusionCaseStr == "1,4,")||(tmpFusionCaseStr == "2,5,")||(tmpFusionCaseStr == "7,8,"))
	{
		if(tmpMapRangeMax_1 < existingMapRangeMax_1)
			updatedMapRangeMax_1 = tmpMapRangeMax_1;
		else
			updatedMapRangeMax_1 = existingMapRangeMax_1;
	}
	else if((tmpFusionCaseStr == "4,")||(tmpFusionCaseStr == "5,")||(tmpFusionCaseStr == "10,")||(tmpFusionCaseStr == "11,")
		||(tmpFusionCaseStr == "10,11,"))
	{
		if(tmpMapRangeMax_1 > existingMapRangeMax_1)
			updatedMapRangeMax_1 = tmpMapRangeMax_1;
		else
			updatedMapRangeMax_1 = existingMapRangeMax_1;
	}
	else
	{
		cout << "error in updateFusionJuncMapRange, tmpFusionCaseStr: " << tmpFusionCaseStr << endl;
		exit(1);
	}
	// update mapRangeMax_2
	if((tmpFusionCaseStr == "1,")||(tmpFusionCaseStr == "2,")||(tmpFusionCaseStr == "10,")||(tmpFusionCaseStr == "11,")
		||(tmpFusionCaseStr == "1,4,")||(tmpFusionCaseStr == "2,5,")||(tmpFusionCaseStr == "10,11,"))
	{
		if(tmpMapRangeMax_2 > existingMapRangeMax_2)
			updatedMapRangeMax_2 = tmpMapRangeMax_2;
		else
			updatedMapRangeMax_2 = existingMapRangeMax_2;		
	}
	else if((tmpFusionCaseStr == "4,")||(tmpFusionCaseStr == "5,")||(tmpFusionCaseStr == "7,")||(tmpFusionCaseStr == "8,")
		||(tmpFusionCaseStr == "7,8,"))
	{
		if(tmpMapRangeMax_2 < existingMapRangeMax_2)
			updatedMapRangeMax_2 = tmpMapRangeMax_2;
		else
			updatedMapRangeMax_2 = existingMapRangeMax_2;
	}
	else
	{
		cout << "error in updateFusionJuncMapRange, tmpFusionCaseStr: " << tmpFusionCaseStr << endl;
		exit(1);
	}
}

void insertUpdataFusionJuncMapVecVecAndFusionInfoVecWithDetectedFusionJunc(
	vector < vector< FusionJuncPosPair2fusionInfoVecIndexMap > >& fusionJuncMapVecVec, 
	int detectedChrNameInt_1, int detectedChrNameInt_2,
	int detectedBreakPoint_1, int detectedBreakPoint_2,
	string& detectedFusionJuncStrand_1, string& detectedFusionJuncStrand_2,
	string& detectedFusionJuncFlankString,
	int detectedAnchorLength_1, int detectedAnchorLength_2,
	int XM, 
	string& fusionJuncCaseStr,
	int mapRangeMax_1, int mapRangeMax_2,
	//string& nonStrandedFusionJuncCaseStr,
	vector<int>& fusionJuncSupNumVec,
	vector<string>& fusionJuncStrandVec_1,
	vector<string>& fusionJuncStrandVec_2,
	vector<string>& fusionJuncFlankStringVec, 
	vector<int>& fusionJuncAnchorLengthVec_1,
	vector<int>& fusionJuncAnchorLengthVec_2,
	vector<int>& fusionJuncXMvec_min,
	vector<int>& fusionJuncXMvec_max,
	vector<string>& fusionJuncCaseStrVec,
	vector<int>& fusionJuncMapRangeMaxVec_1,
	vector<int>& fusionJuncMapRangeMaxVec_2,
	//vector<string>& nonStrandedFusionJuncCaseStrVec,
	Index_Info* indexInfo)
{
	FusionJuncPosPair2fusionInfoVecIndexMap::iterator tmp1stMapIter 
		= ((fusionJuncMapVecVec[detectedChrNameInt_1])[detectedChrNameInt_2]).find(detectedBreakPoint_1);
	if(tmp1stMapIter == ((fusionJuncMapVecVec[detectedChrNameInt_1])[detectedChrNameInt_2]).end())
	{
		int tmpFusionInfoJuncSize = fusionJuncSupNumVec.size();
		FusionJuncEndPos2fusionInfoVecIndexMap tmp2ndMap;
		tmp2ndMap.insert(pair<int,int>(detectedBreakPoint_2, tmpFusionInfoJuncSize));
		((fusionJuncMapVecVec[detectedChrNameInt_1])[detectedChrNameInt_2]).insert(
			pair<int,FusionJuncEndPos2fusionInfoVecIndexMap>(detectedBreakPoint_1, tmp2ndMap));
		fusionJuncSupNumVec.push_back(1);
		fusionJuncStrandVec_1.push_back(detectedFusionJuncStrand_1);
		fusionJuncStrandVec_2.push_back(detectedFusionJuncStrand_2);
		fusionJuncFlankStringVec.push_back(detectedFusionJuncFlankString);
		fusionJuncAnchorLengthVec_1.push_back(detectedAnchorLength_1);
		fusionJuncAnchorLengthVec_2.push_back(detectedAnchorLength_2);
		fusionJuncXMvec_min.push_back(XM);
		fusionJuncXMvec_max.push_back(XM);
		//nonStrandedFusionJuncCaseStrVec.push_back(nonStrandedFusionJuncCaseStr);
		fusionJuncCaseStrVec.push_back(fusionJuncCaseStr);
		fusionJuncMapRangeMaxVec_1.push_back(mapRangeMax_1);
		fusionJuncMapRangeMaxVec_2.push_back(mapRangeMax_2);
	}
	else
	{
		FusionJuncEndPos2fusionInfoVecIndexMap::iterator tmp2ndMapIter
			= (tmp1stMapIter->second).find(detectedBreakPoint_2);
		if(tmp2ndMapIter == (tmp1stMapIter->second).end())
		{	
			int tmpFusionInfoJuncSize = fusionJuncSupNumVec.size();
			(tmp1stMapIter->second).insert(pair<int,int>(detectedBreakPoint_2, tmpFusionInfoJuncSize));
			fusionJuncSupNumVec.push_back(1);
			fusionJuncStrandVec_1.push_back(detectedFusionJuncStrand_1);
			fusionJuncStrandVec_2.push_back(detectedFusionJuncStrand_2);
			fusionJuncFlankStringVec.push_back(detectedFusionJuncFlankString);			
			fusionJuncAnchorLengthVec_1.push_back(detectedAnchorLength_1);
			fusionJuncAnchorLengthVec_2.push_back(detectedAnchorLength_2);
			fusionJuncXMvec_min.push_back(XM);
			fusionJuncXMvec_max.push_back(XM);
			//nonStrandedFusionJuncCaseStrVec.push_back(nonStrandedFusionJuncCaseStr);
			fusionJuncCaseStrVec.push_back(fusionJuncCaseStr);
			fusionJuncMapRangeMaxVec_1.push_back(mapRangeMax_1);
			fusionJuncMapRangeMaxVec_2.push_back(mapRangeMax_2);
		}
		else
		{
			int tmpIndexInFusionInfoVec = (tmp2ndMapIter->second);
			fusionJuncSupNumVec[tmpIndexInFusionInfoVec]++;
			int existingAnchorLength_1 = fusionJuncAnchorLengthVec_1[tmpIndexInFusionInfoVec];
			int existingAnchorLength_2 = fusionJuncAnchorLengthVec_2[tmpIndexInFusionInfoVec];
			int tmpXM_min = fusionJuncXMvec_min[tmpIndexInFusionInfoVec];
			int tmpXM_max = fusionJuncXMvec_max[tmpIndexInFusionInfoVec];
			if(detectedAnchorLength_1 > existingAnchorLength_1)
			{
				fusionJuncAnchorLengthVec_1[tmpIndexInFusionInfoVec] = detectedAnchorLength_1;
			}
			if(detectedAnchorLength_2 > existingAnchorLength_2)
			{
				fusionJuncAnchorLengthVec_2[tmpIndexInFusionInfoVec] = detectedAnchorLength_2;
			}
			if(XM < tmpXM_min)
				fusionJuncXMvec_min[tmpIndexInFusionInfoVec] = XM;
			if(XM > tmpXM_max)
				fusionJuncXMvec_max[tmpIndexInFusionInfoVec] = XM;
			string tmpFusionCaseStr = fusionJuncCaseStrVec[tmpIndexInFusionInfoVec];
			int existingMapRangeMax_1 = fusionJuncMapRangeMaxVec_1[tmpIndexInFusionInfoVec];
			int existingMapRangeMax_2 = fusionJuncMapRangeMaxVec_2[tmpIndexInFusionInfoVec];
			int updatedFusionJuncMapRangeMax_1, updatedFusionJuncMapRangeMax_2;
			updateFusionJuncMapRange(tmpFusionCaseStr, existingMapRangeMax_1, existingMapRangeMax_2,
				mapRangeMax_1, mapRangeMax_2, updatedFusionJuncMapRangeMax_1, updatedFusionJuncMapRangeMax_2);
			fusionJuncMapRangeMaxVec_1[tmpIndexInFusionInfoVec] = updatedFusionJuncMapRangeMax_1;
			fusionJuncMapRangeMaxVec_2[tmpIndexInFusionInfoVec] = updatedFusionJuncMapRangeMax_2;
		}
	}	
}


/*
bool checkMatchTroughSeqSimilarity(string& strand_1, string& strand_2, 
	int chrNameInt_1, int chrNameInt_2, int breakPointPos_1, int breakPointPos_2,
	int anchorLength_1, int anchorLength_2, int toCheckAnchorLengthMax, 
	Index_Info* indexInfo, int& validAnchorLength_1, int& validAnchorLength_2,
	int& penalty_1, int& penalty_2)
{
	if((strand_1 == "N") || (strand_2 == "N"))
		return false;
	validAnchorLength_1 = anchorLength_1;
	validAnchorLength_2 = anchorLength_2;
	if(validAnchorLength_1 > toCheckAnchorLengthMax)
		validAnchorLength_1 = toCheckAnchorLengthMax;
	if(validAnchorLength_2 > toCheckAnchorLengthMax)
		validAnchorLength_2 = toCheckAnchorLengthMax;

	string matchThroughSeq_gene1;
	string anchorSeq_gene2;

	string matchThroughSeq_gene2;
	string anchorSeq_gene1;
	// cout << "validAnchorLength_1: " << validAnchorLength_1 << endl;
	// cout << "validAnchorLength_2: " << validAnchorLength_2 << endl;
	// cout << "strand_1: " << strand_1 << endl;
	// cout << "strand_2: " << strand_2 << endl;
	if((strand_1 == "+")&&(strand_2 == "+"))
	{
		matchThroughSeq_gene1 = indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1+1, validAnchorLength_2);
		anchorSeq_gene2 = indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2, validAnchorLength_2);
		matchThroughSeq_gene2 = indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 - 1 - validAnchorLength_1 + 1, validAnchorLength_1);
		anchorSeq_gene1 = indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1 - validAnchorLength_1 + 1, validAnchorLength_1);
	}
	else if((strand_1 == "-")&&(strand_2 == "-"))
	{
		matchThroughSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1 - 1 - validAnchorLength_2 + 1, validAnchorLength_2));
		anchorSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 - validAnchorLength_2 + 1, validAnchorLength_2));
		matchThroughSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 + 1, validAnchorLength_1));
		anchorSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1, validAnchorLength_1));
	}
	else if((strand_1 == "+")&&(strand_2 == "-"))
	{
		matchThroughSeq_gene1 = indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1+1, validAnchorLength_2);
		anchorSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 - validAnchorLength_2 + 1, validAnchorLength_2));
		matchThroughSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 + 1, validAnchorLength_1));
		anchorSeq_gene1 = indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1 - validAnchorLength_1 + 1, validAnchorLength_1);
	}
	else
	{
		matchThroughSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1 - 1 - validAnchorLength_2 + 1, validAnchorLength_2));
		anchorSeq_gene2 = indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2, validAnchorLength_2);
		matchThroughSeq_gene2 = indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 - 1 - validAnchorLength_1 + 1, validAnchorLength_1);
		anchorSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1, validAnchorLength_1));
	}
	//cout << "matchThroughSeq_gene1: " << matchThroughSeq_gene1 << endl;
	//cout << "matchThroughSeq_gene2: " << matchThroughSeq_gene2 << endl;
	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_1;
	tmpFixNWDPinfo_1.doNWDP_withMismatchJumpCode(
		matchThroughSeq_gene1, anchorSeq_gene2);
	penalty_1 = tmpFixNWDPinfo_1.getPenalty();
	//delete tmpFixNWDPinfo_1;
	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_2;
	tmpFixNWDPinfo_2.doNWDP_withMismatchJumpCode(
		matchThroughSeq_gene2, anchorSeq_gene1);
	penalty_2 = tmpFixNWDPinfo_2.getPenalty();
	//delete tmpFixNWDPinfo_2;
	return true;
}*/

int main(int argc, char** argv)
{
	if((argc != 4)&&(argc != 5))
	{
		cout << "Executable inputIndexInfoPath inputAdjustedFusionBreakPointFilePath outputFusionJuncPath" << endl;
		//cout << " (inputNormalAlignInferJuncHashWithOutAnchorSizePath)" << endl;
		exit(1);
	}
	int toCheckAnchorLengthMax = 30;

	bool checkAlterSpliceBool;
	if(argc == 4)
		checkAlterSpliceBool = false;
	else
		checkAlterSpliceBool = true;
	/*int toCheckAlterSpliceAnchorLengthMax = 10;
	string inputNormalAlignInferJuncHashPath;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();	
	SJhash_Info* SJhashInfo = new SJhash_Info();
	// insert SJ into SJmap
	if(checkAlterSpliceBool)
	{
		inputNormalAlignInferJuncHashPath = argv[4];
		int chromNum = indexInfo->returnChromNum();
		cout << "start to initaite alignInferJunctionHashInfo " << endl;
		alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);		
		cout << "start to initaite SJhashInfo" << endl;
		SJhashInfo->initiateAreaAndStringHash(chromNum);
		cout << "start to insert SJ into alignInferJunctionHashInfo" << endl;
		alignInferJunctionHashInfo->insertJuncFromJuncFile(inputNormalAlignInferJuncHashPath, indexInfo);
		cout << "start to convert 2 SJhashInfo" << endl;
		alignInferJunctionHashInfo->convert2SJhashInfo(SJhashInfo, indexInfo);
	}*/

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	//settings_log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	//settings_log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;
	cout << "end of initiating indexInfo" << endl;

	string inputFuionReadWithBreakPointInfoFilePath = argv[2];
	ifstream fusionBreakPointResults_ifs(inputFuionReadWithBreakPointInfoFilePath.c_str());

	string outputPath = argv[3];
	ofstream fusionJunc_ofs(outputPath.c_str());
	cout << "chromNum: " << chromNum << endl;

	vector< vector< FusionJuncPosPair2fusionInfoVecIndexMap > > fusionJuncMapVecVec;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++ )
	{
		vector<FusionJuncPosPair2fusionInfoVecIndexMap> tmpFusionJuncMapVec;
		for(int tmpChr2 = 0; tmpChr2 < chromNum; tmpChr2 ++)
		{	
			FusionJuncPosPair2fusionInfoVecIndexMap tmpFusionJuncPosPair2fusionInfoVecIndexMap;
			tmpFusionJuncMapVec.push_back(tmpFusionJuncPosPair2fusionInfoVecIndexMap);
		}
		fusionJuncMapVecVec.push_back(tmpFusionJuncMapVec);
	}

	vector<int> fusionJuncSupNumVec;
	vector<string> fusionJuncStrandVec_1;
	vector<string> fusionJuncStrandVec_2;
	vector<string> fusionJuncFlankStringVec;
	vector<int> fusionJuncAnchorLengthVec_1;
	vector<int> fusionJuncAnchorLengthVec_2;
	vector<int> fusionJuncXMvec_min;
	vector<int> fusionJuncXMvec_max;
	vector<string> fusionJuncCaseStrVec;
	vector<int> fusionJuncMapRangeMaxVec_1;
	vector<int> fusionJuncMapRangeMaxVec_2;

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
		string detectedFusionJuncStrand_1;
		string detectedFusionJuncStrand_2;
		string detectedFusionJuncFlankString;
		int detectedFusionAnchorLength_1;
		int detectedFusionAnchorLength_2;
		int XM;
		//string nonStrandedFusionJuncCaseStr = "NULL";
		string detectedFusionCaseStr;
		int detectedFusionMapRangeMax_1;
		int detectedFusionMapRangeMax_2;

		extractFusionInfoFromAdjustedFusionBreakPointStr(
			detectedChrNameInt_1, detectedChrNameInt_2,
			detectedBreakPoint_1, detectedBreakPoint_2, 
			detectedFusionJuncStrand_1, detectedFusionJuncStrand_2,
			detectedFusionJuncFlankString,
			detectedFusionAnchorLength_1, detectedFusionAnchorLength_2,
			XM, detectedFusionCaseStr, 
			detectedFusionMapRangeMax_1, detectedFusionMapRangeMax_2,
			tmpFusionStr, indexInfo);
		//cout << "detectedChrNameInt_1: " << detectedChrNameInt_1 << endl;
		//cout << "detectedChrNameInt_2: " << detectedChrNameInt_2 << endl;
		//cout << "detectedBreakPoint_1: " << detectedBreakPoint_1 << endl;
		//cout << "detectedBreakPoint_2: " << detectedBreakPoint_2 << endl;
		insertUpdataFusionJuncMapVecVecAndFusionInfoVecWithDetectedFusionJunc(
			fusionJuncMapVecVec, 
			detectedChrNameInt_1, detectedChrNameInt_2,
			detectedBreakPoint_1, detectedBreakPoint_2,
			detectedFusionJuncStrand_1, detectedFusionJuncStrand_2,
			detectedFusionJuncFlankString,
			detectedFusionAnchorLength_1, detectedFusionAnchorLength_2,
			XM, detectedFusionCaseStr,
			detectedFusionMapRangeMax_1, detectedFusionMapRangeMax_2,
			fusionJuncSupNumVec,
			fusionJuncStrandVec_1, fusionJuncStrandVec_2,
			fusionJuncFlankStringVec, 
			fusionJuncAnchorLengthVec_1, fusionJuncAnchorLengthVec_2,
			fusionJuncXMvec_min, fusionJuncXMvec_max,
			fusionJuncCaseStrVec,
			fusionJuncMapRangeMaxVec_1, fusionJuncMapRangeMaxVec_2,
			indexInfo);
	}


	cout << "start to output Fusion junction info ..." << endl;
	//output Fusion junction info
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		string tmpChrNameStr_1 = indexInfo->returnChrNameStr(tmpChr);
		for(int tmpChr2 = 0; tmpChr2 < chromNum; tmpChr2++)
		{
			string tmpChrNameStr_2 = indexInfo->returnChrNameStr(tmpChr2);
			for(FusionJuncPosPair2fusionInfoVecIndexMap::iterator tmp1stMapIter 
				= (fusionJuncMapVecVec[tmpChr])[tmpChr2].begin();
				tmp1stMapIter != (fusionJuncMapVecVec[tmpChr])[tmpChr2].end();
				tmp1stMapIter ++)
			{
				int tmpJunc_startPos = tmp1stMapIter->first;
				for(FusionJuncEndPos2fusionInfoVecIndexMap::iterator tmp2ndMapIter
					= (tmp1stMapIter->second).begin();
					tmp2ndMapIter != (tmp1stMapIter->second).end();
					tmp2ndMapIter ++)
				{
					int tmpJunc_endPos = tmp2ndMapIter->first;
					int tmpJunc_index = tmp2ndMapIter->second;
					int tmpJunc_supNum = fusionJuncSupNumVec[tmpJunc_index];
					string tmpJunc_strand_1 = fusionJuncStrandVec_1[tmpJunc_index];
					string tmpJunc_strand_2 = fusionJuncStrandVec_2[tmpJunc_index];
					string tmpJunc_flankString = fusionJuncFlankStringVec[tmpJunc_index];
					int tmpJunc_anchorLength_1 = fusionJuncAnchorLengthVec_1[tmpJunc_index];
					int tmpJunc_anchorLength_2 = fusionJuncAnchorLengthVec_2[tmpJunc_index];
					int tmpJunc_XM_min = fusionJuncXMvec_min[tmpJunc_index];
					int tmpJunc_XM_max = fusionJuncXMvec_max[tmpJunc_index];
					string tmpJunc_fusionCase = fusionJuncCaseStrVec[tmpJunc_index];
					int tmpJunc_mapRangeMax_1 = fusionJuncMapRangeMaxVec_1[tmpJunc_index];
					int tmpJunc_mapRangeMax_2 = fusionJuncMapRangeMaxVec_2[tmpJunc_index];
					int tmpJunc_supNum_encompassing = 0;
					fusionJunc_ofs << tmpChrNameStr_1 << "\t" << tmpChrNameStr_2 << "\t"
						<< tmpJunc_startPos << "\t" << tmpJunc_endPos << "\t" 
						<< tmpJunc_strand_1 << "\t" << tmpJunc_strand_2 << "\t" 
						<< tmpJunc_flankString << "\t" 
						<< tmpJunc_anchorLength_1 << "\t" << tmpJunc_anchorLength_2 << "\t"
						<< tmpJunc_supNum << "\t" << int_to_str(tmpJunc_supNum_encompassing) << "\t"
						<< tmpJunc_XM_min << "\t" << int_to_str(tmpJunc_XM_max);
					// if((tmpJunc_strand_1 == "N")&&(tmpJunc_strand_2 == "N"))	
					// 	fusionJunc_ofs << "\t" << nonStrandedFusionJuncCaseStrVec[tmpJunc_index] << endl;
					// else
					// 	fusionJunc_ofs << endl;
					fusionJunc_ofs << "\t" << tmpJunc_fusionCase << "\t" << tmpJunc_mapRangeMax_1 << "\t" << tmpJunc_mapRangeMax_2 << endl;
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