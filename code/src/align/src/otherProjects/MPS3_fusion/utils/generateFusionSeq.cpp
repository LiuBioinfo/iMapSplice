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

string returnFusionSeqInOneEnd(int chrNameInt,
	int breakPointPos, int anchorSize,
	bool gene1_or_gene2_bool, string& strand, 
	int toCheckAnchorLengthMax, Index_Info* indexInfo)
{
	int validAnchorSize = anchorSize;
	if(anchorSize > toCheckAnchorLengthMax)
		validAnchorSize = toCheckAnchorLengthMax;
	int fusionSeqInOneEnd_startPos; 
	int fusionSeqInOneEnd_endPos;
	bool for_or_rev_bool;
	if(gene1_or_gene2_bool) // gene 1 
	{
		if(strand == "+")
		{
			for_or_rev_bool = true;
			fusionSeqInOneEnd_startPos = breakPointPos - validAnchorSize + 1;
			fusionSeqInOneEnd_endPos = breakPointPos;
		}
		else
		{
			for_or_rev_bool = false;
			fusionSeqInOneEnd_startPos = breakPointPos;
			fusionSeqInOneEnd_endPos = breakPointPos + validAnchorSize - 1;
		}
	}
	else // gene2
	{
		if(strand == "+")
		{
			for_or_rev_bool = true;
			fusionSeqInOneEnd_startPos = breakPointPos;
			fusionSeqInOneEnd_endPos = breakPointPos + validAnchorSize - 1;
		}
		else
		{
			for_or_rev_bool = false;
			fusionSeqInOneEnd_startPos = breakPointPos - validAnchorSize + 1;
			fusionSeqInOneEnd_endPos = breakPointPos;
		}
	}
	if(for_or_rev_bool)
		return indexInfo->returnChromStrSubstr(chrNameInt, fusionSeqInOneEnd_startPos, 
			fusionSeqInOneEnd_endPos - fusionSeqInOneEnd_startPos + 1);
	else 
		return convertStringToReverseComplement(
				indexInfo->returnChromStrSubstr(
					chrNameInt, fusionSeqInOneEnd_startPos, 
					fusionSeqInOneEnd_endPos - fusionSeqInOneEnd_startPos + 1));
}

void extractFusionJuncInfoFromStr(
	int& tmpChrNameInt_1, int& tmpChrNameInt_2,
	int& tmpBreakPointPos_1, int& tmpBreakPointPos_2,
	string& tmpStrand_1, string& tmpStrand_2, 
	int& tmpAnchorSize_1, int& tmpAnchorSize_2,
	string& tmpFusionJuncStr, Index_Info* indexInfo)
{
	vector<string> tmpFusionJuncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 9; tmp++)
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
	tmpStrand_1 = tmpFusionJuncFieldVec[4];
	tmpStrand_2 = tmpFusionJuncFieldVec[5];
	tmpAnchorSize_1 = atoi(tmpFusionJuncFieldVec[7].c_str());
	tmpAnchorSize_2 = atoi(tmpFusionJuncFieldVec[8].c_str());
}

/*
string returnFusionSeq(int chrNameInt_1, int chrNameInt_2,
	int breakPointPos_1, int breakPointPos_2,
	int anchorSize_1, int anchorSize_2,
	string& strand_1, string& strand_2, 
	int toCheckAnchorLengthMax, Index_Info* indexInfo)
{
	string fusionSeq_gene1 = returnFusionSeqInOneEnd(
		chrNameInt_1, breakPointPos_1, anchorSize_1, strand_1, toCheckAnchorLengthMax, indexInfo);
	string fusionSeq_gene2 = returnFusionSeqInOneEnd(
		chrNameInt_2, breakPointPos_2, anchorSize_2, strand_2, toCheckAnchorLengthMax, indexInfo);
	return (fusionSeq_gene1 + fusionSeq_gene2);
}*/

int main(int argc, char** argv)
{
	if((argc != 4)&&(argc != 5))
	{
		cout << "Executable inputIndexInfoPath inputFusionJuncPath outputPath" << endl;
		exit(1);
	}
	int toCheckAnchorLengthMax = 30;
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
	indexInfo->readGenome(chrom);
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	string inputFusionJuncPath = argv[2];
	ifstream fusionJunc_ifs(inputFusionJuncPath.c_str());

	string outputPathStr = argv[3];
	ofstream fusionJunc_anchorSeq_ofs(outputPathStr.c_str());

	while(!fusionJunc_ifs.eof())
	{
		string tmpFusionJuncStr;
		getline(fusionJunc_ifs, tmpFusionJuncStr);
		//cout << "tmpFusionStr: " << tmpFusionStr << endl;
		if((fusionJunc_ifs.eof())||(tmpFusionJuncStr == ""))
			break;
		int detectedChrNameInt_1, detectedChrNameInt_2;
		int detectedBreakPoint_1, detectedBreakPoint_2;
		string detectedFusionJuncStrand_1, detectedFusionJuncStrand_2;
		int detectedFusionAnchorLength_1, detectedFusionAnchorLength_2;
		extractFusionJuncInfoFromStr(
			detectedChrNameInt_1, detectedChrNameInt_2,
			detectedBreakPoint_1, detectedBreakPoint_2,
			detectedFusionJuncStrand_1, detectedFusionJuncStrand_2,
			detectedFusionAnchorLength_1, detectedFusionAnchorLength_2,
			tmpFusionJuncStr, indexInfo);
		string fusionSeq_gene1 = returnFusionSeqInOneEnd(detectedChrNameInt_1,
			detectedBreakPoint_1, detectedFusionAnchorLength_1,
			true, detectedFusionJuncStrand_1, toCheckAnchorLengthMax, indexInfo);
		string fusionSeq_gene2 = returnFusionSeqInOneEnd(detectedChrNameInt_2,
			detectedBreakPoint_2, detectedFusionAnchorLength_2,
			false, detectedFusionJuncStrand_2, toCheckAnchorLengthMax, indexInfo);
		string fusionSeq = fusionSeq_gene1 + fusionSeq_gene2;
		fusionJunc_anchorSeq_ofs << tmpFusionJuncStr << endl
			<< fusionSeq_gene1 << endl << fusionSeq_gene2 << endl
			<< fusionSeq << endl;	
	}

	fusionJunc_anchorSeq_ofs.close();
	delete indexInfo;
	free(chrom);
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	return 0;
}