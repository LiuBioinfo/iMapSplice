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

#include "../../general/extractUnmapAlignment2ReadFile.h"
#include "../../phase1/arrayQueue.h"
#include "../../stats_info.h"
#include "../../constantDefinitions.h"
#include "../../general/option_info.h"
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/otherFunc.h"
#include "../../general/index_info.h"
#include "../../general/enhanced_suffix_array_info.h"
#include "../../general/annotation_info.h"
#include "../../phase1/repeatRegion.h"
#include "../../general/segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "../../general/splice_info.h"
#include "../../general/fixGapRelationParameters.h"
#include "../../general/read_info.h"
#include "../../general/seg_info.h"
//#include "general/fixDoubleAnchor_annotation_info.h"
#include "../../general/fixDoubleAnchorNWDP_info.h"
#include "../../general/fixDoubleAnchorMatch_info.h"
#include "../../general/fixDoubleAnchorInsertion_info.h"
#include "../../general/fixDoubleAnchorDeletion_info.h"
#include "../../general/fixDoubleAnchorSplice_complicate_info.h"
#include "../../general/fixDoubleAnchorSplice_info.h"
#include "../../general/fixDoubleAnchorCirRNA_info.h"
#include "../../general/path_info.h"
#include "../../general/gap_info.h"
#include "../../general/align_info.h"
#include "../../general/peAlign_info.h"
#include "../../general/groupSeg_info.h"
#include "../../general/alignInferJunctionHash_info_vec.h"
#include "../../phase2/spliceJunctionHash_info.h"
#include "../../phase2/unmapEnd_info.h"
#include "../../phase2/unfixedHead.h"
#include "../../phase2/unfixedTail.h"
#include "../../phase2/incompleteLongHead.h"
#include "../../phase2/incompleteLongTail.h"
#include "../../phase2/sam2junc.h"
#include "../../fixHeadTail.h"
#include "../../phase2/fixOneEndUnmapped.h"
#include "../../fixPhase1.h"
#include "../../general/readSeqPreProcessing.h"
#include "../../general/headerSection_info.h"
#include "../../general/otherFunc2.h"
#include "../../general/alignmentToJunc.h"
#include "general/uniqueUnpairedAlignment_info.h"
#include "general/fusionBreakPointHash_info.h"
#include "general/incompleteUniquePairedAlignment2detectFusion_info.h"

using namespace std;

int getMax(int a, int b)
{
	if(a <= b)
		return b;
	else
		return a;
}

void getFusionGeneFieldFromFusionJuncStrz_withGeneName(string& targetFJstr,
	string& tmpChrName_1_targetFJ, string& tmpChrName_2_targetFJ, 
	int& tmpChrPos_1_targetFJ, int& tmpChrPos_2_targetFJ,
	string& tmpStrand_1_targetFJ, string& tmpStrand_2_targetFJ, string& tmpFlankString_targetFJ,
	int& tmpAnchorSize_left_targetFJ, int& tmpAnchorSize_right_targetFJ, 
	int& tmpSupNum_targetFJ, int& tmpEncompNum_targetFJ,
	string& tmpGeneName_1_targetFJ, string& tmpGeneName_2_targetFJ)
{
	vector<string> tmpFusStrOriFieldVec;
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = targetFJstr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpFusField = targetFJstr.substr(startLoc, tabLoc-startLoc);
		tmpFusStrOriFieldVec.push_back(tmpFusField);
		startLoc = tabLoc + 1;
	}
	tmpFusStrOriFieldVec.push_back(targetFJstr.substr(startLoc));	
	int tmpFusStrOriFieldVecSize = tmpFusStrOriFieldVec.size();
	tmpChrName_1_targetFJ = tmpFusStrOriFieldVec[0];
	tmpChrName_2_targetFJ = tmpFusStrOriFieldVec[1];
	tmpChrPos_1_targetFJ = atoi(tmpFusStrOriFieldVec[2].c_str());
	tmpChrPos_2_targetFJ = atoi(tmpFusStrOriFieldVec[3].c_str());
	tmpStrand_1_targetFJ = tmpFusStrOriFieldVec[4];
	tmpStrand_2_targetFJ = tmpFusStrOriFieldVec[5];
	tmpFlankString_targetFJ = tmpFusStrOriFieldVec[6];
	tmpAnchorSize_left_targetFJ = atoi(tmpFusStrOriFieldVec[7].c_str());
	tmpAnchorSize_right_targetFJ = atoi(tmpFusStrOriFieldVec[8].c_str());
	tmpSupNum_targetFJ = atoi(tmpFusStrOriFieldVec[9].c_str());
	tmpEncompNum_targetFJ = atoi(tmpFusStrOriFieldVec[10].c_str());
	tmpGeneName_1_targetFJ = tmpFusStrOriFieldVec[tmpFusStrOriFieldVecSize-2];
	tmpGeneName_2_targetFJ = tmpFusStrOriFieldVec[tmpFusStrOriFieldVecSize-1];
}

void getFusionGeneFieldFromFusionJuncStr(string& targetFJstr,
	string& tmpChrName_1_targetFJ, string& tmpChrName_2_targetFJ, 
	int& tmpChrPos_1_targetFJ, int& tmpChrPos_2_targetFJ,
	string& tmpStrand_1_targetFJ, string& tmpStrand_2_targetFJ, string& tmpFlankString_targetFJ,
	int& tmpAnchorSize_left_targetFJ, int& tmpAnchorSize_right_targetFJ, 
	int& tmpSupNum_targetFJ, int& tmpEncompNum_targetFJ)
{
	vector<string> tmpFusStrOriFieldVec;
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = targetFJstr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpFusField = targetFJstr.substr(startLoc, tabLoc-startLoc);
		tmpFusStrOriFieldVec.push_back(tmpFusField);
		startLoc = tabLoc + 1;
	}
	tmpFusStrOriFieldVec.push_back(targetFJstr.substr(startLoc));	
	//int tmpFusStrOriFieldVecSize = tmpFusStrOriFieldVec.size();
	tmpChrName_1_targetFJ = tmpFusStrOriFieldVec[0];
	tmpChrName_2_targetFJ = tmpFusStrOriFieldVec[1];
	tmpChrPos_1_targetFJ = atoi(tmpFusStrOriFieldVec[2].c_str());
	tmpChrPos_2_targetFJ = atoi(tmpFusStrOriFieldVec[3].c_str());
	tmpStrand_1_targetFJ = tmpFusStrOriFieldVec[4];
	tmpStrand_2_targetFJ = tmpFusStrOriFieldVec[5];
	tmpFlankString_targetFJ = tmpFusStrOriFieldVec[6];
	tmpAnchorSize_left_targetFJ = atoi(tmpFusStrOriFieldVec[7].c_str());
	tmpAnchorSize_right_targetFJ = atoi(tmpFusStrOriFieldVec[8].c_str());
	tmpSupNum_targetFJ = atoi(tmpFusStrOriFieldVec[9].c_str());
	tmpEncompNum_targetFJ = atoi(tmpFusStrOriFieldVec[10].c_str());
	//tmpGeneName_1_targetFJ = tmpFusStrOriFieldVec[tmpFusStrOriFieldVecSize-2];
	//tmpGeneName_2_targetFJ = tmpFusStrOriFieldVec[tmpFusStrOriFieldVecSize-1];
}

void mergeFusionJuncResults(
	string& inputTargetFusionJuncPath, 
	string& remap_fusionJuncPath, 
	string& countCompletePair_fusionJuncPath, 
	string& countIncompletePair_fusionJuncPath, 
	string& mergedFusionJuncResults)
{
	ifstream targetFJ_ifs(inputTargetFusionJuncPath.c_str());
	ifstream remapFJ_ifs(remap_fusionJuncPath.c_str());
	ifstream countCompletePairFJ_ifs(countCompletePair_fusionJuncPath.c_str());
	ifstream countIncompletePairFJ_ifs(countIncompletePair_fusionJuncPath.c_str());
	ofstream mergedFJ_ofs(mergedFusionJuncResults.c_str());
	while(!targetFJ_ifs.eof())
	{
		string targetFJstr, remapFJstr, countCompletePairFJstr, countIncompletePairFJstr;
		getline(targetFJ_ifs, targetFJstr);
		if(targetFJstr == "")
			break;
		getline(remapFJ_ifs, remapFJstr);
		getline(countCompletePairFJ_ifs, countCompletePairFJstr);
		getline(countIncompletePairFJ_ifs, countIncompletePairFJstr);

		string tmpChrName_1_targetFJ, tmpChrName_2_targetFJ; 
		int tmpChrPos_1_targetFJ, tmpChrPos_2_targetFJ;
		string tmpStrand_1_targetFJ, tmpStrand_2_targetFJ, tmpFlankString_targetFJ;
		int tmpAnchorSize_left_targetFJ, tmpAnchorSize_right_targetFJ, tmpSupNum_targetFJ, tmpEncompNum_targetFJ;
		string tmpGeneName_1_targetFJ, tmpGeneName_2_targetFJ;
		getFusionGeneFieldFromFusionJuncStrz_withGeneName(targetFJstr,
			tmpChrName_1_targetFJ, tmpChrName_2_targetFJ, tmpChrPos_1_targetFJ, tmpChrPos_2_targetFJ,
			tmpStrand_1_targetFJ, tmpStrand_2_targetFJ, tmpFlankString_targetFJ,
			tmpAnchorSize_left_targetFJ, tmpAnchorSize_right_targetFJ, tmpSupNum_targetFJ, tmpEncompNum_targetFJ,
			tmpGeneName_1_targetFJ, tmpGeneName_2_targetFJ);

		string tmpChrName_1_remapFJ, tmpChrName_2_remapFJ; 
		int tmpChrPos_1_remapFJ, tmpChrPos_2_remapFJ;
		string tmpStrand_1_remapFJ, tmpStrand_2_remapFJ, tmpFlankString_remapFJ;
		int tmpAnchorSize_left_remapFJ, tmpAnchorSize_right_remapFJ, tmpSupNum_remapFJ, tmpEncompNum_remapFJ;
		getFusionGeneFieldFromFusionJuncStr(remapFJstr,
			tmpChrName_1_remapFJ, tmpChrName_2_remapFJ, tmpChrPos_1_remapFJ, tmpChrPos_2_remapFJ,
			tmpStrand_1_remapFJ, tmpStrand_2_remapFJ, tmpFlankString_remapFJ,
			tmpAnchorSize_left_remapFJ, tmpAnchorSize_right_remapFJ, tmpSupNum_remapFJ, tmpEncompNum_remapFJ);

		string tmpChrName_1_countCompletePairFJ, tmpChrName_2_countCompletePairFJ; 
		int tmpChrPos_1_countCompletePairFJ, tmpChrPos_2_countCompletePairFJ;
		string tmpStrand_1_countCompletePairFJ, tmpStrand_2_countCompletePairFJ, tmpFlankString_countCompletePairFJ;
		int tmpAnchorSize_left_countCompletePairFJ, tmpAnchorSize_right_countCompletePairFJ, tmpSupNum_countCompletePairFJ, tmpEncompNum_countCompletePairFJ;
		getFusionGeneFieldFromFusionJuncStr(countCompletePairFJstr,
			tmpChrName_1_countCompletePairFJ, tmpChrName_2_countCompletePairFJ, tmpChrPos_1_countCompletePairFJ, tmpChrPos_2_countCompletePairFJ,
			tmpStrand_1_countCompletePairFJ, tmpStrand_2_countCompletePairFJ, tmpFlankString_countCompletePairFJ,
			tmpAnchorSize_left_countCompletePairFJ, tmpAnchorSize_right_countCompletePairFJ, tmpSupNum_countCompletePairFJ, tmpEncompNum_countCompletePairFJ);

		string tmpChrName_1_countIncompletePairFJ, tmpChrName_2_countIncompletePairFJ; 
		int tmpChrPos_1_countIncompletePairFJ, tmpChrPos_2_countIncompletePairFJ;
		string tmpStrand_1_countIncompletePairFJ, tmpStrand_2_countIncompletePairFJ, tmpFlankString_countIncompletePairFJ;
		int tmpAnchorSize_left_countIncompletePairFJ, tmpAnchorSize_right_countIncompletePairFJ, tmpSupNum_countIncompletePairFJ, tmpEncompNum_countIncompletePairFJ;
		getFusionGeneFieldFromFusionJuncStr(countIncompletePairFJstr,
			tmpChrName_1_countIncompletePairFJ, tmpChrName_2_countIncompletePairFJ, tmpChrPos_1_countIncompletePairFJ, tmpChrPos_2_countIncompletePairFJ,
			tmpStrand_1_countIncompletePairFJ, tmpStrand_2_countIncompletePairFJ, tmpFlankString_countIncompletePairFJ,
			tmpAnchorSize_left_countIncompletePairFJ, tmpAnchorSize_right_countIncompletePairFJ, tmpSupNum_countIncompletePairFJ, tmpEncompNum_countIncompletePairFJ);
		
		int tmpUpdatedAnchorSize_left = getMax(tmpAnchorSize_left_targetFJ, tmpAnchorSize_left_remapFJ);
		int tmpUpdatedAnchorSize_right = getMax(tmpAnchorSize_right_targetFJ, tmpAnchorSize_right_remapFJ);

		mergedFJ_ofs << tmpChrName_1_targetFJ << "\t" << tmpChrName_2_targetFJ << "\t" 
			<< tmpChrPos_1_targetFJ << "\t" << tmpChrPos_2_targetFJ << "\t" 
			<< tmpStrand_1_targetFJ << "\t" << tmpStrand_2_targetFJ << "\t" << tmpFlankString_targetFJ << "\t"
			<< tmpUpdatedAnchorSize_left << "\t" << tmpUpdatedAnchorSize_right << "\t"
			<< tmpSupNum_targetFJ << "\t" << tmpSupNum_remapFJ - tmpSupNum_targetFJ << "\t"
			<< tmpEncompNum_countCompletePairFJ << "\t" << tmpEncompNum_countIncompletePairFJ << "\t" 
			<< tmpSupNum_remapFJ << "\t" << tmpEncompNum_countCompletePairFJ + tmpEncompNum_countIncompletePairFJ << "\t"
			<< tmpSupNum_remapFJ + tmpEncompNum_countCompletePairFJ + tmpEncompNum_countIncompletePairFJ << "\t"
			<< tmpGeneName_1_targetFJ << "\t" << tmpGeneName_2_targetFJ << endl;
	}

	targetFJ_ifs.close();
	remapFJ_ifs.close();
	countCompletePairFJ_ifs.close();
	countIncompletePairFJ_ifs.close();
	mergedFJ_ofs.close();
	return;
}

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cout << "Executable inputIndexFolder inputTargetFusionJuncPath ";
		cout << "inputInitialMapPhase2outputFolder inputNonFusionSamInGlobalMapOutput ";
		cout << "outputFolderPath threads_num fasta_or_fastq strandedOrNot " << endl;
		exit(1);
	}

	string inputIndexFolderPath = argv[1];
	//string inputBinFolderPath = argv[2]
	string inputTargetFusionJuncPath = argv[2];
	string inputInitialMapPhase2outputFolderPath = argv[3];
	string inputInitialMap_completeUnpair = inputInitialMapPhase2outputFolderPath
		+ "/fixHeadTail_complete_unpair.sam";
	string inputInitialMap_incompleteUnpair = inputInitialMapPhase2outputFolderPath
		+ "/fixHeadTail_incomplete_unpair.sam";		
	string inputNonFusionSamInGlobalMapOutputPath = argv[4];
	string outputFolderPath = argv[5];
	string threads_num_str = argv[6];
	string fasta_or_fastq_str = argv[7];
	string strandedOrNotStr = argv[8];

	string mkdir= "mkdir -p " + outputFolderPath;
	system(mkdir.c_str());
	string remapOutputFolderPath = outputFolderPath + "/remapOutput";
	string countCompletePairOutputPath = outputFolderPath + "/countCompletePairOutput";
	string countIncompletePairOutputPath = outputFolderPath + "/countIncompletePairOutput";

	cout << endl << "start to do remap" << endl;
	string cmd_remap = "./remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint "
		+ inputIndexFolderPath + " " + inputTargetFusionJuncPath + " "
		+ inputNonFusionSamInGlobalMapOutputPath + " " + remapOutputFolderPath + " "
		+ threads_num_str + " " + fasta_or_fastq_str;
	system(cmd_remap.c_str());

	cout << endl << "start to do completePair count" << endl;
	string cmd_completePairCount = "./countUniqueUnpairedReadsEncompassingFusionBreakPoint "
		+ inputIndexFolderPath + " " + inputTargetFusionJuncPath + " "
		+ inputInitialMap_completeUnpair + " " + countCompletePairOutputPath + " " 
		+ threads_num_str + " " + fasta_or_fastq_str + " " + strandedOrNotStr;
	system(cmd_completePairCount.c_str());

	cout << endl << "start to do incompletePair count" << endl;
	string cmd_incompletePairCount = "./countUniqueUnpairedReadsEncompassingFusionBreakPoint "
		+ inputIndexFolderPath + " " + inputTargetFusionJuncPath + " "
		+ inputInitialMap_incompleteUnpair + " " + countIncompletePairOutputPath + " "
		+ threads_num_str + " " + fasta_or_fastq_str + " " + strandedOrNotStr;
	system(cmd_incompletePairCount.c_str());

	cout << "start to merge fusion detection results" << endl;
	string remap_fusionJuncPath = remapOutputFolderPath + "/fusionBreakPointHashInfo_updated.txt";
	string countCompletePair_fusionJuncPath = countCompletePairOutputPath + "/fusionBreakPointHashInfo_updated.txt";
	string countIncompletePair_fusionJuncPath = countIncompletePairOutputPath + "/fusionBreakPointHashInfo_updated.txt";
	string mergedFusionJuncResults = outputFolderPath + "/mergedFusionJuncResults.txt";
	mergeFusionJuncResults(inputTargetFusionJuncPath, remap_fusionJuncPath, countCompletePair_fusionJuncPath, 
		countIncompletePair_fusionJuncPath, mergedFusionJuncResults);
	cout << "end of running remapAndCountEncompRead !" << endl;
	return 0;
}