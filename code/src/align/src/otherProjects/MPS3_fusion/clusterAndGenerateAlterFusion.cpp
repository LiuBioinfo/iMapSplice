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
//#include "../../general/otherFunc2.h"
#include "general/fusionBreakPointHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 4)
	{
		cout << "Executable inputIndexFolder outputFolder inputFusionJuncPath_1 (inputFusionJuncPath_2 ...) " << endl;
		exit(1);
	}
	int offset = 100;
	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[2];
	string outputDirStr = outputFolderStr + "/";
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
   	string logStr = outputDirStr + "/log.txt";
   	ofstream log_ofs(logStr.c_str());
	log_ofs << "Command Line:";
	for(int tmp = 0; tmp < argc; tmp++)
		log_ofs << "\t" << argv[tmp];
   	string fusionBreakPointHashInfo_outputPath = outputDirStr + "/fusionBreakPointHash.txt";
   	string fusionBreakPointHashInfo_outputPath_withAlterFusionJunc = outputDirStr + "/fusionBreakPointHash_withAlterFusionJunc.txt";
   	string fusionBreakPointHashInfo_outputPath_withAlterFusionJunc_keptAF = outputDirStr + "/fusionBreakPointHash_withAlterFusionJunc_keptAF.txt";
   	string fusionBreakPointHashInfo_outputPath_withAlterFusionJunc_filterOutAF = outputDirStr + "/fusionBreakPointHash_withAlterFusionJunc_filterOutAF.txt";
	log_ofs << endl << "start to load whole genome" << endl;
	string globalIndexStr = argv[1];
	string indexStr = globalIndexStr + "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;
	vector<string> inputFusionJuncPathVec;
	for(int tmp = 3; tmp < argc; tmp++)
	{
		string tmpInputFusionJuncPath = argv[tmp];
		inputFusionJuncPathVec.push_back(tmpInputFusionJuncPath);
	}
	log_ofs << endl << "inputFusionJuncFile #: " << argc - 3 << endl;
	log_ofs << "start to generate fusionBreakPointHash ..." << endl;
	FusionBreakPointHash_Info* fusionBreakPointHashInfo = new FusionBreakPointHash_Info();
	fusionBreakPointHashInfo->initiateWithChromNum(chromNum);
	fusionBreakPointHashInfo->generateFusionBreakPointHashInfo_fromFuionJuncFileVec(inputFusionJuncPathVec, indexInfo);
	log_ofs << "start to output fusionBreakPointHashInfo ..." << endl;
	fusionBreakPointHashInfo->outputFusionBreakPointHashInfoStr(fusionBreakPointHashInfo_outputPath, indexInfo);
	log_ofs << "start to detect alterFusionJunc ..." << endl;
	fusionBreakPointHashInfo->detectAlterFusionJunc(offset);
	fusionBreakPointHashInfo->outputFusionBreakPointHashInfoStr_withAlterFusionJunc(fusionBreakPointHashInfo_outputPath_withAlterFusionJunc, indexInfo);
	fusionBreakPointHashInfo->memoryFree();
	delete fusionBreakPointHashInfo;
	delete indexInfo;

	log_ofs << "start to separate keptAF and filterOutAF fusions ..." << endl;
	
	string extractKeptAF_cmd = "grep kept_AF " + fusionBreakPointHashInfo_outputPath_withAlterFusionJunc 
		+ " > " + fusionBreakPointHashInfo_outputPath_withAlterFusionJunc_keptAF;
	system(extractKeptAF_cmd.c_str());
	
	string extractfilterOutAF_cmd = "grep filterOut_AF " + fusionBreakPointHashInfo_outputPath_withAlterFusionJunc 
		+ " > " + fusionBreakPointHashInfo_outputPath_withAlterFusionJunc_filterOutAF;	
	system(extractfilterOutAF_cmd.c_str());

	log_ofs << "all jobs done !" << endl;
	cout << "all jobs done !" << endl;
	return 0;
}