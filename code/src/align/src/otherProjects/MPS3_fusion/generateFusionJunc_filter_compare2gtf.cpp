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

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable inputIndexInfoPath inputGTFpath inputAdjustedFusionBreakPointFilePath ";
		cout << "inputSam2alignInferJuncHashPath offsetWhenComparingWithGTF ";
		cout << "checkStrandedOrNotWhenComparingWithGTF outputFolderPath" << endl;
		exit(1);
	}
	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... generateFusionJunc_filter_compare2gtf starts ......" << endl << endl;  
	//log_ofs << endl << "[" << asctime(local) << "... generateFusionJunc_filter_compare2gtf starts ......" << endl << endl; 

	string inputIndexInfoPath = argv[1];
	string inputGTFpath = argv[2];
	string inputAdjustedFusionBreakPointFilePath = argv[3];
	string inputSam2alignInferJuncHashPath = argv[4];
	string offsetWhenComparingWithGTF = argv[5];
	//int offset = atoi(offsetWhenComparingWithGTF.c_str());
	string checkStrandedOrNotWhenComparingWithGTF = argv[6];
	string outputFolderPath = argv[7];
	string fusionJuncPath = outputFolderPath + "/fusionJunc.txt";
	string fusionJuncPath_stranded = outputFolderPath + "/fusionJunc.stranded";
	string fusionJuncPath_nonStranded = outputFolderPath + "/fusionJunc.nonStranded";
	string fusionJuncPath_stranded_filter_folder = outputFolderPath + "/fusionJunc_stranded_filter";
	string fusionJuncPath_nonStranded_filter_folder = outputFolderPath + "/fusionJunc_nonStranded_filter";

	string mkdir= "mkdir -p " + outputFolderPath;
	system(mkdir.c_str());
	string log_path = outputFolderPath + "/log.txt";
	ofstream log_ofs(log_path.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	log_ofs << endl << "[" << asctime(local) << "... generateFusionJunc_filter_compare2gtf starts ......" << endl << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to do generateFusionJuncInfoFromAdjustedBreakPointFile !" << endl;
	log_ofs << endl << "[" << asctime(local) << "start to do generateFusionJuncInfoFromAdjustedBreakPointFile !" << endl;
	string cmd_generateFusionJunc = "./generateFusionJuncInfoFromAdjustedBreakPointFile "
		+ inputIndexInfoPath + " " + inputAdjustedFusionBreakPointFilePath + " "
		+ fusionJuncPath;
	system(cmd_generateFusionJunc.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to do separateFusionJuncStrandedOrNot !" << endl;
	log_ofs << endl << "[" << asctime(local) << "start to do separateFusionJuncStrandedOrNot !" << endl;
	string cmd_separateFusionJuncStrandedorNot = "./separateFusionJuncStrandedOrNot "
		+ fusionJuncPath + " " + outputFolderPath + "/fusionJunc";
	system(cmd_separateFusionJuncStrandedorNot.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to do filterFusionJunc_anchorSeqSimilarity !" << endl;
	log_ofs << endl << "[" << asctime(local) << "start to do filterFusionJunc_anchorSeqSimilarity !" << endl;
	string cmd_filterFusionJunc_stranded = "./filterFusionJunc_anchorSeqSimilarity "
		+ inputIndexInfoPath + " " + fusionJuncPath_stranded + " " 
		+ fusionJuncPath_stranded_filter_folder + " " + inputSam2alignInferJuncHashPath;
	system(cmd_filterFusionJunc_stranded.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to do filterNonStrandedFusionJunc_anchorSeqSimilarity !" << endl;
	log_ofs << endl << "[" << asctime(local) << "start to do filterNonStrandedFusionJunc_anchorSeqSimilarity !" << endl;
	string cmd_filterFusionJunc_nonStranded = "./filterNonStrandedFusionJunc_anchorSeqSimilarity "
		+ inputIndexInfoPath + " " + fusionJuncPath_nonStranded + " "
		+ fusionJuncPath_nonStranded_filter_folder + " " + inputSam2alignInferJuncHashPath;
	system(cmd_filterFusionJunc_nonStranded.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to cat Passed Fusion Junc !" << endl;
	log_ofs << endl << "[" << asctime(local) << "start to cat Passed Fusion Junc !" << endl;
	string cmd_catPassedFusionJunc = "cat " + fusionJuncPath_stranded_filter_folder 
		+ "/fusionJunc_classified_pass.txt " + fusionJuncPath_nonStranded_filter_folder
		+ "/fusionJunc_classified_pass.txt > " + outputFolderPath + "/fusionJunc_pass.txt";
	system(cmd_catPassedFusionJunc.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to cat compare2gtf !" << endl;
	log_ofs << endl << "[" << asctime(local) << "start to cat compare2gtf !" << endl;	
	string compare2gtf_folder = outputFolderPath + "/fusionJunc_pass_compare2gtf";
	string cmd_compare2gtf = "./compareFusionJuncResultsWithGTFannotation "
		+ inputIndexInfoPath + " " + inputGTFpath + " " + outputFolderPath + "/fusionJunc_pass.txt "
		+ compare2gtf_folder + " " + offsetWhenComparingWithGTF + " " + checkStrandedOrNotWhenComparingWithGTF;
	system(cmd_compare2gtf.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "end of generating fusion junc, junc filtering and comapre 2 gtf !" << endl;
	log_ofs << endl << "[" << asctime(local) << "end of generating fusion junc, junc filtering and comapre 2 gtf !" << endl;	

	log_ofs << endl << "input command: " << endl;
	cout << endl << "input command: " << endl;
	for(int tmp = 0; tmp < argc; tmp++)
	{
		log_ofs << "command " << tmp+1 << ": " << argv[tmp] << endl;
		cout << "command " << tmp+1 << ": " << argv[tmp] << endl;
	}	

	log_ofs << endl << "final fusion junc file: " << compare2gtf_folder << "/fusionJunc_geneInfo_interGene.txt" << endl;
	cout << endl << "final fusion junc file: " << compare2gtf_folder << "/fusionJunc_geneInfo_interGene.txt" << endl;
	log_ofs.close();
	return 0;
}