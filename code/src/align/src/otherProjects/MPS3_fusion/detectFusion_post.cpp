// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//Note: this program is used to stitch up all the sub-programs for fusion detection
// Note: no header section in input SAM file
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
#include "general/mapSplice_fusion_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputIndexFolderPath inputInitialMapResultsFolderPath outputFolderPath threads_num fasta_or_fastq min_fusion_distance" << endl;
		exit(1);
	}

	MapSplice_Fusion_Info mapSpliceFusionInfo;

	time_t nowtime;
	struct tm *local;

	// get command line info
	string inputIndexFolderPath = argv[1];
	string inputInitialMapResultsFolderPath = argv[2];
	string outputFolderPath = argv[3];
	string threads_num_str = argv[4];
	string fasta_or_fastq_str = argv[5];
	string min_fusion_distance_str = argv[6];

	string threadNumStr = argv[4];
	int threads_num = atoi(threadNumStr.c_str());
	omp_set_num_threads(threads_num);

	int min_fusion_distance = atoi(min_fusion_distance_str.c_str());

	// creat output folder, and log file
	cout << "creating results folder ...." << endl;
	string mkdir= "mkdir -p " + outputFolderPath;
	system(mkdir.c_str());
   	string logStr = outputFolderPath + "/log.txt";
   	string settingsLogStr = outputFolderPath + "/settings.txt";
   	ofstream log_ofs(logStr.c_str());
   	ofstream settings_ofs(settingsLogStr.c_str());

	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... MPS-fusion starts ......" << endl << endl;  
	log_ofs << endl << "[" << asctime(local) << "... MPS-fusion starts ......" << endl << endl;  
	log_ofs << "Command Line: " << endl;
	for(int tmp = 0; tmp < argc; tmp++)
		log_ofs << "\t" << argv[tmp];
	log_ofs << endl << "************************" << endl;

	bool fasta_or_fastq_bool = true;
	if((fasta_or_fastq_str == "fasta")||(fasta_or_fastq_str == "Fasta"))
		fasta_or_fastq_bool = true;
	else if((fasta_or_fastq_str == "fastq")||(fasta_or_fastq_str == "Fastq"))
		fasta_or_fastq_bool = false;
	else
	{
		cout << "Please set the correct format, fasta or fastq" << endl;
		exit(1);
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 	

	//string index = inputIndexFolderPath + "/";
    string preIndexArrayPreStr = inputIndexFolderPath; 	
    
	string preIndexMapLengthArrayStr = preIndexArrayPreStr + "_MapLength";
	ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr = preIndexArrayPreStr + "_IntervalStart"; 
	ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr = preIndexArrayPreStr + "_IntervalEnd";
	ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);
	int* preIndexMapLengthArray; 
	preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int)); 
	preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray; 
	preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); 
	preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray; 
	preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); 
	preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));
	cout << "finish loading preIndex ..." << endl;
 	log_ofs << "finish loading preIndex ..." << endl;
 	string indexStr = inputIndexFolderPath + "/";
	string chrom_bit_file = indexStr + "_chrom"; 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr + "_parameter";
	ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	log_ofs << "index: " << indexStr << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	log_ofs << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	log_ofs << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;

	string SA_file = indexStr; SA_file.append("_SA"); 
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); 
	string childTab_file = indexStr; childTab_file.append("_childTab"); 
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); 	
    unsigned int *sa; sa = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	unsigned int *childTab; childTab = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	BYTE *lcpCompress; lcpCompress = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 
	BYTE *verifyChild; verifyChild = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 	

	ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl;

	log_ofs << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl;
	log_ofs << "start to load SA" << endl;
	SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load lcpCompress" << endl;
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_ofs << "start to load childTab " << endl;
	childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load detChild" << endl;
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	////////////   global_map_outer_soft_clip_paired_alignment_to_detect_fusion starts ...//////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to do global mapping ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "start to do global mapping ......" << endl << endl;

	int offsetForAdjustingFusionBreakPoint_max = 5;
	mapSpliceFusionInfo.initiateInputOutputFilePath_globalMapPairedAlignment();
	mapSpliceFusionInfo.globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion(
		indexInfo, sa, lcpCompress, childTab, chrom, verifyChild, preIndexMapLengthArray, 
		preIndexIntervalStartArray, preIndexIntervalEndArray, //repeatRegionInfo,
		offsetForAdjustingFusionBreakPoint_max, threads_num, fasta_or_fastq_bool, 
		min_fusion_distance, log_ofs);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to generateFusionJuncFromAdjustedFusionBreakPointFile ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "start to generateFusionJuncFromAdjustedFusionBreakPointFile ......" << endl << endl;

	mapSpliceFusionInfo.initiateInputOutputFilePath_generateFusionjuncFromAdjustedFusionBreakPointFile();
	mapSpliceFusionInfo.generateFusionJuncFromAdjustedFusionBreakPointFile(
		threads_num, fasta_or_fastq_bool, indexInfo, log_ofs);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to do remapping against fusion break points ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "start to do remapping against fusion break points ......" << endl << endl;	

	mapSpliceFusionInfo.initiateInputOutputFilePath_remapPairedAlignment();
	mapSpliceFusionInfo.remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint(
		threads_num, fasta_or_fastq_bool, indexInfo, log_ofs);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "start to countUniqueUnpairedReadsEncompassingFusionBreakPoint ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "start to countUniqueUnpairedReadsEncompassingFusionBreakPoint ......" << endl << endl;		
	mapSpliceFusionInfo.initiateInputOutputFilePath_countUniqueUnpairedReadsEncompassingFusionBreakPoint();
	mapSpliceFusionInfo.countEncompassingReads(threads_num, fasta_or_fastq_bool, indexInfo, log_ofs);

	free(preIndexMapLengthArray);
	free(preIndexIntervalStartArray);
	free(preIndexIntervalEndArray);
	free(sa); 
	free(childTab);
	free(lcpCompress);
	free(verifyChild);
	SA_file_ifs.close();
	lcpCompress_file_ifs.close();
	childTab_file_ifs.close();
	verifyChild_file_ifs.close();
	preIndexMapLengthArray_ifs.close();
	preIndexIntervalStartArray_ifs.close();
	preIndexIntervalEndArray_ifs.close();
	parameter_file_ifs.close();
	chrom_bit_file_ifs.close();
	log_ofs.close();
	settings_ofs.close();
	return 0;
}