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

#include "general/extractUnmapAlignment2ReadFile.h"
#include "phase1/arrayQueue.h"
#include "stats_info.h"
#include "constantDefinitions.h"
#include "general/option_info.h"
#include "general/read_block_test.h"
#include "general/bwtmap_info.h"
#include "general/DoubleAnchorScore.h"
#include "general/sbndm.h"
#include "general/otherFunc.h"
#include "general/index_info.h"
#include "general/enhanced_suffix_array_info.h"
#include "general/annotation_info.h"
#include "phase1/repeatRegion.h"
#include "general/segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "general/splice_info.h"
#include "general/fixGapRelationParameters.h"
#include "general/read_info.h"
#include "general/seg_info.h"
//#include "general/fixDoubleAnchor_annotation_info.h"
#include "general/fixDoubleAnchorNWDP_info.h"
#include "general/fixDoubleAnchorMatch_info.h"
#include "general/fixDoubleAnchorInsertion_info.h"
#include "general/fixDoubleAnchorDeletion_info.h"
#include "general/fixDoubleAnchorSplice_complicate_info.h"
#include "general/fixDoubleAnchorSplice_info.h"
#include "general/fixDoubleAnchorCirRNA_info.h"
#include "general/path_info.h"
#include "general/gap_info.h"
#include "otherProjects/incorporateGenomicVariants/general/syntheticSNPtransSeq_info.h"
#include "general/align_info.h"
#include "general/peAlign_info.h"
#include "general/groupSeg_info.h"
#include "general/alignInferJunctionHash_info_vec.h"
#include "phase2/spliceJunctionHash_info.h"
#include "phase2/unmapEnd_info.h"
#include "phase2/unfixedHead.h"
#include "phase2/unfixedTail.h"
#include "phase2/incompleteLongHead.h"
#include "phase2/incompleteLongTail.h"
#include "phase2/sam2junc.h"
#include "fixHeadTail.h"
#include "phase2/fixOneEndUnmapped.h"
#include "fixPhase1.h"
#include "general/readSeqPreProcessing.h"
#include "general/headerSection_info.h"
#include "general/otherFunc2.h"
//#include "general/alignmentToJunc.h"
//#define PreIndexSize 268435456

using namespace std;  


#ifdef CAL_TIME
clock_t input_begin, input_end, input_cost = 0,
	output_begin, output_end, output_cost = 0,
	process_begin, process_end, process_cost = 0,
	
	getReadInfo_begin, getReadInfo_end, getReadInfo_cost = 0, 
	getSegInfo_begin, getSegInfo_end, getSegInfo_cost = 0,
	updateSegInfo_begin, updateSegInfo_end, updateSegInfo_cost = 0,
	getPathInfo_begin, getPathInfo_end, getPathInfo_cost = 0,
	getGapInfo_begin, getGapInfo_end, getGapInfo_cost  = 0,
	getAlignInfo_begin, getAlignInfo_end, getAlignInfo_cost = 0,
	selectAlignInfo_begin, selectAlignInfo_end, selectAlignInfo_cost = 0,
	printAlignInfo_begin, printAlignInfo_end, printAlignInfo_cost = 0,

	updateSegInfo_greedyMap_begin, updateSegInfo_greedyMap_end, updateSegInfo_greedyMap_cost = 0,
	updateSegInfo_targetMap_begin, updateSegInfo_targetMap_end, updateSegInfo_targetMap_cost = 0;
#endif

//unsigned int PairedReadNum = 0, BothUnmappedReadNum = 0, BothUnmappedReadNum_mappedToRepeatRegion = 0, UnpairedReadNum = 0;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void mmapcopy(int fd, int size, void* dest, ofstream& log_ofs)
{
	log_ofs << "start to mmapCopy" << endl;
	void* start_addr = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
	log_ofs << errno << endl;
	log_ofs << "finish mmap" << endl;
	if (start_addr == (void *)-1)
	{
		log_ofs << "mmap failed" << endl;
		fprintf(stderr, "mmap: %s\n", strerror(errno));
	}
	//cout << "before memory copy " << endl;
	memcpy(dest, start_addr, size);
	munmap(start_addr, size);
	return;
}


int main(int argc, char**argv)
{
    Option_Info* optionInfo = new Option_Info();
    optionInfo->getOpt_long(argc, argv);
    string outputDirStr = optionInfo->outputFolder_path; //argv[3];
   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
	string mkdirOutputCommand_phase1 = "mkdir -p " + outputDirStr + "/phase1_output";
	system(mkdirOutputCommand_phase1.c_str());
	string mkdirOutputCommand_repeatRegionFile = mkdirOutputCommand_phase1 + "/repeat_region";
   	system(mkdirOutputCommand_repeatRegionFile.c_str());
	string mkdirOutputCommand_tmpAlignCompleteRead = mkdirOutputCommand_phase1 + "/completePair";
	system(mkdirOutputCommand_tmpAlignCompleteRead.c_str());
	string mkdirOutputCommand_tmpAlignOneEndUnmapped = mkdirOutputCommand_phase1 + "/oneEndUnmapped";
	system(mkdirOutputCommand_tmpAlignOneEndUnmapped.c_str());	
	string mkdirOutputCommand_phase2 = "mkdir -p " + outputDirStr + "/phase2_output";
	system(mkdirOutputCommand_phase2.c_str());

	int SNPlocInSyntheticSNPseq = 101;
    bool checkQualSeqForShortAnchorSeqToTargetMap = false;//true;
    cout << "Attention! checkQualSeqForShortAnchorSeqToTargetMap true or not: " 
    	<< checkQualSeqForShortAnchorSeqToTargetMap << endl;
    bool checkQualSeqForReadSegSeq = false;//true;
	cout << "Attention! checkQualSeqForReadSegSeq true or not: " 
		<< checkQualSeqForReadSegSeq << endl;
	bool outputUnpairedSAM_bool = false;
	bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool = false;
	//outputUnpairedSAM_bool = true;

	#ifdef OUTPUT_UNPAIRED_REGULAR
	outputUnpairedSAM_bool = true;
	outputUnpairedSAM_bothEndsUniqueMappedOnly_bool = false;
	#endif	

	#ifdef OUTPUTUNIQUEUNPAIREDBOTHENDSMAPPEDSAM
	outputUnpairedSAM_bool = true;
	outputUnpairedSAM_bothEndsUniqueMappedOnly_bool = true;
	#endif 

	#ifdef DETECT_CIRCULAR_RNA
	if(outputUnpairedSAM_bool)
	{
		cout << "for now, MPS3_circRNA only supports SE reads." << endl;
		cout << "MPS3_circRNA does not support 'outputUnpairedSAM' " << endl;
		exit(1);
	}
	#endif 

    /////////////////   get option from command line ////////////////////
    bool extractUnmapAlignment2ReadFile_bool 
    	= optionInfo->return_extractUnmapAlignment2ReadFile_bool();
    string inputSamFilePathStr = optionInfo->returnInputSamFilePath();
    if(extractUnmapAlignment2ReadFile_bool)
    {
    	unmappedAlignemnt2ReadFile(inputSamFilePathStr);
    }

    ////////////////////////////// check SE or PE reads  ///////////////////////////////
    bool SE_or_PE_bool = optionInfo->returnSEorPE_bool();
    
	//////////////////////////////////////////////////
   	string settingsLogStr = outputDirStr + "/settings.log";
   	ofstream settings_log_ofs(settingsLogStr.c_str());
   	string progressLogStr = outputDirStr + "/process.log";
   	ofstream log_ofs(progressLogStr.c_str());

   	string runtimeLogStr = outputDirStr + "/runtime.log";
   	ofstream runtime_log_ofs(runtimeLogStr.c_str());
   	string statsStr = outputDirStr + "/stats.txt";
   	ofstream stats_ofs(statsStr.c_str());

   	optionInfo->outputOptStr(settings_log_ofs);

	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   switches of seperate processes    ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
   	bool annotation_provided_bool = optionInfo->annotation_provided_bool;
   	bool Do_annotation_only_bool = false;//annotation_provided_bool;

	bool Do_Phase1_Only = optionInfo->Do_phase1_only_bool;
	bool outputAlignInfoAndSamForAllPairedAlignmentBool = false;
	//outputAlignInfoAndSamForAllPairedAlignmentBool = true;
	bool removeAllIntermediateFilesBool = false;
	//removeAllIntermediateFilesBool = true;

	bool DoSam2JuncBool = false;
	DoSam2JuncBool = true;
	bool load2ndLevelIndexBool = false;
	load2ndLevelIndexBool = true;
	bool load2ndLevelIndexBool_compressedSize = false;
	load2ndLevelIndexBool_compressedSize = true;
	bool DoRemappingOnUnmapEndReadsBool = false;
	DoRemappingOnUnmapEndReadsBool = true;
	bool DoRemappingOnUnfixedHeadTailAlignmentBool = false;
	DoRemappingOnUnfixedHeadTailAlignmentBool = true;
	bool outputDirectlyBool_Phase1Only = true;
	outputDirectlyBool_Phase1Only = false;
	bool Do_cirRNA = true;
	Do_cirRNA = false;

	//bool Do_extendHeadTail = true;
	//Do_extendHeadTail = false;

	bool Do_extendHeadTail_phase1 = true;
	//Do_extendHeadTail_phase1 = false;
	bool Do_extendHeadTail_fixOneEndUnmapped = true;
	//Do_extendHeadTail_fixOneEndUnmapped = false;
	bool Do_extendHeadTail_fixHeadTail = true;
	//Do_extendHeadTail_fixHeadTail = false;	

	bool Do_fixHeadTail_remapping = true;
	//Do_fixHeadTail_remapping = false;
	bool Do_fixHeadTail_greedyMapping = true;
	//Do_fixHeadTail_greedyMapping = false;
	bool Do_fixHeadTail_remappingAndTargetMapping = true;
	//Do_fixHeadTail_remappingAndTargetMapping = false;
	bool Do_fixHeadTail_remappingAgain = true;
	Do_fixHeadTail_remappingAgain = false;
	bool Do_fixHeadTail_extend2end_finalStep = true;
	//Do_fixHeadTail_extend2end_finalStep = false;
	bool Do_fixHeadTail_extend2end_fixIndel = true;
	//Do_fixHeadTail_extend2end_fixIndel = false;

	// Fix ME:
	// #ifdef DETECT_CIRCULAR_RNA
	// Do_fixHeadTail_remapping = false;
	// Do_fixHeadTail_greedyMapping = false;
	// Do_fixHeadTail_remappingAndTargetMapping = false;
	// #endif
	if(Do_Phase1_Only)
	{
		settings_log_ofs << "Do_Phase1 only!" << endl;
		DoSam2JuncBool = false;
		load2ndLevelIndexBool = false;
		load2ndLevelIndexBool_compressedSize = false;
		DoRemappingOnUnmapEndReadsBool = false;
		DoRemappingOnUnfixedHeadTailAlignmentBool = false;
	}	
	else
	{
		settings_log_ofs << "Do_Phase1_Phase2! " << endl;
		DoSam2JuncBool = true;//false;
		load2ndLevelIndexBool = true;//false;
		load2ndLevelIndexBool_compressedSize = true;//false;
		DoRemappingOnUnmapEndReadsBool = true;//false;
		DoRemappingOnUnfixedHeadTailAlignmentBool = true;//false;
	}

	int normalRecordNum_1stMapping = 500000;
	int normalRecordNum_fixOneEndUnmapped = 500000;
	int normalRecordNum_fixHeadTail = 500000;

	int readTotalNum = 0;

	optionInfo->outputSwitchInfo(Do_Phase1_Only, outputAlignInfoAndSamForAllPairedAlignmentBool,
		removeAllIntermediateFilesBool, Do_cirRNA, outputDirectlyBool_Phase1Only, 
		normalRecordNum_1stMapping, normalRecordNum_fixOneEndUnmapped,
		normalRecordNum_fixHeadTail, Do_extendHeadTail_phase1, 
		Do_extendHeadTail_fixOneEndUnmapped, Do_extendHeadTail_fixHeadTail, 
		Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
		Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, settings_log_ofs);
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... MPS starts ......" << endl << endl;  
	log_ofs << endl << "[" << asctime(local) << "... MPS starts ......" << endl << endl; 

	//////////////////////////////////////////////        LOAD INDEX         ////////////////////////////////////////////////////////////////

    string InputReadFile = optionInfo->read_file_path_1;//read sample, exacted from fastq file every time
    if(SE_or_PE_bool)
    	InputReadFile = optionInfo->read_file_path_SE;
    string InputReadFile_PE = optionInfo->read_file_path_2;// another end read for pair-end reads

	int threads_num = optionInfo->threads_num;//atoi(threadsNumStr.c_str());

	#ifdef CAL_TIME	
	threads_num = 1;
	#endif

	#ifdef MAP_INFO
	threads_num = 1;
	#endif

	omp_set_num_threads(threads_num);

	bool InputAsFastq = (!(optionInfo->fasta_or_fastq_bool));
	bool fasta_or_fastq_bool = optionInfo->fasta_or_fastq_bool;
	if(checkQualSeqForShortAnchorSeqToTargetMap && fasta_or_fastq_bool)
	{
		cout << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		log_ofs << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		exit(1);
	}
	/////////////////////////////////////          LOAD INDEX         ////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 

    //cout << "start to load preIndex ..." << endl;
    log_ofs << "start to load preIndex ..." << endl;
	   

    string preIndexArrayPreStr;
    string indexStr;// = "/data/homes/lxauky/adSA_table/mm9_table/testAll2_table/testAll2Index";
    string chromDirStr;
    string secondLevelIndexStr;

    indexStr = optionInfo->global_index_file_path_prefix; //argv[6];
    preIndexArrayPreStr = indexStr;
    chromDirStr = optionInfo->chromsome_file_path_prefix; //argv[8];
    secondLevelIndexStr = optionInfo->local_index_file_path_prefix; //argv[7];

    preIndexArrayPreStr.append("/");
    indexStr.append("/");
    chromDirStr.append("/");
    secondLevelIndexStr.append("/");

	string preIndexMapLengthArrayStr = preIndexArrayPreStr; preIndexMapLengthArrayStr.append("_MapLength"); ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr = preIndexArrayPreStr; preIndexIntervalStartArrayStr.append("_IntervalStart"); ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr = preIndexArrayPreStr; preIndexIntervalEndArrayStr.append("_IntervalEnd"); ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);
	int* preIndexMapLengthArray; preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int)); preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray; preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray; preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));
 	log_ofs << "finish loading preIndex ..." << endl;

	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, settings_log_ofs);

	settings_log_ofs << "index: " << indexStr << endl;
	/////////////////////////////////////// 
	log_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	settings_log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	settings_log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
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

	//int *lcpCompress_test; lcpCompress_test = (int*)malloc((indexInfo->returnIndexSize()) * sizeof(int)); 	
	//int *verifyChild_test; verifyChild_test = (int*)malloc((indexInfo->returnIndexSize()) * sizeof(int)); 	

	ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl;

	log_ofs << "start to load SA" << endl;
	SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load lcpCompress" << endl;
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_ofs << "start to load childTab " << endl;
	childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	log_ofs << "start to load detChild" << endl;
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_ofs << "All index files loaded" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl;
	runtime_log_ofs << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl;	

	//////////////////////////////////////////////////
	HeaderSection_Info* headerSectionInfo = new HeaderSection_Info(indexInfo);	
	/////////////////////////////////////   start to load annotation  /////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load annotation file (SJs)......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load annotation file (SJs) ......" << endl << endl; 	
	string annotation_file_path = optionInfo->annotation_file_path; // junction files
	ifstream annotation_ifs(annotation_file_path.c_str());
	Annotation_Info* annotationInfo = new Annotation_Info();
	if(annotation_provided_bool)
		annotationInfo->initiateAndReadAnnotationFile(indexInfo, annotation_ifs);
	/////////////////////////////////////   finish loading annotation  /////////////////////////////////////	
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	bool SJalignInferHash_provided_bool = optionInfo->spliceJunctionAlignInferHash_provided_bool;
	cout << endl << "start to initaite alignInferJunctionHashInfo " << endl;
	log_ofs << endl << "start to initaite alignInferJunctionHashInfo " << endl;
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	if(SJalignInferHash_provided_bool)
	{	
		cout << "start to insert SJ into SJmap" << endl;
		string tmpInputJuncFile = optionInfo->spliceJunctionAlignInferHash_file_path;
		alignInferJunctionHashInfo->insertJuncFromJuncFile_onlyChrNamePos(tmpInputJuncFile, indexInfo);
		cout << "start to output SJ map" << endl;
	}
	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	#ifdef PERSONALIZED_CHR_SEQ
	bool SNP_provided_bool = optionInfo->return_SNP_provided_bool();
	if(!SNP_provided_bool)
	{
		cout << "Under PERSONALIZED mode, but SNP file is not provided with -P option, exiting ......" << endl;
		log_ofs << "Under PERSONALIZED mode, but SNP file is not provided with -P option, exiting ......" << endl;
		exit(1);
	}
	string SNPfilePath = optionInfo->SNPfilePath;
	indexInfo->insertSNP2chromStr(SNPfilePath, log_ofs);

	bool SNP_seq_index_provided_bool = optionInfo->return_SNP_seq_index_provided_bool();
	string SNP_seq_index_path;
	if(SNP_seq_index_provided_bool)
	{
		SNP_seq_index_path = optionInfo->SNP_seq_index_path;
	}
	cout << "start to load indexes" << endl;
	string indexStr_SNP = SNP_seq_index_path;
	indexStr_SNP.append("/");
	string SA_file_SNP = indexStr_SNP; SA_file_SNP.append("_SA"); ifstream SA_file_ifs_SNP(SA_file_SNP.c_str(),ios::binary); 
	string lcpCompress_file_SNP = indexStr_SNP; lcpCompress_file_SNP.append("_lcpCompress"); ifstream lcpCompress_file_ifs_SNP(lcpCompress_file_SNP.c_str(),ios::binary);
	string childTab_file_SNP = indexStr_SNP; childTab_file_SNP.append("_childTab"); ifstream childTab_file_ifs_SNP(childTab_file_SNP.c_str(),ios::binary);
	string verifyChild_file_SNP = indexStr_SNP; verifyChild_file_SNP.append("_detChild"); ifstream verifyChild_file_ifs_SNP(verifyChild_file_SNP.c_str(),ios::binary);	
	string chrom_bit_file_SNP = indexStr_SNP; chrom_bit_file_SNP.append("_chrom"); ifstream chrom_bit_file_ifs_SNP(chrom_bit_file_SNP.c_str(),ios::binary);
	string parameter_file_SNP = indexStr_SNP; parameter_file_SNP.append("_parameter"); ifstream parameter_file_ifs_SNP(parameter_file_SNP.c_str(),ios::binary);
	Index_Info* indexInfo_SNP = new Index_Info(parameter_file_ifs_SNP, log_ofs);
	char *chrom_SNP; chrom_SNP = (char*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs_SNP.read((char*)chrom_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(char)); 
	indexInfo_SNP->readGenome(chrom_SNP);
	indexInfo_SNP->initiate();	
	indexInfo_SNP->initiateChrNameIndexArray(1000);
    unsigned int *sa_SNP; sa_SNP = (unsigned int*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); SA_file_ifs_SNP.read((char*)sa_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));
	BYTE *lcpCompress_SNP; lcpCompress_SNP = (BYTE*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); lcpCompress_file_ifs_SNP.read((char*)lcpCompress_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));
	unsigned int *childTab_SNP; childTab_SNP = (unsigned int*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); childTab_file_ifs_SNP.read((char*)childTab_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));
	BYTE *verifyChild_SNP; verifyChild_SNP = (BYTE*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); verifyChild_file_ifs_SNP.read((char*)verifyChild_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));
	cout << "SyntheticSNPtransSeq index files loaded" << endl;
	#endif
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	Stats_Info* statsInfo = new Stats_Info();
	if(SE_or_PE_bool)
		statsInfo->initiate_stats_info_SE(threads_num);
	else
		statsInfo->initiate_stats_info_PE(threads_num);

	//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////
	/* align main*/
   	string tmpHeadSectionInfo = outputDirStr + "/headSectionInfo";
   	ofstream tmpHeadSectionInfo_ofs(tmpHeadSectionInfo.c_str());    	
   	tmpHeadSectionInfo_ofs << headerSectionInfo->returnHeaderSectionInfoStr() << endl;
	tmpHeadSectionInfo_ofs.close();

	string repeatRegionFile = outputDirStr + "/phase1_output/repeat_region/repeatRegion";
	ofstream repeatRegionFile_ofs(repeatRegionFile.c_str());
	string tmpAlignCompleteRead_SE = outputDirStr + "/phase1_output/completePair/complete_SE.sam";
	ofstream tmpAlignCompleteRead_SE_ofs(tmpAlignCompleteRead_SE.c_str());
	string tmpAlignCompleteRead = outputDirStr + "/phase1_output/completePair/completePair.sam";
	ofstream tmpAlignCompleteRead_ofs(tmpAlignCompleteRead.c_str());
	string tmpAlignOneEndUnmapped = outputDirStr + "/phase1_output/oneEndUnmapped/oneEndUnmapped";
	if(Do_Phase1_Only)
		tmpAlignOneEndUnmapped += ".sam";	
	else
		tmpAlignOneEndUnmapped += ".alignInfo";
	ofstream tmpAlignOneEndUnmapped_ofs(tmpAlignOneEndUnmapped.c_str());


	//string mkdirOutputCommand_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile = mkdirOutputCommand_phase1 + "/bothEndsUnmapped_mappedToRepeatRegion";
	//system(mkdirOutputCommand_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());
	string tmpAlignUnmapped_mappedToRepeatRegionFile_SE = outputDirStr + "/phase1_output/repeat_region/unmapped_mappedToRepeatRegion_SE.sam";
	ofstream tmpAlignUnmapped_mappedToRepeatRegionFile_SE_ofs(tmpAlignUnmapped_mappedToRepeatRegionFile_SE.c_str());
	string tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile 
		= outputDirStr + "/phase1_output/repeat_region/bothEndsUnmapped_mappedToRepeatRegion.sam";
	ofstream tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs(
		tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());

	string mkdirOutputCommand_tmpAlignBothEndsUnmapped = mkdirOutputCommand_phase1 + "/bothEndsUnmapped";
	system(mkdirOutputCommand_tmpAlignBothEndsUnmapped.c_str());
	string tmpAlignUnmapped_SE = outputDirStr + "/phase1_output/bothEndsUnmapped/unmapped_SE.sam";
	ofstream tmpAlignUnmapped_SE_ofs(tmpAlignUnmapped_SE.c_str());
	string tmpAlignBothEndsUnmapped = outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped.sam";
	ofstream tmpAlignBothEndsUnmapped_ofs(tmpAlignBothEndsUnmapped.c_str());
	string tmpAlignUnmapped_lowScore_SE = outputDirStr + "/phase1_output/bothEndsUnmapped/unmapped_lowScore_SE.sam";
	ofstream tmpAlignUnmapped_lowScore_SE_ofs(tmpAlignUnmapped_lowScore_SE.c_str());
	string tmpAlignBothEndsUnmapped_lowScore = outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped_lowScore.sam";
	ofstream tmpAlignBothEndsUnmapped_lowScore_ofs(tmpAlignBothEndsUnmapped_lowScore.c_str());

	string mkdirOutputCommand_tmpAlignIncompletePair = mkdirOutputCommand_phase1 + "/incomplete";
	system(mkdirOutputCommand_tmpAlignIncompletePair.c_str());
	string tmpAlignIncomplete_SE = outputDirStr + "/phase1_output/incomplete/incomplete_SE.alignInfo"; 
	ofstream tmpAlignIncomplete_SE_ofs(tmpAlignIncomplete_SE);
	string tmpAlignIncompletePair = outputDirStr + "/phase1_output/incomplete/incomplete.alignInfo"; 
	ofstream tmpAlignIncompletePair_ofs(tmpAlignIncompletePair.c_str());
	string tmpAlignIncomplete_SE_SAM = outputDirStr + "/phase1_output/incomplete/incomplete_SE.sam";
	ofstream tmpAlignIncomplete_SE_SAM_ofs(tmpAlignIncomplete_SE_SAM.c_str());
	string tmpAlignIncompletePair_SAM = outputDirStr + "/phase1_output/incomplete/incompletePair.sam"; 
	ofstream tmpAlignIncompletePair_SAM_ofs(tmpAlignIncompletePair_SAM.c_str());	

	string tmpIntermediateJunctionFile = outputDirStr + "/phase2_output/inter.junc";

	ifstream inputRead_ifs(InputReadFile.c_str());
	ifstream inputRead_PE_ifs(InputReadFile_PE.c_str());
	//ifstream inputRead_SE_ifs(InputReadFile_SE.c_str());

	if(SE_or_PE_bool)
	{
		//inputRead_ifs.close();
		inputRead_PE_ifs.close();
	}
	else
	{
		//inputRead_SE_ifs.close();
	}

    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE;
    string line2_afterProcess, line2_PE_afterProcess;

	int normalRecordNum = normalRecordNum_1stMapping; //1000000;//1500000;

	#ifdef DEBUG_INFO
	normalRecordNum = 1;
	#endif

	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;

	int readPairNum = 0;

	// *************** needed for both SE and PE *************//
	vector<string> readName1Vec(normalRecordNum);
	vector<string> readSeq1Vec(normalRecordNum);
	vector<string> readQualSeq1Vec(normalRecordNum);
	// *******************************************************//
	// ***************  only needed for PE  ******************//
	vector<string> readName2Vec(normalRecordNum);
	vector<string> readSeq2Vec(normalRecordNum);
	vector<string> readQualSeq2Vec(normalRecordNum);

	vector<string> PeAlignSamStrVec_complete(normalRecordNum);
	vector<string> PeAlignInfoStrVec_inCompletePair(normalRecordNum);
	vector<string> PeAlignInfoStrVec_oneEndUnmapped(normalRecordNum);
	vector<string> PeAlignSamStrVec_bothEndsUnmapped(normalRecordNum);
	vector<string> PeAlignSamStrVec_bothEndsUnmapped_lowScore(normalRecordNum);
	vector<string> PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion(normalRecordNum);
	vector<string> PeAlignSamStrVec_inCompletePair(normalRecordNum);
	vector<string> PeAlignInfoStrVec_completePaired(normalRecordNum);
	// ******************************************************//
	// ***************  only needed for SE  ******************//
	vector<string> SeAlignSamStrVec_complete(normalRecordNum);
	vector<string> SeAlignInfoStrVec_inComplete(normalRecordNum);
	vector<string> SeAlignSamStrVec_inComplete(normalRecordNum);
	vector<string> SeAlignSamStrVec_unmapped(normalRecordNum);
	vector<string> SeAlignSamStrVec_unmapped_lowScore(normalRecordNum);
	vector<string> SeAlignSamStrVec_unmapped_mappedToRepeatRegion(normalRecordNum);
	// *******************************************************//
	vector< RepeatRegion_Info* > repeatRegionInfoVec;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		RepeatRegion_Info* repeatRegionInfo = new RepeatRegion_Info();
		repeatRegionInfoVec.push_back(repeatRegionInfo);
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... 1st mapping process starts ......" << endl << endl; 

	InputReadPreProcess* readPreProcessInfo = new InputReadPreProcess();

	for(tmpTurn = 0; 
		//tmpTurn <= 300     //used to control # of rounds to process
		; tmpTurn++)
	{
		if(EndOfRecord)
			break;
		int recordNum = normalRecordNum;
		//cout << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		//cout << endl << "[" << asctime(local) 
		//	<< "start to read Fasta file, turn: " << tmpTurn + 1 << endl;			
		realRecordNum = normalRecordNum;

		#ifdef CAL_TIME
		input_begin = clock();
		#endif

		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}

    		getline(inputRead_ifs, line1); // readName_1
    		//cout << "readName_1: " << endl << line1 << endl;
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		readName1Vec[recordNumTmp] = line1.substr(1);
    		getline(inputRead_ifs, line2); // readSeq_1
    		//cout << "readSeq_ori_1: " << endl << line2 << endl;
    		line2_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2);		
    		readSeq1Vec[recordNumTmp] = line2_afterProcess;
    		//cout << "afterProcessing: " << endl << line2_afterProcess << endl;
    		if(InputAsFastq)
    		{
    			getline(inputRead_ifs, line3);
    			getline(inputRead_ifs, line4);
    			readQualSeq1Vec[recordNumTmp] = line4;   
    		}
    			
    		if(!SE_or_PE_bool)
    		{	
	    		getline(inputRead_PE_ifs, line1_PE); // readName_2
	    		readName2Vec[recordNumTmp] = line1_PE.substr(1);
	    		getline(inputRead_PE_ifs, line2_PE); // readSeq_2
	    		line2_PE_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2_PE);
	    		readSeq2Vec[recordNumTmp] = line2_PE_afterProcess;
	    		if(InputAsFastq)
	    		{
	    			getline(inputRead_PE_ifs, line3_PE);
	    			getline(inputRead_PE_ifs, line4_PE);
	    			readQualSeq2Vec[recordNumTmp] = line4_PE;
	    		}
    		}
		}

		#ifdef CAL_TIME
		input_end = clock();
		input_cost += (input_end - input_begin);
		process_begin = clock();
		#endif

		readTotalNum += realRecordNum;
		runtime_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);		
		runtime_log_ofs << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		runtime_log_ofs << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;
		//cout << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		//cout << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;
		omp_set_num_threads(threads_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			int threadNO = omp_get_thread_num();
			#ifdef CAL_TIME
			getReadInfo_begin = clock();
			#endif
			PE_Read_Info readInfo; //= new PE_Read_Info();
			readInfo.initiateReadInfo(readName1Vec[tmpOpenMP], readName2Vec[tmpOpenMP],
				readSeq1Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP],
				readQualSeq1Vec[tmpOpenMP], readQualSeq2Vec[tmpOpenMP], fasta_or_fastq_bool, SE_or_PE_bool);
			// //bool peReadsOverlapped_bool = readInfo.checkPEreadsOverlapped();
			// cout << "###############################" << endl;	
			// cout << endl << "read_name_1: " << readName1Vec[tmpOpenMP] << endl;
			// cout << "read_name_2: " << readName2Vec[tmpOpenMP] << endl;
			// cout << "readSeq_1: " << readSeq1Vec[tmpOpenMP] << endl;
			// cout << "readSeq_2: " << readSeq2Vec[tmpOpenMP] << endl;	
			#ifdef CAL_TIME
			getReadInfo_end = clock();
			getReadInfo_cost += (getReadInfo_end - getReadInfo_begin);
			getSegInfo_begin = clock();
			#endif
    		FixPhase1Info fixPhase1Info;
    		//cout << "start to do segInfo" << endl;
			fixPhase1Info.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
				preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray,
				readInfo, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, SE_or_PE_bool);
			// cout << "\nsegInfo_Nor1: " << endl; fixPhase1Info.coutSegInfo_Nor1(indexInfo);
			// cout << "segInfo_Rcm1: " << endl; fixPhase1Info.coutSegInfo_Rcm1(indexInfo);
			// cout << "segInfo_Nor2: " << endl; fixPhase1Info.coutSegInfo_Nor2(indexInfo);
			// cout << "segInfo_Rcm2: " << endl; fixPhase1Info.coutSegInfo_Rcm2(indexInfo);
			#ifdef CAL_TIME
			getSegInfo_end = clock();
			getSegInfo_cost += (getSegInfo_end - getSegInfo_begin);
			updateSegInfo_begin = clock();
			#endif
			
			#ifdef CAL_TIME
			updateSegInfo_greedyMap_begin = clock();
			#endif

			#ifdef PERSONALIZED_CHR_SEQ
			# ifdef VARY_SNP_MER
			fixPhase1Info.fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo_varySNPmer(
				sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, 
				readInfo, indexInfo);
			# else
			fixPhase1Info.fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo(
				sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, 
				readInfo, indexInfo, SNPlocInSyntheticSNPseq);
			# endif

			#ifdef CAL_TIME
			updateSegInfo_greedyMap_end = clock();
			updateSegInfo_greedyMap_cost += (updateSegInfo_greedyMap_end - updateSegInfo_greedyMap_begin);
			updateSegInfo_targetMap_begin = clock();
			#endif
			
			// cout << "###############################" << endl << "after updateing with snpSeq greedyMapping map results" << endl;
			// cout << "\nsegInfo_Nor1: " << endl; fixPhase1Info.coutSegInfo_Nor1(indexInfo);
			// cout << "segInfo_Rcm1: " << endl; fixPhase1Info.coutSegInfo_Rcm1(indexInfo);
			// cout << "segInfo_Nor2: " << endl; fixPhase1Info.coutSegInfo_Nor2(indexInfo);
			// cout << "segInfo_Rcm2: " << endl; fixPhase1Info.coutSegInfo_Rcm2(indexInfo);			

			# ifdef VARY_SNP_MER
			fixPhase1Info.fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo_varySNPmer(
				sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, 
				readInfo, indexInfo);
			# else
			fixPhase1Info.fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo(
				sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, 
				readInfo, indexInfo, SNPlocInSyntheticSNPseq);
			# endif			

			#ifdef CAL_TIME
			updateSegInfo_targetMap_end = clock();
			updateSegInfo_targetMap_cost += (updateSegInfo_targetMap_end - updateSegInfo_targetMap_begin);
			#endif

			// cout << "###############################" << endl << "after updateing with snpSeq targetMapping map results" << endl;
			// cout << "\nsegInfo_Nor1: " << endl; fixPhase1Info.coutSegInfo_Nor1(indexInfo);
			// cout << "segInfo_Rcm1: " << endl; fixPhase1Info.coutSegInfo_Rcm1(indexInfo);
			// cout << "segInfo_Nor2: " << endl; fixPhase1Info.coutSegInfo_Nor2(indexInfo);
			// cout << "segInfo_Rcm2: " << endl; fixPhase1Info.coutSegInfo_Rcm2(indexInfo);
			#endif

			#ifdef CAL_TIME
			updateSegInfo_end = clock();
			updateSegInfo_cost += (updateSegInfo_end - updateSegInfo_begin);
			getPathInfo_begin = clock();
			#endif

			fixPhase1Info.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, SE_or_PE_bool);

			#ifdef MAP_INFO
			cout << "getPossiPathInfo_Nor1: " << endl; fixPhase1Info.coutPossiPathStr_Nor1();
			cout << "getPossiPathInfo_Rcm1: " << endl; fixPhase1Info.coutPossiPathStr_Rcm1();
			#endif

			#ifdef CAL_TIME
			getPathInfo_end = clock();
			getPathInfo_cost += (getPathInfo_end - getPathInfo_begin);
			getGapInfo_begin = clock();
			#endif

			fixPhase1Info.fixPhase1_gapInfo(readInfo, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, SE_or_PE_bool);

			#ifdef MAP_INFO
			cout << "output map info: " << endl;
			fixPhase1Info.coutDebugInfo(readInfo, indexInfo, SE_or_PE_bool);
			cout << "finish output map info ..." << endl;
			#endif

			#ifdef CAL_TIME
			getGapInfo_end = clock();
			getGapInfo_cost += (getGapInfo_end - getGapInfo_begin);
			getAlignInfo_begin = clock();
			#endif

			PE_Read_Alignment_Info peAlignInfo;
			peAlignInfo.initiatePeAlignInfo( fixPhase1Info.pathInfo_Nor1, fixPhase1Info.pathInfo_Rcm1, 
				fixPhase1Info.pathInfo_Nor2, fixPhase1Info.pathInfo_Rcm2, indexInfo, SE_or_PE_bool);

			#ifdef CAL_TIME
			getAlignInfo_end = clock();
			getAlignInfo_cost += (getAlignInfo_end - getAlignInfo_begin);
			selectAlignInfo_begin = clock();
			#endif
			//cout << "###############################" << endl; cout << "peAlignInfo: \n" << peAlignInfo.returnPeAlignInfoStr() << endl;			
			//#ifdef MAP_INFO
			//cout << endl << endl << "peAlignInfo: \n" << peAlignInfo.returnPeAlignInfoStr() << endl << endl;
			//cout << "start to do alignments filtering and selection" << endl;
			//#endif
			if(Do_Phase1_Only)
			{
				if(SE_or_PE_bool)
					peAlignInfo.chooseBestAlignment_final_SE();
				else
					peAlignInfo.chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();
			}
			else
			{	
				if(SE_or_PE_bool)
					peAlignInfo.alignmentFilter_fixPhase1_SE(readInfo.returnReadLength_SE());
				else
					peAlignInfo.alignmentFilter_fixPhase1_SJpenalty(readInfo.returnReadSeqLength_1(),
						readInfo.returnReadSeqLength_2());
			}		
			//cout << "after filtering and choosing best alignments: " << endl; peAlignInfo.cout_finalAlignPair();	
			int overlapLength_top2candiAlignment = peAlignInfo.getMaxOverlapLengthInAlignmentPair(
				readInfo.returnReadSeqLength_1(), readInfo.returnReadSeqLength_2());

			#ifdef MAP_INFO
			cout << "alignments filtering ends ..." << endl;
			cout << "starts to record alignment ..." << endl; 
			#endif

			#ifdef CAL_TIME
			selectAlignInfo_end = clock();
			selectAlignInfo_cost += (selectAlignInfo_end - selectAlignInfo_begin);
			printAlignInfo_begin = clock();
			#endif			
			
			if(SE_or_PE_bool)
			{	
				// **************** SE ******************* //
				SeAlignSamStrVec_complete[tmpOpenMP] = "";
				SeAlignInfoStrVec_inComplete[tmpOpenMP] = "";
				SeAlignSamStrVec_inComplete[tmpOpenMP] = "";
				SeAlignSamStrVec_unmapped[tmpOpenMP] = "";
				SeAlignSamStrVec_unmapped_lowScore[tmpOpenMP] = "";
				SeAlignSamStrVec_unmapped_mappedToRepeatRegion[tmpOpenMP] = "";
			}
			else
			{
				// ***************  PE *******************//
				PeAlignSamStrVec_complete[tmpOpenMP] = "";			
				PeAlignInfoStrVec_inCompletePair[tmpOpenMP] = "";
				PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = "";
				PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = "";
				PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmpOpenMP] = "";
				PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmpOpenMP] = "";
				PeAlignSamStrVec_inCompletePair[tmpOpenMP] = "";
			}

			if(SE_or_PE_bool)
			{
				peAlignInfo.output_phase1_SE(Do_Phase1_Only, SeAlignSamStrVec_complete, SeAlignInfoStrVec_inComplete, 
					SeAlignSamStrVec_inComplete, SeAlignSamStrVec_unmapped, SeAlignSamStrVec_unmapped_lowScore, 
					SeAlignSamStrVec_unmapped_mappedToRepeatRegion,
					repeatRegionInfoVec, readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo, fasta_or_fastq_bool);
			}
			else
			{	
				peAlignInfo.output_phase1(outputDirectlyBool_Phase1Only, Do_Phase1_Only, PeAlignSamStrVec_complete,
					PeAlignInfoStrVec_inCompletePair, PeAlignInfoStrVec_oneEndUnmapped, PeAlignSamStrVec_bothEndsUnmapped,
					PeAlignSamStrVec_bothEndsUnmapped_lowScore, PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion,
					PeAlignSamStrVec_inCompletePair, repeatRegionInfoVec, readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo,
					fasta_or_fastq_bool, overlapLength_top2candiAlignment);
			}
			fixPhase1Info.memoryFree();
			peAlignInfo.memoryFree();

			#ifdef CAL_TIME
			printAlignInfo_end = clock();
			printAlignInfo_cost += (printAlignInfo_end - printAlignInfo_begin);
			#endif	

		} // read file end
		
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		runtime_log_ofs << "start to output ... turn: " << tmpTurn+1 << endl;
		//cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		//cout << "start to output ... turn: " << tmpTurn+1 << endl;

		#ifdef CAL_TIME
		process_end = clock();
		process_cost += (process_end - process_begin);
		output_begin = clock();
		#endif

		if(SE_or_PE_bool)
		{
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				if(SeAlignSamStrVec_complete[tmp] != "")
					tmpAlignCompleteRead_SE_ofs << SeAlignSamStrVec_complete[tmp] << endl;
				if(SeAlignInfoStrVec_inComplete[tmp] != "")
					tmpAlignIncomplete_SE_ofs << SeAlignInfoStrVec_inComplete[tmp] << endl;
				if(SeAlignSamStrVec_inComplete[tmp] != "")
					tmpAlignIncomplete_SE_SAM_ofs << SeAlignSamStrVec_inComplete[tmp] << endl;
				if(SeAlignSamStrVec_unmapped[tmp] != "")
					tmpAlignUnmapped_SE_ofs << SeAlignSamStrVec_unmapped[tmp] << endl;
				if(SeAlignSamStrVec_unmapped_lowScore[tmp] != "")
					tmpAlignUnmapped_lowScore_SE_ofs << SeAlignSamStrVec_unmapped_lowScore[tmp] << endl;
				if(SeAlignSamStrVec_unmapped_mappedToRepeatRegion[tmp] != "")
					tmpAlignUnmapped_mappedToRepeatRegionFile_SE_ofs << SeAlignSamStrVec_unmapped_mappedToRepeatRegion[tmp] << endl;
			}
		}
		else
		{
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{	
				if(PeAlignSamStrVec_complete[tmp] != "")
					tmpAlignCompleteRead_ofs << PeAlignSamStrVec_complete[tmp] << endl;		
				if(PeAlignInfoStrVec_inCompletePair[tmp] != "")
					tmpAlignIncompletePair_ofs << PeAlignInfoStrVec_inCompletePair[tmp] << endl;
				if(PeAlignInfoStrVec_oneEndUnmapped[tmp] != "")
					tmpAlignOneEndUnmapped_ofs << PeAlignInfoStrVec_oneEndUnmapped[tmp] << endl;
				if(PeAlignSamStrVec_bothEndsUnmapped[tmp] != "")
					tmpAlignBothEndsUnmapped_ofs << PeAlignSamStrVec_bothEndsUnmapped[tmp] << endl;
				if(PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmp] != "")
					tmpAlignBothEndsUnmapped_lowScore_ofs << PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmp] << endl;
				if(PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] != "")
					tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs << PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] << endl;
				if(PeAlignSamStrVec_inCompletePair[tmp] != "")
					tmpAlignIncompletePair_SAM_ofs << PeAlignSamStrVec_inCompletePair[tmp] << endl;
			}
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) << "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;

		#ifdef CAL_TIME
		output_end = clock();
		output_cost += (output_end - output_begin);
		#endif

	}
	#ifdef CAL_TIME
	cout << "input_cost: " << (double)input_cost/CLOCKS_PER_SEC << endl;
	cout << "process_cost: " << (double)process_cost/CLOCKS_PER_SEC << endl;
	cout << "output_cost: " << (double)output_cost/CLOCKS_PER_SEC << endl << endl;

	cout << "getReadInfo_cost: " << (double)getReadInfo_cost/CLOCKS_PER_SEC << endl;
	cout << "getSegInfo_cost: " << (double)getSegInfo_cost/CLOCKS_PER_SEC << endl;
	cout << "updateSegInfo_cost: " << (double)updateSegInfo_cost/CLOCKS_PER_SEC << endl;
	cout << "getPathInfo_cost: " << (double)getPathInfo_cost/CLOCKS_PER_SEC << endl;
	cout << "getGapInfo_cost: " << (double)getGapInfo_cost/CLOCKS_PER_SEC << endl;
	cout << "getAlignInfo_cost: " << (double)getAlignInfo_cost/CLOCKS_PER_SEC << endl;
	cout << "selectAlignInfo_cost: " << (double)selectAlignInfo_cost/CLOCKS_PER_SEC << endl;
	cout << "printAlignInfo_cost: " << (double)printAlignInfo_cost/CLOCKS_PER_SEC << endl << endl;

	cout << "updateSegInfo_greedyMap_cost: " << (double)updateSegInfo_greedyMap_cost/CLOCKS_PER_SEC << endl;
	cout << "updateSegInfo_targetMap_cost: " << (double)updateSegInfo_targetMap_cost/CLOCKS_PER_SEC << endl;
	#endif
	//log_ofs << "perfectMatch_pair #: " << perfectMatch_pair << endl;
	repeatRegionFile_ofs << "Repeat Region Info: size = " << repeatRegionInfoVec.size() << endl;
	for(int tmpThread = 0; tmpThread < threads_num; tmpThread++)
		repeatRegionInfoVec[tmpThread]->outputRepeatRegion(tmpThread+1, indexInfo, sa, 100, repeatRegionFile_ofs);

	settings_log_ofs << "readTotalNum: " << readTotalNum << endl;

	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	tmpAlignCompleteRead_ofs.close();
	tmpAlignOneEndUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_lowScore_ofs.close();
	tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs.close();
	tmpAlignIncompletePair_SAM_ofs.close();
	if(Do_Phase1_Only)
		tmpAlignIncompletePair_ofs.close();

	tmpAlignCompleteRead_SE_ofs.close();
	tmpAlignUnmapped_SE_ofs.close();
	tmpAlignUnmapped_lowScore_SE_ofs.close();
	tmpAlignUnmapped_mappedToRepeatRegionFile_SE_ofs.close();
	tmpAlignIncomplete_SE_SAM_ofs.close();
	if(Do_Phase1_Only)
		tmpAlignIncomplete_SE_ofs.close();

	free(preIndexMapLengthArray); free(preIndexIntervalStartArray); free(preIndexIntervalEndArray);
	free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
	free(childTab);
	free(verifyChild);
	free(chrom);
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
	log_ofs << endl << "**********************************" << endl << "**********************************";
	runtime_log_ofs << endl << "**********************************" << endl << "**********************************";

	if(SE_or_PE_bool)
	{
		statsInfo->getPhase1Stats_SE();
		statsInfo->outputAllStats_SE_phase1(stats_ofs, readTotalNum);
	}

	log_ofs << endl << "**********************************" << endl << "**********************************" << endl;
	//cout << endl << "totalReadNum = " << read_num << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////    	Load Second Level Index      ////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load 2nd level index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index starts ......" << endl << endl; 

	vector<char*> secondLevelChrom;
	vector<unsigned int*> secondLevelSa;

	vector<BYTE*> secondLevelLcpCompress;
	vector<unsigned int*> secondLevelChildTab;
	vector<BYTE*> secondLevelDetChild;

	if(load2ndLevelIndexBool)
	{
		log_ofs << "start to load second-level index ..." << endl;
		
		int secondLevelIndexNO = 0;
		for(int tmpChrNO = 0; tmpChrNO < indexInfo->returnChromNum(); tmpChrNO ++)
		{
			for(int tmpSecondLevelIndexNO = 1; tmpSecondLevelIndexNO <= (indexInfo->returnSecondLevelIndexPartsNum(tmpChrNO)); tmpSecondLevelIndexNO ++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpSecondLevelIndexNO);
				string tmpFileNumStr = tmpFileNumChar;
				
				string inputIndexFileStr = secondLevelIndexStr + "/" + 
					//indexInfo->chrNameStr[tmpChrNO] 
					indexInfo->returnChrNameStr(tmpChrNO)
					+ "/" 
					//+ indexInfo->chrNameStr[tmpChrNO] + "_part."
					+ tmpFileNumStr + "/";//"." + "test3_";


				string secondLevelIndexFileChromStr = inputIndexFileStr + "chrom"; 
				ifstream secondLevelChrom_file_ifs(secondLevelIndexFileChromStr.c_str(), ios::binary);
				string secondLevelIndexFileSaStr = inputIndexFileStr + "SA";
				ifstream secondLevelSA_file_ifs(secondLevelIndexFileSaStr.c_str(), ios::binary);

					string secondLevelIndexFileLcpCompressStr = inputIndexFileStr + "_lcpCompress";	
					ifstream secondLevelLcpCompress_file_ifs(secondLevelIndexFileLcpCompressStr.c_str(), ios::binary);	
					string secondLevelIndexFileChildTabStr = inputIndexFileStr + "childTab";	
					ifstream secondLevelChildTab_file_ifs(secondLevelIndexFileChildTabStr.c_str(), ios::binary);
					string secondLevelIndexFileDetChildStr = inputIndexFileStr + "detChild";	
					ifstream secondLevelDetChild_file_ifs(secondLevelIndexFileDetChildStr.c_str(), ios::binary);					

				int sizeOfIndex = indexInfo->returnSecondLevelIndexNormalSize() + 1;
				char* tmpSecondLevelChrom = (char*)malloc(sizeOfIndex * sizeof(char));
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					tmpSecondLevelChrom[tmpMallocSpace] = '0';
				}
				secondLevelChrom_file_ifs.read((char*)tmpSecondLevelChrom, sizeOfIndex * sizeof(char));
				if(tmpSecondLevelChrom[sizeOfIndex-1] != 'X')
				{
					//(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
				}

				bool No_ATGC_Bool = true;
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					char ch = tmpSecondLevelChrom[tmpMallocSpace];
					if((ch == 'A')||(ch == 'T')||(ch == 'G')||(ch == 'C'))
					{
						No_ATGC_Bool = false;
						break;
					}
				}				
				if(No_ATGC_Bool)
				{
					//(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
				}	

				secondLevelChrom.push_back(tmpSecondLevelChrom);
				
				unsigned int* tmpSecondLevelSa = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelSA_file_ifs.read((char*)tmpSecondLevelSa, sizeOfIndex * sizeof(unsigned int));
				secondLevelSa.push_back(tmpSecondLevelSa);

				BYTE* tmpSecondLevelLcpCompress = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelLcpCompress_file_ifs.read((char*)tmpSecondLevelLcpCompress, sizeOfIndex * sizeof(BYTE));
				secondLevelLcpCompress.push_back(tmpSecondLevelLcpCompress);
					
				unsigned int* tmpSecondLevelChildTab = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelChildTab_file_ifs.read((char*)tmpSecondLevelChildTab, sizeOfIndex * sizeof(unsigned int));
				secondLevelChildTab.push_back(tmpSecondLevelChildTab);

				BYTE* tmpSecondLevelDetChild = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelDetChild_file_ifs.read((char*)tmpSecondLevelDetChild, sizeOfIndex * sizeof(BYTE));
				secondLevelDetChild.push_back(tmpSecondLevelDetChild);

				secondLevelChrom_file_ifs.close();
				secondLevelSA_file_ifs.close();

				secondLevelLcpCompress_file_ifs.close();
				secondLevelChildTab_file_ifs.close();
				secondLevelDetChild_file_ifs.close();							

				secondLevelIndexNO ++;
			}
			log_ofs << "finish loading 2nd-level index of " << indexInfo->returnChrNameStr(tmpChrNO) << endl; 
		}
		log_ofs << "finish loading ALL 2nd-level index !" << endl;
		log_ofs << indexInfo->getInvalidSecondLevelIndexNOstr() << endl;
		//loadIndex_end = clock(); 
	}

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... load 2nd level index ends ......" << endl << endl ; 		
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index ends ......" << endl << endl ; 	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On one end unmapped Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << endl << "[" << asctime(local) 
		<< "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 	 
	log_ofs << endl << endl << "[" << asctime(local) 
		<< "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 
	runtime_log_ofs << endl << endl << "[" << asctime(local) 
		<< "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 

	string OutputSamFile_oneEndMapped = outputDirStr + "/phase2_output/oneEndUnmapped.pairedComplete.sam";
	ofstream OutputSamFile_oneEndMapped_ofs(OutputSamFile_oneEndMapped.c_str());

	string OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore = outputDirStr + "/phase2_output/oneEndUnmapped.bothEndsUnmapped_lowScore.sam";
	ofstream OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs(OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore.c_str());

	string OutputSamFile_oneEndMapped_unpair = outputDirStr + "/phase2_output/oneEndUnmapped.unpaired.sam";
	ofstream OutputSamFile_oneEndMapped_unpair_ofs(OutputSamFile_oneEndMapped_unpair.c_str());

	string OutputSamFile_oneEndMapped_alignInfo = outputDirStr + "/phase2_output/oneEndUnmapped.pairedComplete.sam_alignInfo";
	ofstream OutputSamFile_oneEndMapped_alignInfo_ofs(OutputSamFile_oneEndMapped_alignInfo.c_str());	

	int tmpRecordNum_oneEndUnmapped = 0;

	if(SE_or_PE_bool)
		DoRemappingOnUnmapEndReadsBool = false;

	if(DoRemappingOnUnmapEndReadsBool)
	{
		nowtime = time(NULL);
		local = localtime(&nowtime);
		log_ofs << endl << endl << "[" << asctime(local) << "start doing remapping on unmapped end reads" << endl;
		runtime_log_ofs << endl << endl << "[" << asctime(local) << "start doing remapping on unmapped end reads" << endl;
		cout << endl << endl << "[" << asctime(local) << "start doing remapping on unmapped end reads" << endl;

		string oneEndMappedFileStr = tmpAlignOneEndUnmapped;
		ifstream inputRecord_ifs(oneEndMappedFileStr.c_str());
		
		int normalRecordNum = normalRecordNum_fixOneEndUnmapped; //1000000;
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;// = normalRecordNum;
		//getline(inputRecord_ifs, line11);

		for(tmpTurn = 0; /*tmpTurn < TurnNum*/; tmpTurn++)
		{
			vector<string> line1StrVec(normalRecordNum);
			vector<string> line2StrVec(normalRecordNum);
			vector<string> line3StrVec(normalRecordNum);
			vector<string> line4StrVec(normalRecordNum);
			vector<string> line5StrVec(normalRecordNum);
			vector<string> line6StrVec(normalRecordNum);
			vector<string> line7StrVec(normalRecordNum);
			vector<string> line8StrVec(normalRecordNum);
			vector<string> line9StrVec(normalRecordNum);
			vector<string> line10StrVec(normalRecordNum);
			vector<string> peAlignInfoVec_fixUnpair(normalRecordNum);
			vector<string> peAlignSamVec_fixUnpair(normalRecordNum);

			// if outputUnpairedSAM_bool = false, store all unpaired alignments
			// else, store all pairs of alignments within which all mapped alignments are complete (no partial alignments)
			vector<string> peAlignSamVec_unpair_fixUnpair(normalRecordNum);
			vector<string> peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair(normalRecordNum);

			if(EndOfRecord)
				break;
			int recordNum = normalRecordNum;
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << endl << "[" << asctime(local) 
				<< "start to input oneEndUnmapped records, turn: " << tmpTurn+1 << endl;
			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				string line1, line2, line3, line4, line5, line6, line7, 
					line8, line9, line10, line11;
				if(inputRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}				
				getline(inputRecord_ifs, line11);
				if(inputRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}					
				getline(inputRecord_ifs, line1);
				getline(inputRecord_ifs, line2);
				getline(inputRecord_ifs, line3);
				getline(inputRecord_ifs, line4);
				getline(inputRecord_ifs, line5);
				getline(inputRecord_ifs, line6);
				getline(inputRecord_ifs, line7);
				getline(inputRecord_ifs, line8);
				getline(inputRecord_ifs, line9);
				getline(inputRecord_ifs, line10);
				//getline(inputRecord_ifs, line11);

				line1StrVec[recordNumTmp] = line1;
				line2StrVec[recordNumTmp] = line2;
				line3StrVec[recordNumTmp] = line3;
				line4StrVec[recordNumTmp] = line4;
				line5StrVec[recordNumTmp] = line5;
				line6StrVec[recordNumTmp] = line6;
				line7StrVec[recordNumTmp] = line7;
				line8StrVec[recordNumTmp] = line8;
				line9StrVec[recordNumTmp] = line9;
				line10StrVec[recordNumTmp] = line10;
			}	
			runtime_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local) 
				<< "finish reading record, turn: " << tmpTurn+1 << endl;
			runtime_log_ofs << endl << "[" << asctime(local) 
				<< "start to fix oneEndUnmapped, turn: " << tmpTurn+1 << endl;	
			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for schedule(dynamic)
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				int threadNO = omp_get_thread_num();
				int multiMapSeg_maxLength = 0;

				PE_Read_Info peReadInfo;// = new PE_Read_Info();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_toFixOneEndUnmapped_getline(
					line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], line3StrVec[tmpOpenMP],
					line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP], line6StrVec[tmpOpenMP],
					line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
					line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo, indexInfo, 
					fasta_or_fastq_bool, SE_or_PE_bool, multiMapSeg_maxLength);
				
				//cout << "readName: " << peReadInfo.returnReadName_1() << endl;
				FixOneEndUnmappedInfo* fixOneEndUnmappedInfo = new FixOneEndUnmappedInfo();
				fixOneEndUnmappedInfo->fixOneEndUnmapped(peReadInfo, peAlignInfo,
					secondLevelChrom,
					secondLevelSa,
					secondLevelLcpCompress,
					secondLevelChildTab,
					secondLevelDetChild,
					indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE2,
					checkQualSeqForReadSegSeq);

				#ifdef PERSONALIZED_CHR_SEQ
				# ifdef VARY_SNP_MER 
				fixOneEndUnmappedInfo->fixOneEndUnmapped_includeSNPseqMap_varySNPmer(peReadInfo, peAlignInfo,
					secondLevelChrom,
					secondLevelSa,
					secondLevelLcpCompress,
					secondLevelChildTab,
					secondLevelDetChild,
					indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE2,
					checkQualSeqForReadSegSeq,
					sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP,
					verifyChild_SNP, indexInfo_SNP);//, SNPlocInSyntheticSNPseq);
				# else	
				fixOneEndUnmappedInfo->fixOneEndUnmapped_includeSNPseqMap(peReadInfo, peAlignInfo,
					secondLevelChrom,
					secondLevelSa,
					secondLevelLcpCompress,
					secondLevelChildTab,
					secondLevelDetChild,
					indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE2,
					checkQualSeqForReadSegSeq,
					sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP,
					verifyChild_SNP, indexInfo_SNP, SNPlocInSyntheticSNPseq);						
				# endif
				#endif

				peAlignInfo->alignmentFilter_fixOneEndUnmapped_SJpenalty(peReadInfo.returnReadSeqLength_1(),
					peReadInfo.returnReadSeqLength_2());

				bool pairExistsBool = peAlignInfo->finalPairExistsBool();
				bool allAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();
				bool allUnpairedAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();
				bool unique_bool = peAlignInfo->checkUniqueOrMulti();

				//cout << "pairExistsBool: " << pairExistsBool << endl;
				//cout << "allAlignmentCompleteBool: " << allAlignmentCompleteBool << endl;
				//cout << "allUnpairedAlignmentCompleteBool: " << allUnpairedAlignmentCompleteBool << endl;
				//cout << "unique_bool: " << unique_bool << endl;

				statsInfo->increNum_fixUnpaired(threadNO, pairExistsBool, allAlignmentCompleteBool, 
					allUnpairedAlignmentCompleteBool, unique_bool);

				peAlignSamVec_fixUnpair[tmpOpenMP] = "";
				peAlignInfoVec_fixUnpair[tmpOpenMP] = "";
				peAlignSamVec_unpair_fixUnpair[tmpOpenMP] = "";
				peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmpOpenMP] = "";

				if(pairExistsBool && allAlignmentCompleteBool) // some pair exists, all completed, print out paired SAM info
				{
					bool completeAlign_lowScore_bool = peAlignInfo->alignPairScoreTooLow_bool(peReadInfo);
					if(completeAlign_lowScore_bool)
					{
						statsInfo->increLowScoreComplete_fixUnpair(threadNO, unique_bool);
						peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmpOpenMP] 
							= peAlignInfo->getSAMformatForBothEndsUnmapped(peReadInfo, fasta_or_fastq_bool);					
					}
					else	
						peAlignSamVec_fixUnpair[tmpOpenMP]
							= peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo,
								fasta_or_fastq_bool, multiMapSeg_maxLength);
				}
				else if(pairExistsBool && (!allAlignmentCompleteBool)) // pair exists, incomplete
					peAlignInfoVec_fixUnpair[tmpOpenMP] 
						= peAlignInfo->getTmpAlignInfoForFinalPair(
						 		peReadInfo.returnReadName_1(), peReadInfo.returnReadName_2(), 
								peReadInfo.returnReadSeq_1(), peReadInfo.returnReadSeq_2(),
								peReadInfo.returnReadQual_1(), peReadInfo.returnReadQual_2(),
								fasta_or_fastq_bool, multiMapSeg_maxLength);
				else
				{
					//cout << "outputUnpairedSAM_bool: " << outputUnpairedSAM_bool << endl;
					if(!outputUnpairedSAM_bool)
						peAlignSamVec_unpair_fixUnpair[tmpOpenMP]
							= peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(
								peReadInfo, fasta_or_fastq_bool);
					else
					{
						bool partialAlignmentExistsOrNot_bool = peAlignInfo->partialAlignmentExists();
						//cout << "partialAlignmentExistsOrNot_bool " << endl;
						if(partialAlignmentExistsOrNot_bool)
						{
							peAlignInfoVec_fixUnpair[tmpOpenMP]
								= peAlignInfo->getTmpAlignInfo(
									peReadInfo.returnReadName_1(), peReadInfo.returnReadName_2(),
									peReadInfo.returnReadSeq_1(), peReadInfo.returnReadSeq_2(),
									peReadInfo.returnReadQual_1(), peReadInfo.returnReadQual_2(),
									fasta_or_fastq_bool, multiMapSeg_maxLength);
							//cout << "tmpAlignInfo: " << peAlignInfoVec_fixUnpair[tmpOpenMP] << endl;
						}
						else
						{
							peAlignInfo->chooseBestAlignment_final_PEasSE();
							if(outputUnpairedSAM_bothEndsUniqueMappedOnly_bool)
								peAlignSamVec_unpair_fixUnpair[tmpOpenMP]
									= peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_outputUniqueUnpairedBothEndsMappedOnly(
										peReadInfo, fasta_or_fastq_bool, multiMapSeg_maxLength);
							else 
								peAlignSamVec_unpair_fixUnpair[tmpOpenMP]
									= peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_allCases(
										peReadInfo, fasta_or_fastq_bool, multiMapSeg_maxLength);
							//cout << "tmpAlignSAM: " << peAlignSamVec_unpair_fixUnpair[tmpOpenMP] << endl;
						}
					}				
				}
				delete fixOneEndUnmappedInfo;
				peAlignInfo->memoryFree();
				delete peAlignInfo;
			}

			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish fixing oneEndUnmapped, turn: " << tmpTurn+1 << endl;// << endl;
			runtime_log_ofs << "start to output ... turn: " << tmpTurn+1 << endl;
			// cout << endl << "[" << asctime(local)
			// 	<< "finish fixing oneEndUnmapped, turn: " << tmpTurn+1 << endl;// << endl;
			// cout << "start to output ... turn: " << tmpTurn+1 << endl;
						
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				if(peAlignSamVec_fixUnpair[tmp] != "")
					OutputSamFile_oneEndMapped_ofs << peAlignSamVec_fixUnpair[tmp] << endl;
				if(peAlignInfoVec_fixUnpair[tmp] != "")
					tmpAlignIncompletePair_ofs << peAlignInfoVec_fixUnpair[tmp] << endl;
				if(peAlignSamVec_unpair_fixUnpair[tmp] != "")
					OutputSamFile_oneEndMapped_unpair_ofs << peAlignSamVec_unpair_fixUnpair[tmp] << endl;
				if(peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmp] != "")
					OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs << peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmp] << endl;
			}		
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish output, turn: " << tmpTurn+1 << endl << endl;// << endl;
			//cout << endl << "[" << asctime(local)
			//	<< "finish output, turn: " << tmpTurn+1 << endl << endl;// << endl;
		}		
		inputRecord_ifs.close();
	}

	OutputSamFile_oneEndMapped_ofs.close();
	OutputSamFile_oneEndMapped_unpair_ofs.close();
	tmpAlignIncompletePair_ofs.close();
	OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs.close();

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl;  
	log_ofs << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl; 
	
	log_ofs << endl << "**********************************************************************************" << endl << endl;
	runtime_log_ofs << endl << "**********************************************************************************" << endl << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Sam 2 Junc   ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... generating SJhashInfo starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... generating SJhashInfo starts ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... generating SJhashInfo starts ......" << endl << endl; 

	string juncfile = tmpIntermediateJunctionFile;
	string juncfile_alignInferHash = juncfile + ".alignInferHash";
	if(DoSam2JuncBool)
	{
		// generate SJ from already mapped reads
		vector<string> tmpAlignmentFileVec;
		if(SE_or_PE_bool)
		{
			tmpAlignmentFileVec.push_back(tmpAlignCompleteRead_SE);
			tmpAlignmentFileVec.push_back(tmpAlignIncomplete_SE_SAM);
		}
		else
		{
			tmpAlignmentFileVec.push_back(tmpAlignCompleteRead);
			tmpAlignmentFileVec.push_back(tmpAlignIncompletePair_SAM);
			tmpAlignmentFileVec.push_back(OutputSamFile_oneEndMapped);
		}

		AlignInferJunctionHash_Info_Vec* alignInferJunctionHashInfoVec 
			= new AlignInferJunctionHash_Info_Vec();
		int alignInferJuncHashInfoVecSize = threads_num;
		alignInferJunctionHashInfoVec->initiateAlignInferJunctionHashInfoVec(
			alignInferJuncHashInfoVecSize, chromNum);		
		alignInferJunctionHashInfoVec->insertJuncFromAlignmentFileVec_chrNamePos_supportNum_parallel(
			tmpAlignmentFileVec, indexInfo, alignInferJuncHashInfoVecSize, log_ofs);
		alignInferJunctionHashInfoVec->mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum(
			alignInferJunctionHashInfo, indexInfo);
		alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_chrNamePos_supportNum(
			indexInfo, juncfile_alignInferHash);
		alignInferJunctionHashInfoVec->freeMemory();
		delete alignInferJunctionHashInfoVec;
	}

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... generating SJhashInfo ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... generating SJhashInfo ends ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... generating SJhashInfo ends ......" << endl << endl; 
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On unfixed head/tail Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ; 

	string OutputSamFile_fixHeadTail_complete_pair = outputDirStr + "/phase2_output/fixHeadTail_complete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_pair_ofs(OutputSamFile_fixHeadTail_complete_pair.c_str());

	string OutputSamFile_fixHeadTail_incomplete_pair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_ofs(OutputSamFile_fixHeadTail_incomplete_pair.c_str());

	string OutputSamFile_fixHeadTail_complete_unpair = outputDirStr + "/phase2_output/fixHeadTail_complete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_unpair_ofs(OutputSamFile_fixHeadTail_complete_unpair.c_str());

	string OutputSamFile_fixHeadTail_incomplete_unpair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_unpair_ofs(OutputSamFile_fixHeadTail_incomplete_unpair.c_str());	

	string OutputSamFile_fixHeadTail_pair_lowScore = outputDirStr + "/phase2_output/fixHeadTail_pair_lowScore.sam";
	ofstream OutputSamFile_fixHeadTail_pair_lowScore_ofs(OutputSamFile_fixHeadTail_pair_lowScore.c_str());	

	string OutputSamFile_fixHeadTail_complete_SE = outputDirStr + "/phase2_output/fixHeadTail_complete_SE.sam";
	ofstream OutputSamFile_fixHeadTail_complete_SE_ofs(OutputSamFile_fixHeadTail_complete_SE.c_str());

	string OutputSamFile_fixHeadTail_incomplete_SE = outputDirStr + "/phase2_output/fixHeadTail_incomplete_SE.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_SE_ofs(OutputSamFile_fixHeadTail_incomplete_SE.c_str());

	string OutputSamFile_fixHeadTail_lowScore_SE = outputDirStr + "/phase2_output/fixHeadTail_lowScore_SE.sam";
	ofstream OutputSamFile_fixHeadTail_lowScore_SE_ofs(OutputSamFile_fixHeadTail_lowScore_SE.c_str());

	// start to read splice junction
	int junctionNum = 0;
	int junctionNum_in_alignInferJuncHashInfo = 0;
	int junctionNum_in_annotation = 0;

	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	

	string InputSpliceJunction = juncfile_alignInferHash;
	settings_log_ofs << "InputSpliceJunction: " << InputSpliceJunction << endl; 
	log_ofs << "InputSpliceJunction: " << InputSpliceJunction << endl; 
	cout << "InputSpliceJunction: " << InputSpliceJunction << endl;
	///////////////////////////////////////////////////////////////////
	if(DoRemappingOnUnfixedHeadTailAlignmentBool)
	{
    	log_ofs << "start to build spliceJunction Hash" << endl;
    	cout << "start to build spliceJunction Hash" << endl;
    	bool spliceJunctionHashExists = true;

		string entryString;
		int tabLocation1, tabLocation2, tabLocation3, tabLocation4, tabLocation5;
		//char entry[500];
		int chrInt;
		int spliceStartPos;
		int spliceEndPos;
		string chrIntString;
		string spliceStartPosString;
		string spliceEndPosString;

		/////////////////////////////////////////////////////////////////////////////
		//////////////////////  string hash /////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////

		if(!Do_annotation_only_bool)
		{	
	    	log_ofs << "start to load SJs in alignments" << endl;
    		cout << "start to load SJs in alignments" << endl;
			alignInferJunctionHashInfo->convert2SJhashInfo(SJ, indexInfo);
			junctionNum_in_alignInferJuncHashInfo = alignInferJunctionHashInfo->returnAlignInferInfoVecSize();
		}
		if(annotation_provided_bool)
		{	
			// loading SJs in annotation file (if provided)
	    	log_ofs << "start to load SJs in annotation" << endl;
    		cout << "start to load SJs in annotation" << endl;
			//FILE *fp_annotatedSJ_file = fopen(annotation_file_path.c_str(), "r");
			//fgets(entry, sizeof(entry), fp_annotatedSJ_file);
			ifstream annotatedSJ_ifs(annotation_file_path.c_str());
			while(!annotatedSJ_ifs.eof())
			{
				// fgets(entry, sizeof(entry), fp_annotatedSJ_file);
				// if(feof(fp_annotatedSJ_file))
				// 	break;
				getline(annotatedSJ_ifs, entryString);
				if(entryString == "")
					break;
				junctionNum_in_annotation ++;
				//entryString = entry;
				tabLocation1 = entryString.find('\t', 0);
				tabLocation2 = entryString.find('\t', tabLocation1+1);
				tabLocation3 = entryString.find('\t', tabLocation2+1);
				chrIntString = entryString.substr(0, tabLocation1);
				spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
				if(tabLocation3 == string::npos)
					spliceEndPosString = entryString.substr(tabLocation2+1);
				else
					spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
				chrInt = indexInfo->convertStringToInt(chrIntString);
				if(chrInt >= 0)
				{	
					spliceStartPos = atoi(spliceStartPosString.c_str());
					spliceEndPos = atoi(spliceEndPosString.c_str());	
					SJ->insert2AreaAndStringHash(chrInt, spliceStartPos, spliceEndPos, indexInfo);
				}
			}
			annotatedSJ_ifs.close();
			//fclose(fp_annotatedSJ_file);
		}
		junctionNum = junctionNum_in_alignInferJuncHashInfo + junctionNum_in_annotation;
		if(junctionNum == 0)
			spliceJunctionHashExists = false;

		log_ofs << "After inserting SJs generated from alignments and annotation" << endl;
		log_ofs << "\tjunctionNum in alignments = " << junctionNum_in_alignInferJuncHashInfo << endl;
		log_ofs << "\tjunctionNum in annotation = " << junctionNum_in_annotation << endl;
		log_ofs << "finish building spliceJunction Hash" << endl;		
		log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		settings_log_ofs << "After inserting SJs generated from alignments and annotation" << endl;
		settings_log_ofs << "\tjunctionNum in alignments = " << junctionNum_in_alignInferJuncHashInfo << endl;
		settings_log_ofs << "\tjunctionNum in annotation = " << junctionNum_in_annotation << endl;
		settings_log_ofs << "finish building spliceJunction Hash" << endl;		
		settings_log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		cout << "After inserting SJs generated from alignments and annotation" << endl;
		cout << "\tjunctionNum in alignments = " << junctionNum_in_alignInferJuncHashInfo << endl;
		cout << "\tjunctionNum in annotation = " << junctionNum_in_annotation << endl;
		cout << "finish building spliceJunction Hash" << endl;		
		cout << "start doing remapping on unfixed head/tail alignments" << endl;
		
		string headTailSoftClippingFile;
		if(SE_or_PE_bool)
			headTailSoftClippingFile = tmpAlignIncomplete_SE;
		else
			headTailSoftClippingFile = tmpAlignIncompletePair;// + ".all";

		if(outputDirectlyBool_Phase1Only)
			headTailSoftClippingFile += ".all";

		ifstream inputUnfixedHeadTailRecord_ifs(headTailSoftClippingFile.c_str());

		int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;// = normalRecordNum;

		for(tmpTurn = 0; /*tmpTurn < TurnNum*/; tmpTurn++)
		{
			vector<string> line1StrVec(normalRecordNum);
			vector<string> line2StrVec(normalRecordNum);
			vector<string> line3StrVec(normalRecordNum);
			vector<string> line4StrVec(normalRecordNum);
			vector<string> line5StrVec(normalRecordNum);
			vector<string> line6StrVec(normalRecordNum);
			vector<string> line7StrVec(normalRecordNum);
			vector<string> line8StrVec(normalRecordNum);
			vector<string> line9StrVec(normalRecordNum);
			vector<string> line10StrVec(normalRecordNum);
			vector<string> peAlignSamVec_complete_pair(normalRecordNum);
			vector<string> peAlignSamVec_incomplete_pair(normalRecordNum);
			vector<string> peAlignSamVec_complete_unpair(normalRecordNum);
			vector<string> peAlignSamVec_incomplete_unpair(normalRecordNum);
			vector<string> peAlignSamVec_pair_lowScore(normalRecordNum);

			vector<string> seAlignSamVec_complete(normalRecordNum);
			vector<string> seAlignSamVec_incomplete(normalRecordNum);
			vector<string> seAlignSamVec_lowScore(normalRecordNum);

			if(EndOfRecord)
				break;

			int recordNum = normalRecordNum;

			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "start to read Head/Tail file record, turn: " << tmpTurn+1 << endl;
			// cout << endl << "[" << asctime(local)
			// 	<< "start to read Head/Tail file record, turn: " << tmpTurn+1 << endl;

			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				string line1, line2, line3, line4, line5, line6, line7, 
					line8, line9, line10, line11;
				if(inputUnfixedHeadTailRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}				
				getline(inputUnfixedHeadTailRecord_ifs, line11);
				if(inputUnfixedHeadTailRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}		
				getline(inputUnfixedHeadTailRecord_ifs, line1);
				getline(inputUnfixedHeadTailRecord_ifs, line2);
				getline(inputUnfixedHeadTailRecord_ifs, line3);
				getline(inputUnfixedHeadTailRecord_ifs, line4);
				getline(inputUnfixedHeadTailRecord_ifs, line5);
				getline(inputUnfixedHeadTailRecord_ifs, line6);
				getline(inputUnfixedHeadTailRecord_ifs, line7);
				getline(inputUnfixedHeadTailRecord_ifs, line8);
				getline(inputUnfixedHeadTailRecord_ifs, line9);
				getline(inputUnfixedHeadTailRecord_ifs, line10);
				//getline(inputUnfixedHeadTailRecord_ifs, line11);

				line1StrVec[recordNumTmp] = line1;
				line2StrVec[recordNumTmp] = line2;
				line3StrVec[recordNumTmp] = line3;
				line4StrVec[recordNumTmp] = line4;
				line5StrVec[recordNumTmp] = line5;
				line6StrVec[recordNumTmp] = line6;
				line7StrVec[recordNumTmp] = line7;
				line8StrVec[recordNumTmp] = line8;
				line9StrVec[recordNumTmp] = line9;
				line10StrVec[recordNumTmp] = line10;
			}	

			runtime_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl;
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;					
			// cout << endl << "[" << asctime(local)
			// 	<< "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl;
			// cout << endl << "[" << asctime(local)
			// 	<< "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;					
			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for schedule(dynamic)
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				int threadNO = omp_get_thread_num();

				int multiMapSeg_maxLength;
				PE_Read_Info peReadInfo;// = new PE_Read_Info();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_toFixIncompleteAlignment_getline(
					line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], line3StrVec[tmpOpenMP],
					line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP], line6StrVec[tmpOpenMP],
					line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
					line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo, 
					indexInfo, fasta_or_fastq_bool, SE_or_PE_bool, multiMapSeg_maxLength);

				//cout << "multiMapSeg_maxLength: " << multiMapSeg_maxLength << endl;
if(Debug_SNPmap_Bool)
{			
			 	//#ifdef MAP_INFO
				cout << endl << endl << "readName_1: " << peReadInfo.returnReadName_1() << endl;
				cout << "readName_2: " << peReadInfo.returnReadName_2() << endl;
				cout << "start fixHeadTail: " << endl;
				cout << "PeAlignInfo:" << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
				//#endif
}
				FixHeadTailInfo* fixHeadTailInfo = new FixHeadTailInfo();
				if(Do_fixHeadTail_remapping)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail, SE_or_PE_bool);

					// fixHeadTailInfo->fixHeadTail_shortAnchorRemappingOnly_withAlignInfer(
					// 	peReadInfo, peAlignInfo, SJ, spliceJunctionHashExists,
					// 	indexInfo, alignInferJunctionHashInfo);

				}
if(Debug_SNPmap_Bool)
{	
				//#ifdef MAP_INFO
				cout << "after remapping:" << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
				//#endif
}						
				if(Do_fixHeadTail_greedyMapping)
				{
					#ifdef PERSONALIZED_CHR_SEQ
					# ifdef VARY_SNP_MER
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_greedyMappingOnly_includeSNPseqMap_varySNPmer(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
						MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq, SE_or_PE_bool,
						sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP,
						verifyChild_SNP, indexInfo_SNP);//, SNPlocInSyntheticSNPseq);
					# else
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_greedyMappingOnly_includeSNPseqMap(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
						MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq, SE_or_PE_bool,
						sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP,
						verifyChild_SNP, indexInfo_SNP, SNPlocInSyntheticSNPseq);	
					# endif
					#else
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_greedyMappingOnly(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
						MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq, SE_or_PE_bool);
					#endif
				}
if(Debug_SNPmap_Bool)
{			
				//#ifdef MAP_INFO
				cout << "after greedyMapping:" << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
				//#endif
}				
				if(Do_fixHeadTail_remappingAndTargetMapping)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingAndTargetMapping(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail,
						checkQualSeqForShortAnchorSeqToTargetMap, SE_or_PE_bool);	
				}
if(Debug_SNPmap_Bool)
{	
				//#ifdef MAP_INFO
				cout << "after remapping And Target Mapping:" << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;				
				//#endif
}
				if(Do_fixHeadTail_remappingAgain)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail, SE_or_PE_bool);	
				}

				#ifdef MAP_INFO
				cout << "start to do fixHeadTail_extend2end_finalStepForAligner:" << endl;
				#endif

				if(Do_fixHeadTail_extend2end_finalStep)
				{				
					fixHeadTailInfo->fixHeadTail_extend2end_finalStepForAligner(peReadInfo, peAlignInfo, indexInfo, SE_or_PE_bool);
				}

				#ifdef MAP_INFO
				cout << "start to do fixHeadTail_extend2end_fixIndel:" << endl;
				#endif				

				if(Do_fixHeadTail_extend2end_fixIndel)
				{
					fixHeadTailInfo->fixHeadTail_extend2end_fixIndel(peReadInfo, peAlignInfo, indexInfo, SE_or_PE_bool);	
				}
if(Debug_SNPmap_Bool)
{	
				//#ifdef MAP_INFO
				cout << "after extending: " << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;								
				//#endif
}
				// remove duplicate mismatch
				peAlignInfo->removeDuplicateMismatch(SE_or_PE_bool);

				if(SE_or_PE_bool)
					peAlignInfo->chooseBestAlignment_final_SE();
				else
					peAlignInfo->chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();

				if(SE_or_PE_bool)
				{
					seAlignSamVec_complete[tmpOpenMP] = "";
					seAlignSamVec_incomplete[tmpOpenMP] = "";
					seAlignSamVec_lowScore[tmpOpenMP] = "";					
					bool allFinalAlignComplete_bool = peAlignInfo->allFinalAlignmentComplete_SE();
					bool unique_bool = peAlignInfo->checkUniqueOrMulti_SE();
					bool alignmentScore_tooLow = peAlignInfo->alignScoreTooLow_bool_SE(peReadInfo);
				    statsInfo->increNum_fixHeadTail_SE(threadNO, allFinalAlignComplete_bool, unique_bool, alignmentScore_tooLow);
				    if(alignmentScore_tooLow)
				    	seAlignSamVec_lowScore[tmpOpenMP] 
				    		= peAlignInfo->getSAMformatForUnmapped_SE(peReadInfo, fasta_or_fastq_bool);
				    else if(allFinalAlignComplete_bool)
				    	seAlignSamVec_complete[tmpOpenMP] 
				    		= peAlignInfo->getSAMformatForFinalAlignment_SE(peReadInfo, fasta_or_fastq_bool);
				    else
				    	seAlignSamVec_incomplete[tmpOpenMP] 
				    		= peAlignInfo->getSAMformatForFinalAlignment_SE(peReadInfo, fasta_or_fastq_bool);
				}				
				else
				{	
					bool pairExistsBool = peAlignInfo->finalPairExistsBool();									
					peAlignSamVec_complete_pair[tmpOpenMP] = "";
					peAlignSamVec_incomplete_pair[tmpOpenMP] = "";
					peAlignSamVec_complete_unpair[tmpOpenMP] = "";
					peAlignSamVec_incomplete_unpair[tmpOpenMP] = "";
					peAlignSamVec_pair_lowScore[tmpOpenMP] = "";
					//cout << "pairExistsBool: " << pairExistsBool << endl;
					if(pairExistsBool) // some pair exists, all completed, print out paired SAM info
					{
						bool allFinalPairAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();	
						bool unique_bool = peAlignInfo->checkUniqueOrMulti();
						statsInfo->increPairedNum_fixHeadTail(threadNO, allFinalPairAlignmentCompleteBool, unique_bool);
						bool align_lowScore_bool = peAlignInfo->alignPairScoreTooLow_bool(peReadInfo);
						if(align_lowScore_bool)
						{
							statsInfo->increLowScoreComplete_fixHeadTail(
								threadNO, allFinalPairAlignmentCompleteBool, unique_bool);
							peAlignSamVec_pair_lowScore[tmpOpenMP] 
								= peAlignInfo->getSAMformatForBothEndsUnmapped(peReadInfo, fasta_or_fastq_bool);
						}
						else if(allFinalPairAlignmentCompleteBool)
							peAlignSamVec_complete_pair[tmpOpenMP] 
								= peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo, fasta_or_fastq_bool, multiMapSeg_maxLength);
						else
							peAlignSamVec_incomplete_pair[tmpOpenMP] 
								= peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo, fasta_or_fastq_bool, multiMapSeg_maxLength);		
					}
					else //if((!pairExistsBool) && (allAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
					{
						if(!outputUnpairedSAM_bool)
						{
							bool allUnpairAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();			
							//cout << "allUnpairAlignmentCompleteBool: " << allUnpairAlignmentCompleteBool << endl; 
							statsInfo->increUnpairedNum_fixHeadTail(threadNO, allUnpairAlignmentCompleteBool);							
							if(allUnpairAlignmentCompleteBool)
								peAlignSamVec_complete_unpair[tmpOpenMP] 
									= peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(
										peReadInfo, fasta_or_fastq_bool);
							else
								peAlignSamVec_incomplete_unpair[tmpOpenMP] 
									= peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(
										peReadInfo, fasta_or_fastq_bool);
						}
						else
						{
							//cout << "start to chooseBestAlignment_final_PEasSE ...." << endl;
							peAlignInfo->chooseBestAlignment_final_PEasSE();
							bool allUnpairAlignmentCompleteBool_final = peAlignInfo->allFinalUnpairedAlignmentCompleted();			
							//cout << "allUnpairAlignmentCompleteBool: " << allUnpairAlignmentCompleteBool << endl; 
							statsInfo->increUnpairedNum_fixHeadTail(threadNO, allUnpairAlignmentCompleteBool_final);
							if(outputUnpairedSAM_bothEndsUniqueMappedOnly_bool)
							{	
								if(allUnpairAlignmentCompleteBool_final)
									peAlignSamVec_complete_unpair[tmpOpenMP] 
										= peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_outputUniqueUnpairedBothEndsMappedOnly(
											peReadInfo, fasta_or_fastq_bool, multiMapSeg_maxLength);
								else
									peAlignSamVec_incomplete_unpair[tmpOpenMP] 
										= peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_outputUniqueUnpairedBothEndsMappedOnly(
											peReadInfo, fasta_or_fastq_bool, multiMapSeg_maxLength);							
							}
							else
							{
								if(allUnpairAlignmentCompleteBool_final)
									peAlignSamVec_complete_unpair[tmpOpenMP] 
										= peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_allCases(
											peReadInfo, fasta_or_fastq_bool, multiMapSeg_maxLength);
								else
									peAlignSamVec_incomplete_unpair[tmpOpenMP] 
										= peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_allCases(
											peReadInfo, fasta_or_fastq_bool, multiMapSeg_maxLength);																
							}
						}
					}
				}
				#ifdef MAP_INFO
				cout << "start to do memory free ..." << endl;
				#endif
				delete fixHeadTailInfo;
				peAlignInfo->memoryFree();
				delete peAlignInfo;		
			}
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish fixing Head/Tail, turn: " << tmpTurn+1 << endl;// << endl;
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "start to output ... turn: " << tmpTurn+1 << endl;
			// cout << endl << "[" << asctime(local)
			// 	<< "finish fixing Head/Tail, turn: " << tmpTurn+1 << endl;// << endl;
			// cout << endl << "[" << asctime(local)
			// 	<< "start to output ... turn: " << tmpTurn+1 << endl;
			if(SE_or_PE_bool)
			{
				for(int tmp = 0; tmp < realRecordNum; tmp++)
				{
					if(seAlignSamVec_complete[tmp] != "")
						OutputSamFile_fixHeadTail_complete_SE_ofs << seAlignSamVec_complete[tmp] << endl;				
					if(seAlignSamVec_incomplete[tmp] != "")
						OutputSamFile_fixHeadTail_incomplete_SE_ofs << seAlignSamVec_incomplete[tmp] << endl;
					if(seAlignSamVec_lowScore[tmp] != "")
						OutputSamFile_fixHeadTail_lowScore_SE_ofs << seAlignSamVec_lowScore[tmp] << endl;
				}				
			}
			else
			{	
				for(int tmp = 0; tmp < realRecordNum; tmp++)
				{
					if(peAlignSamVec_complete_pair[tmp] != "")
						OutputSamFile_fixHeadTail_complete_pair_ofs << peAlignSamVec_complete_pair[tmp] << endl;
					if(peAlignSamVec_incomplete_pair[tmp] != "")	
						OutputSamFile_fixHeadTail_incomplete_pair_ofs << peAlignSamVec_incomplete_pair[tmp] << endl;
					if(peAlignSamVec_complete_unpair[tmp] != "")
						OutputSamFile_fixHeadTail_complete_unpair_ofs << peAlignSamVec_complete_unpair[tmp] << endl;
					if(peAlignSamVec_incomplete_unpair[tmp] != "")	
						OutputSamFile_fixHeadTail_incomplete_unpair_ofs << peAlignSamVec_incomplete_unpair[tmp] << endl;
					if(peAlignSamVec_pair_lowScore[tmp] != "")
						OutputSamFile_fixHeadTail_pair_lowScore_ofs << peAlignSamVec_pair_lowScore[tmp] << endl;
				}		
			}
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish output, turn: " << tmpTurn+1 << endl << endl;
			// cout << endl << "[" << asctime(local)
			// 	<< "finish output, turn: " << tmpTurn+1 << endl << endl;
		}
		inputUnfixedHeadTailRecord_ifs.close();
	}
	delete(SJ);

	OutputSamFile_fixHeadTail_complete_pair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_ofs.close();
	OutputSamFile_fixHeadTail_complete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_pair_lowScore_ofs.close();

	OutputSamFile_fixHeadTail_complete_SE_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_SE_ofs.close();
	OutputSamFile_fixHeadTail_lowScore_SE_ofs.close();
	
	nowtime = time(NULL);
	local = localtime(&nowtime);

	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	runtime_log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  

	if(SE_or_PE_bool)
	{
		//statsInfo->getPhase1Stats_SE();
		statsInfo->getFixHeadTailStats_SE();
		statsInfo->outputAllStats_SE_fixHeadTail(stats_ofs, readTotalNum);
		statsInfo->outputFinalStats_SE(stats_ofs, readTotalNum);
	}	
	else
	{	
		statsInfo->getPhase1Stats();
		statsInfo->getFixUnpairedStats();
		statsInfo->getFixHeadTailStats();
		//statsInfo->outputAllStats(log_ofs, readTotalNum);
		statsInfo->outputAllStats(stats_ofs, Do_Phase1_Only, readTotalNum);

		statsInfo->outputFinalStats(stats_ofs, Do_Phase1_Only, readTotalNum);	
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);

	cout << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  
	runtime_log_ofs << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  

	string finalOutputSam = outputDirStr + "/output.sam";
	if(Do_Phase1_Only)
	{
		string cat_cmd;
		if(SE_or_PE_bool)
		{
			cat_cmd = "cat "
				+ tmpHeadSectionInfo
				+ " " + tmpAlignCompleteRead_SE
				+ " " + tmpAlignIncomplete_SE_SAM
				+ " " + tmpAlignUnmapped_mappedToRepeatRegionFile_SE
				+ " " + tmpAlignUnmapped_lowScore_SE
				+ " " + tmpAlignUnmapped_SE
				+ " > " + finalOutputSam;			
		}
		else	
		{
			cat_cmd = "cat "
				+ tmpHeadSectionInfo
				+ " " + tmpAlignCompleteRead
				+ " " + tmpAlignIncompletePair_SAM 
				+ " " + tmpAlignOneEndUnmapped
				+ " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile
				+ " " + tmpAlignBothEndsUnmapped_lowScore
				+ " " + tmpAlignBothEndsUnmapped
				+ " > " + finalOutputSam;
		}
		system(cat_cmd.c_str()); 
	}
	else
	{
		cout << "start to merge all alignment files into one...." << endl;
		log_ofs << "start to merge all alignment files into one...." << endl;
		string cat_cmd;
		if(!SE_or_PE_bool)
		{	
			cat_cmd = "cat " 
				+ tmpHeadSectionInfo
				+ " " + tmpAlignCompleteRead
				+ " " + OutputSamFile_oneEndMapped 
				+ " " + OutputSamFile_fixHeadTail_complete_pair
				+ " " + OutputSamFile_fixHeadTail_incomplete_pair
				+ " " + OutputSamFile_fixHeadTail_complete_unpair
				+ " " + OutputSamFile_fixHeadTail_incomplete_unpair
				+ " " + OutputSamFile_oneEndMapped_unpair
				+ " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile
				+ " " + tmpAlignBothEndsUnmapped_lowScore
				+ " " + OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore
				+ " " + OutputSamFile_fixHeadTail_pair_lowScore
				+ " " + tmpAlignBothEndsUnmapped 
				+ " > " + finalOutputSam;
		}
		else
		{
			cat_cmd = "cat "
				+ tmpHeadSectionInfo
				+ " " + tmpAlignCompleteRead_SE
				+ " " + OutputSamFile_fixHeadTail_complete_SE
				+ " " + OutputSamFile_fixHeadTail_incomplete_SE
				+ " " + tmpAlignUnmapped_mappedToRepeatRegionFile_SE
				+ " " + tmpAlignUnmapped_lowScore_SE
				+ " " + OutputSamFile_fixHeadTail_lowScore_SE
				+ " " + tmpAlignUnmapped_SE
				+ " > " + finalOutputSam;
		}
		system(cat_cmd.c_str()); 	
	}


	if(removeAllIntermediateFilesBool)
	{
		remove(InputSpliceJunction.c_str());
		//remove(juncInsFile.c_str());
		remove(tmpAlignIncompletePair.c_str());
		remove(OutputSamFile_fixHeadTail_complete_pair.c_str());
		remove(OutputSamFile_fixHeadTail_incomplete_pair.c_str());
		remove(OutputSamFile_fixHeadTail_complete_unpair.c_str());
		remove(OutputSamFile_fixHeadTail_incomplete_unpair.c_str());		
		remove(OutputSamFile_oneEndMapped.c_str());
		remove(OutputSamFile_oneEndMapped_unpair.c_str());
		remove(tmpAlignOneEndUnmapped.c_str());
		remove(tmpAlignIncompletePair_SAM.c_str());
		remove(tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());
		remove(tmpAlignBothEndsUnmapped.c_str());
		remove(tmpAlignIncompletePair.c_str());
		//remove(tmpAlignCompleteRead_alignInfo.c_str());
		remove(OutputSamFile_oneEndMapped_alignInfo.c_str());
		remove(tmpAlignCompleteRead.c_str());
		remove(tmpHeadSectionInfo.c_str());

		remove(tmpAlignBothEndsUnmapped_lowScore.c_str());
		remove(OutputSamFile_fixHeadTail_pair_lowScore.c_str());
		remove(OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore.c_str());

		// SE reads alignments
		remove(tmpAlignCompleteRead_SE.c_str());
		remove(tmpAlignIncomplete_SE_SAM.c_str());
		remove(tmpAlignUnmapped_mappedToRepeatRegionFile_SE.c_str());
		remove(tmpAlignUnmapped_lowScore_SE.c_str());
		remove(tmpAlignUnmapped_SE.c_str());
		remove(tmpAlignIncomplete_SE.c_str());	
	}

	annotation_ifs.close();
	delete alignInferJunctionHashInfo;
	delete annotationInfo;	
	delete indexInfo;

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  
	runtime_log_ofs << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl ;  

    return 0;
} //end main