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
#include "phase1/arrayQueue_phase1.h"
#include "phase2/arrayQueue_phase2.h"
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
#include "phase1/phase1_parallelProcesses.h"
#include "phase2/fixOneEndUnmapped_parallelProcesses.h"
#include "phase2/fixHeadTail_parallelProcesses.h"
//#define PreIndexSize 268435456

using namespace std;  

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
	int SNPlocInSyntheticSNPseq = 101;
    bool checkQualSeqForShortAnchorSeqToTargetMap = false;//true;
    cout << "Attention! checkQualSeqForShortAnchorSeqToTargetMap true or not: " 
    	<< checkQualSeqForShortAnchorSeqToTargetMap << endl;
    bool checkQualSeqForReadSegSeq = false;//true;
	cout << "Attention! checkQualSeqForReadSegSeq true or not: " 
		<< checkQualSeqForReadSegSeq << endl;
	
	#ifdef CAL_TIME
    overall_begin = clock();
    #endif
    /////////////////   get option from command line ////////////////////
    Option_Info* optionInfo = new Option_Info();
    optionInfo->getOpt_long(argc, argv);

    bool extractUnmapAlignment2ReadFile_bool 
    	= optionInfo->return_extractUnmapAlignment2ReadFile_bool();
    string inputSamFilePathStr = optionInfo->returnInputSamFilePath();
    if(extractUnmapAlignment2ReadFile_bool)
    {
    	unmappedAlignemnt2ReadFile(inputSamFilePathStr);
    }

    ////////////////////////////// check SE or PE reads  ///////////////////////////////
    bool SE_or_PE_bool = optionInfo->returnSEorPE_bool();
	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   initiate log files    ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
    string outputDirStr = optionInfo->outputFolder_path; //argv[3];
    string outputDirStr_logs = outputDirStr + "/logs";
    string outputDirStr_logs_phase1 = outputDirStr_logs + "/phase1_log";
    string outputDirStr_logs_phase2 = outputDirStr_logs + "/phase2_log";
    string outputDirStr_logs_phase2_fixOneEndUnmapped = outputDirStr_logs_phase2 + "/fixOneEndUnmapped";
    string outputDirStr_logs_phase2_fixHeadTail = outputDirStr_logs_phase2 + "/fixHeadTail";

   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());	
   	string mkdirOutputCommand_log = "mkdir -p " + outputDirStr_logs;
   	system(mkdirOutputCommand_log.c_str());
   	string mkdirOutputCommand_log_phase1 = "mkdir -p " + outputDirStr_logs_phase1;
   	system(mkdirOutputCommand_log_phase1.c_str());   	
   	string mkdirOutputCommand_log_phase2 = "mkdir -p " + outputDirStr_logs_phase2;
   	system(mkdirOutputCommand_log_phase2.c_str());
   	string mkdirOutputCommand_log_phase2_fixOneEndUnmapped = "mkdir -p " + outputDirStr_logs_phase2_fixOneEndUnmapped;
   	system(mkdirOutputCommand_log_phase2_fixOneEndUnmapped.c_str());
   	string mkdirOutputCommand_log_phase2_fixHeadTail = "mkdir -p " + outputDirStr_logs_phase2_fixHeadTail;
   	system(mkdirOutputCommand_log_phase2_fixHeadTail.c_str());

   	string settingsLogStr = outputDirStr_logs + "/settings.log";
   	ofstream settings_log_ofs(settingsLogStr.c_str());
   	string progressLogStr = outputDirStr_logs + "/process.log";
   	ofstream log_ofs(progressLogStr.c_str());
   	string runtimeLogStr = outputDirStr_logs + "/runtime.log";
   	ofstream runtime_log_ofs(runtimeLogStr.c_str());
   	string statsStr = outputDirStr + "/stats.txt";
   	ofstream stats_ofs(statsStr.c_str());

   	string inputLogStr_phase1 = outputDirStr_logs_phase1 + "/input.log";
   	ofstream input_log_ofs_phase1(inputLogStr_phase1.c_str());
   	string outputLogStr_phase1 = outputDirStr_logs_phase1 + "/output.log";
   	ofstream output_log_ofs_phase1(outputLogStr_phase1.c_str());
   	string mappingLogStr_phase1 = outputDirStr_logs_phase1 + "/mapping.log";
   	ofstream mapping_log_ofs_phase1(mappingLogStr_phase1.c_str());   	
   	
   	string inputLogStr_phase2_fixOneEndUnmapped = outputDirStr_logs_phase2_fixOneEndUnmapped + "/input.log";
   	ofstream input_log_ofs_phase2_fixOneEndUnmapped(inputLogStr_phase2_fixOneEndUnmapped.c_str());
   	string outputLogStr_phase2_fixOneEndUnmapped = outputDirStr_logs_phase2_fixOneEndUnmapped + "/output.log";
   	ofstream output_log_ofs_phase2_fixOneEndUnmapped(outputLogStr_phase2_fixOneEndUnmapped.c_str());
   	string mappingLogStr_phase2_fixOneEndUnmapped = outputDirStr_logs_phase2_fixOneEndUnmapped + "/mapping.log";
   	ofstream mapping_log_ofs_phase2_fixOneEndUnmapped(mappingLogStr_phase2_fixOneEndUnmapped.c_str());   

   	string inputLogStr_phase2_fixHeadTail = outputDirStr_logs_phase2_fixHeadTail + "/input.log";
   	ofstream input_log_ofs_phase2_fixHeadTail(inputLogStr_phase2_fixHeadTail.c_str());
   	string outputLogStr_phase2_fixHeadTail = outputDirStr_logs_phase2_fixHeadTail + "/output.log";
   	ofstream output_log_ofs_phase2_fixHeadTail(outputLogStr_phase2_fixHeadTail.c_str());
   	string mappingLogStr_phase2_fixHeadTail = outputDirStr_logs_phase2_fixHeadTail + "/mapping.log";
   	ofstream mapping_log_ofs_phase2_fixHeadTail(mappingLogStr_phase2_fixHeadTail.c_str());


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

	int normalRecordNum_1stMapping = 5000;//000;
	int normalRecordNum_fixOneEndUnmapped = 5000;//000;
	int normalRecordNum_fixHeadTail = 5000;//000;

	//int readTotalNum = 0;

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
	//////////////////////////////////////////////        LOAD INDEX         //////////////////////////////////////////
    string InputReadFile = optionInfo->read_file_path_1;//read sample, exacted from fastq file every time
    if(SE_or_PE_bool)
    	InputReadFile = optionInfo->read_file_path_SE;
    string InputReadFile_PE = optionInfo->read_file_path_2;// another end read for pair-end reads

	int threads_num = optionInfo->threads_num;//atoi(threadsNumStr.c_str());
	#ifdef MAP_INFO
	threads_num = 1;
	#endif

	bool InputAsFastq = (!(optionInfo->fasta_or_fastq_bool));
	bool fasta_or_fastq_bool = optionInfo->fasta_or_fastq_bool;
	if(checkQualSeqForShortAnchorSeqToTargetMap && fasta_or_fastq_bool)
	{
		cout << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		log_ofs << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		exit(1);
	}
	/////////////////////////////////////          LOAD INDEX         ////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
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
	int *lcpCompress_test; lcpCompress_test = (int*)malloc((indexInfo->returnIndexSize()) * sizeof(int)); 	
	int *verifyChild_test; verifyChild_test = (int*)malloc((indexInfo->returnIndexSize()) * sizeof(int)); 	

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
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
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
		alignInferJunctionHashInfo->insertJuncFromJuncFile(tmpInputJuncFile, indexInfo);
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
		SNP_seq_index_path = optionInfo->SNP_seq_index_path;
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

	string mkdirOutputCommand_phase1 = "mkdir -p " + outputDirStr + "/phase1_output";
	system(mkdirOutputCommand_phase1.c_str());

	string mkdirOutputCommand_repeatRegionFile = mkdirOutputCommand_phase1 + "/repeat_region";
   	system(mkdirOutputCommand_repeatRegionFile.c_str());
	string repeatRegionFile = outputDirStr + "/phase1_output/repeat_region/repeatRegion";
	ofstream repeatRegionFile_ofs(repeatRegionFile.c_str());

	string mkdirOutputCommand_tmpAlignCompleteRead = mkdirOutputCommand_phase1 + "/completePair";
	system(mkdirOutputCommand_tmpAlignCompleteRead.c_str());
	string tmpAlignCompleteRead_SE = outputDirStr + "/phase1_output/completePair/complete_SE.sam";
	ofstream tmpAlignCompleteRead_SE_ofs(tmpAlignCompleteRead_SE.c_str());
	string tmpAlignCompleteRead = outputDirStr + "/phase1_output/completePair/completePair.sam";
	ofstream tmpAlignCompleteRead_ofs(tmpAlignCompleteRead.c_str());

	string mkdirOutputCommand_tmpAlignOneEndUnmapped = mkdirOutputCommand_phase1 + "/oneEndUnmapped";
	system(mkdirOutputCommand_tmpAlignOneEndUnmapped.c_str());
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
		inputRead_PE_ifs.close();

    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE;
    string line2_afterProcess, line2_PE_afterProcess;

	//int normalRecordNum = normalRecordNum_1stMapping; //1000000;//1500000;

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
	Read_Array_Queue* readArrayQueue = new Read_Array_Queue();
	Result_Array_Queue* resultArrayQueue = new Result_Array_Queue();
	
	int tmpInputReadNumInBatchArray_phase1 = normalRecordNum_1stMapping;//ReadNumInReadArray_Phase1;
	int tmpInputTimeWeight_phase1 = InputTimeWeight_Phase1;
	int tmpOutputTimeWeigth_phase1 = OutputTimeWeight_Phase1; 
	bool endOfFile_bool = false;
	bool endOfProcessing_bool = false;

	omp_set_num_threads(2);
	omp_set_nested(1);
#pragma omp parallel
	{
#pragma omp sections
		{
#pragma omp section
			io_stage_phase1(inputRead_ifs, inputRead_PE_ifs, readArrayQueue, resultArrayQueue, endOfFile_bool, endOfProcessing_bool, 
				tmpInputReadNumInBatchArray_phase1, tmpInputTimeWeight_phase1, tmpOutputTimeWeigth_phase1, log_ofs, readPreProcessInfo,
				tmpAlignCompleteRead_ofs, tmpAlignIncompletePair_ofs, tmpAlignOneEndUnmapped_ofs, tmpAlignBothEndsUnmapped_ofs,
				tmpAlignBothEndsUnmapped_lowScore_ofs, tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
				tmpAlignIncompletePair_SAM_ofs, input_log_ofs_phase1, output_log_ofs_phase1, fasta_or_fastq_bool, SE_or_PE_bool);
#pragma omp section
			process_stage_phase1_main(readArrayQueue, resultArrayQueue, endOfFile_bool, endOfProcessing_bool, threads_num-1, 
				sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, preIndexMapLengthArray, preIndexIntervalStartArray,
				preIndexIntervalEndArray, repeatRegionInfoVec, Do_cirRNA, Do_extendHeadTail_phase1, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, outputDirectlyBool_Phase1Only, Do_Phase1_Only,	statsInfo, 
				fasta_or_fastq_bool, mapping_log_ofs_phase1, checkQualSeqForReadSegSeq, SE_or_PE_bool, sa_SNP, lcpCompress_SNP,
				childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, SNPlocInSyntheticSNPseq);
		}
	}	



	//log_ofs << "perfectMatch_pair #: " << perfectMatch_pair << endl;
	repeatRegionFile_ofs << "Repeat Region Info: size = " << repeatRegionInfoVec.size() << endl;
	for(int tmpThread = 0; tmpThread < threads_num; tmpThread++)
		repeatRegionInfoVec[tmpThread]->outputRepeatRegion(tmpThread+1, indexInfo, sa, 100, repeatRegionFile_ofs);

	int readTotalNum = 0;
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
			for(int tmpSecondLevelIndexNO = 1; tmpSecondLevelIndexNO 
				<= (indexInfo->returnSecondLevelIndexPartsNum(tmpChrNO)); tmpSecondLevelIndexNO ++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpSecondLevelIndexNO);
				string tmpFileNumStr = tmpFileNumChar;				
				string inputIndexFileStr = secondLevelIndexStr + "/" + indexInfo->returnChrNameStr(tmpChrNO) + "/" + tmpFileNumStr + "/";
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
					tmpSecondLevelChrom[tmpMallocSpace] = '0';

				secondLevelChrom_file_ifs.read((char*)tmpSecondLevelChrom, sizeOfIndex * sizeof(char));
				if(tmpSecondLevelChrom[sizeOfIndex-1] != 'X')
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);

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
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
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
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On one end unmapped Reads    ///////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 	 
	log_ofs << endl << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 
	runtime_log_ofs << endl << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads starts ......" << endl << endl ; 


	string mkdirOutputCommand_phase2 = "mkdir -p " + outputDirStr + "/phase2_output";
	system(mkdirOutputCommand_phase2.c_str());
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

		AlignInfoInput_Array_Queue* alignInfoInputQueue = new AlignInfoInput_Array_Queue();
		Result_FixOneEndUnmapped_Array_Queue* fixOneEndUnmappedResultQueue = new Result_FixOneEndUnmapped_Array_Queue();

		string oneEndMappedFileStr = tmpAlignOneEndUnmapped;
		ifstream inputRecord_ifs(oneEndMappedFileStr.c_str());		
		bool endOfFile_bool = false;
		bool endOfProcessing_bool = false;

		int tmpInputReadNumInBatchArray_fixOneEndUnmapped = normalRecordNum_fixOneEndUnmapped;//ReadNumInReadArray_FixOneEndUnmapped;
		int tmpInputTimeWeight_fixOneEndUnmapped = InputTimeWeight_FixOneEndUnmapped;
		int tmpOutputTimeWeight_fixOneEndUnmapped = OutputTimeWeight_FixOneEndUnmapped;

		omp_set_num_threads(2);
		omp_set_nested(1);
#pragma omp parallel
		{
#pragma omp sections
			{
#pragma omp section
				io_stage_fixOneEndUnmapped(
					inputRecord_ifs,
					alignInfoInputQueue,
					fixOneEndUnmappedResultQueue,
					endOfFile_bool,
					endOfProcessing_bool,
					tmpInputReadNumInBatchArray_fixOneEndUnmapped,
					tmpInputTimeWeight_fixOneEndUnmapped,
					tmpOutputTimeWeight_fixOneEndUnmapped,					

					log_ofs,

					OutputSamFile_oneEndMapped_ofs,
					tmpAlignIncompletePair_ofs,
					OutputSamFile_oneEndMapped_unpair_ofs,
					OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs,

					input_log_ofs_phase2_fixOneEndUnmapped,
					output_log_ofs_phase2_fixOneEndUnmapped
					);
#pragma omp section
				//#ifdef PERSONALIZED_CHR_SEQ
				process_stage_fixOneEndUnmapped_main(
					alignInfoInputQueue,
					fixOneEndUnmappedResultQueue,
					endOfFile_bool,
					endOfProcessing_bool,
					threads_num-1,
					fasta_or_fastq_bool,
					statsInfo,
					secondLevelChrom,
					secondLevelSa,
					secondLevelLcpCompress,
					secondLevelChildTab,
					secondLevelDetChild,					
					indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE2,
					checkQualSeqForReadSegSeq,
					mapping_log_ofs_phase2_fixOneEndUnmapped,
					SE_or_PE_bool,
					sa_SNP, 
					lcpCompress_SNP,
					childTab_SNP, 
					chrom_SNP,
					verifyChild_SNP, 
					indexInfo_SNP,
					SNPlocInSyntheticSNPseq);
			}
		}
		alignInfoInputQueue->free();
		delete alignInfoInputQueue;
		fixOneEndUnmappedResultQueue->free();
		delete fixOneEndUnmappedResultQueue;
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
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	

	////////  take annotation as reference for remapping //////////////
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
		char entry[500];
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
			alignInferJunctionHashInfo->convert2SJhashInfo(SJ, indexInfo);
		}
		if(annotation_provided_bool)
		{	
			// loading SJs in annotation file (if provided)
			FILE *fp_annotatedSJ_file = fopen(annotation_file_path.c_str(), "r");
			fgets(entry, sizeof(entry), fp_annotatedSJ_file);
			while(!feof(fp_annotatedSJ_file))
			{
				fgets(entry, sizeof(entry), fp_annotatedSJ_file);
				if(feof(fp_annotatedSJ_file))
					break;
				junctionNum ++;
				entryString = entry;
				tabLocation1 = entryString.find('\t', 0);
				tabLocation2 = entryString.find('\t', tabLocation1+1);
				tabLocation3 = entryString.find('\t', tabLocation2+1);
				chrIntString = entryString.substr(0, tabLocation1);
				spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
				spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
				chrInt = indexInfo->convertStringToInt(chrIntString);
				spliceStartPos = atoi(spliceStartPosString.c_str());
				spliceEndPos = atoi(spliceEndPosString.c_str());	
				SJ->insert2AreaAndStringHash(chrInt, spliceStartPos, spliceEndPos, indexInfo);
			}
			fclose(fp_annotatedSJ_file);
		}
		if(junctionNum == 0)
		{
			spliceJunctionHashExists = false;
		}

		settings_log_ofs << "finish building spliceJunction Hash" << endl;
		settings_log_ofs << "After inserting SJs generated from alignments and annotation, junctionNum = " 
			<< junctionNum << endl;
		settings_log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		log_ofs << "finish building spliceJunction Hash" << endl;		
		log_ofs << "After inserting SJs generated from alignments and annotation, junctionNum = " 
			<< junctionNum << endl;
		log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		cout << "finish building spliceJunction Hash" << endl;		
		cout << "After inserting SJs generated from alignments and annotation, junctionNum = " 
			<< junctionNum << endl;
		cout << "start doing remapping on unfixed head/tail alignments" << endl;		
		
		string headTailSoftClippingFile;
		if(SE_or_PE_bool)
			headTailSoftClippingFile = tmpAlignIncomplete_SE;
		else
			headTailSoftClippingFile = tmpAlignIncompletePair;// + ".all";

		AlignInfoInput_Array_Queue* fixHeadTailAlignInfoInputQueue = new AlignInfoInput_Array_Queue();
		Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue = new Result_FixHeadTail_Array_Queue();
		if(outputDirectlyBool_Phase1Only)
			headTailSoftClippingFile += ".all";
		ifstream inputUnfixedHeadTailRecord_ifs(headTailSoftClippingFile.c_str());
		//int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;
		bool endOfFile_bool = false;
		bool endOfProcessing_bool = false;

		int tmpInputReadNumInBatchArray_fixHeadTail = normalRecordNum_fixHeadTail;//ReadNumInReadArray_FixHeadTail;
		int tmpInputTimeWeight_fixHeadTail = InputTimeWeight_FixHeadTail;
		int tmpOutputTimeWeight_fixHeadTail = OutputTimeWeight_FixHeadTail;

		omp_set_num_threads(2);
		omp_set_nested(1);
#pragma omp parallel
		{
#pragma omp sections
			{
#pragma omp section
				io_stage_fixHeadTail(
					inputUnfixedHeadTailRecord_ifs,
					fixHeadTailAlignInfoInputQueue,
					fixHeadTailResultQueue,
					endOfFile_bool,
					endOfProcessing_bool,
					tmpInputReadNumInBatchArray_fixHeadTail,
					tmpInputTimeWeight_fixHeadTail,
					tmpOutputTimeWeight_fixHeadTail,
					
					log_ofs,

					OutputSamFile_fixHeadTail_complete_pair_ofs,
					OutputSamFile_fixHeadTail_incomplete_pair_ofs,
					OutputSamFile_fixHeadTail_complete_unpair_ofs,
					OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
					OutputSamFile_fixHeadTail_pair_lowScore_ofs,					

					input_log_ofs_phase2_fixHeadTail,
					output_log_ofs_phase2_fixHeadTail
					);
#pragma omp section
				//#ifdef PERSONALIZED_CHR_SEQ
				process_stage_fixHeadTail_main(
					fixHeadTailAlignInfoInputQueue,
					fixHeadTailResultQueue,
					endOfFile_bool,
					endOfProcessing_bool,
					threads_num-1,

					fasta_or_fastq_bool,
					statsInfo,
					secondLevelChrom,
					secondLevelSa,
					secondLevelLcpCompress,
					secondLevelChildTab,
					secondLevelDetChild,
					indexInfo,
					SJ, 
					Do_extendHeadTail_fixHeadTail,
					annotation_provided_bool, 
					Do_annotation_only_bool, 
					annotationInfo,
					checkQualSeqForReadSegSeq,
					checkQualSeqForShortAnchorSeqToTargetMap,
					spliceJunctionHashExists,

					Do_fixHeadTail_remapping,
					Do_fixHeadTail_greedyMapping,
					Do_fixHeadTail_remappingAndTargetMapping,
					Do_fixHeadTail_remappingAgain,

					mapping_log_ofs_phase2_fixHeadTail,
					SE_or_PE_bool,
			
					sa_SNP, 
					lcpCompress_SNP,
					childTab_SNP, 
					chrom_SNP,
					verifyChild_SNP, 
					indexInfo_SNP,
					SNPlocInSyntheticSNPseq);
				// #else
				// process_stage_fixHeadTail(
				// 	fixHeadTailAlignInfoInputQueue,
				// 	fixHeadTailResultQueue,
				// 	endOfFile_bool,
				// 	endOfProcessing_bool,
				// 	threads_num-1,

				// 	fasta_or_fastq_bool,
				// 	statsInfo,
				// 	secondLevelChrom,
				// 	secondLevelSa,
				// 	secondLevelLcpCompress,
				// 	secondLevelChildTab,
				// 	secondLevelDetChild,
				// 	indexInfo,
				// 	SJ, 
				// 	Do_extendHeadTail_fixHeadTail,
				// 	annotation_provided_bool, 
				// 	Do_annotation_only_bool, 
				// 	annotationInfo,
				// 	checkQualSeqForReadSegSeq,
				// 	checkQualSeqForShortAnchorSeqToTargetMap,
				// 	spliceJunctionHashExists,

				// 	Do_fixHeadTail_remapping,
				// 	Do_fixHeadTail_greedyMapping,
				// 	Do_fixHeadTail_remappingAndTargetMapping,
				// 	Do_fixHeadTail_remappingAgain,

				// 	mapping_log_ofs_phase2_fixHeadTail,
				// 	SE_or_PE_bool
				// 	);
				// #endif
			}
		}
		fixHeadTailAlignInfoInputQueue->free();
		delete fixHeadTailAlignInfoInputQueue;
		fixHeadTailResultQueue->free();
		delete fixHeadTailResultQueue;
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
		statsInfo->getFixHeadTailStats_SE();
		statsInfo->outputAllStats_SE_fixHeadTail(stats_ofs, readTotalNum);
		statsInfo->outputFinalStats_SE(stats_ofs, readTotalNum);
	}	
	else
	{	
		statsInfo->getPhase1Stats();
		statsInfo->getFixUnpairedStats();
		statsInfo->getFixHeadTailStats();
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