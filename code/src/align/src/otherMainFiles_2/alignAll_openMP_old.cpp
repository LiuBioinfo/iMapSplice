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

#ifdef CAL_TIME
clock_t read_file_begin, read_file_end, read_file_end2, read_file_end3, 
		align_begin, align_end, align_cost = 0,
    	overall_begin, overall_end, overall_cost = 0,
    	input_begin, input_end, input_cost = 0,
    	output_begin, output_end, output_cost = 0,
    	segMap_begin, segMap_end, segMap_cost = 0,
    	getPath_begin, getPath_end, getPath_cost = 0,
    	fixGap_begin, fixGap_end, fixGap_cost = 0,
    	
    	fixGap_1_begin, fixGap_1_end, fixGap_1_cost = 0,
    	fixGap_2_begin, fixGap_2_end, fixGap_2_cost = 0,

    	getChrLocation_begin, getChrLocation_end, getChrLocation_cost = 0,

    	getInterval_begin, getInterval_end, getInterval_cost = 0,
    	getFirstInterval_begin, getFirstInterval_end, getFirstInterval_cost = 0,   
    	getLcp_begin_1, getLcp_end_1, getLcp_cost_1 = 0,
    	getLcp_begin_2, getLcp_end_2, getLcp_cost_2 = 0,
       	getLcp_begin_3, getLcp_end_3, getLcp_cost_3 = 0,
       	 	
       	score_string_begin, score_string_end, score_string_cost = 0,
       	checkTwoStringMatch_begin, checkTwoStringMatch_end, checkTwoStringMatch_cost = 0,

    	searchPrefix_begin, searchPrefix_end, searchPrefix_cost = 0,
    	searchInSA_begin, searchInSA_end, searchInSA_cost = 0,
    	insert2SegInfo_begin, insert2SegInfo_end, insert2SegInfo_cost = 0,

    	readFromFile_begin, readFromFile_end, readFromFile_cost = 0,
    	readPreProcess_begin, readPreProcess_end, readPreProcess_cost = 0,

    	getPEalignInfo_begin, getPEalignInfo_end, getPEalignInfo_cost = 0,
    	selectBestAlign_begin, selectBestAlign_end, selectBestAlign_cost = 0,
    	getSamFormat_begin, getSamFormat_end, getSamFormat_cost = 0,
		insertSamFormat_begin, insertSamFormat_end, insertSamFormat_cost = 0,
		getReadInfo_begin, getReadInfo_end, getReadInfo_cost = 0,
		freeMem_begin, freeMem_end, freeMem_cost = 0;
#endif 
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
#include "general/path_info.h"
#include "general/gap_info.h"
#include "general/align_info.h"
#include "general/peAlign_info.h"
#include "general/alignInferJunctionHash_info.h"
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
#include "general/alignmentToJunc.h"
#include "phase1/phase1_parallelProcesses.h"
#include "phase2/fixOneEndUnmapped_parallelProcesses.h"
#include "phase2/fixHeadTail_parallelProcesses.h"

using namespace std;  


int main(int argc, char**argv)
{
    bool checkQualSeqForShortAnchorSeqToTargetMap = false;
    cout << "Attention! checkQualSeqForShortAnchorSeqToTargetMap true or not: " 
    	<< checkQualSeqForShortAnchorSeqToTargetMap << endl;
    bool checkQualSeqForReadSegSeq = false;
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
	//////////////////////////////////////////////////
    bool SE_or_PE_bool = optionInfo->returnSEorPE_bool();

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

	//int normalRecordNum_1stMapping = 2000000;
	//int normalRecordNum_fixOneEndUnmapped = 2000000;
	//int normalRecordNum_fixHeadTail = 2000000;

	//int readTotalNum = 0;

	optionInfo->outputSwitchInfo(Do_Phase1_Only, outputAlignInfoAndSamForAllPairedAlignmentBool,
		removeAllIntermediateFilesBool, Do_cirRNA, outputDirectlyBool_Phase1Only, 
		//normalRecordNum_1stMapping, normalRecordNum_fixOneEndUnmapped,
		//normalRecordNum_fixHeadTail, 
		Do_extendHeadTail_phase1, 
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
    string InputReadFile_PE = optionInfo->read_file_path_2;// another end read for pair-end reads

	//string threadsNumStr = argv[4];
	int threads_num = optionInfo->threads_num;//atoi(threadsNumStr.c_str());

	#ifdef CAL_TIME	
	threads_num = 1;
	#endif

	#ifdef MAP_INFO
	threads_num = 1;
	#endif

	bool InputAsFastq = (!(optionInfo->fasta_or_fastq_bool));
	bool fasta_or_fastq_bool = optionInfo->fasta_or_fastq_bool;
	// if((checkQualSeqForShortAnchorSeqToTargetMap && fasta_or_fastq_bool)||(checkQualSeqForReadSegSeq && fasta_or_fastq_bool))
	// {
	// 	cout << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
	// 	log_ofs << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
	// 	cout << "if checkQualSeqForReadSegSeq, fasta_or_fastq_bool must be true" << endl;
	// 	log_ofs << "if checkQualSeqForReadSegSeq, fasta_or_fastq_bool must be true" << endl;
	// 	exit(1);
	// }
	if(checkQualSeqForShortAnchorSeqToTargetMap && fasta_or_fastq_bool)
	{
		cout << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		log_ofs << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		exit(1);
	}	
	/////////////////////////////////////          LOAD INDEX         ////////////////////////////////////////////////////////////////
	#ifdef CAL_TIME
	read_file_begin = clock();
	#endif
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
	/////////////////////////////////////   start to load annotation  /////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load annotation file (SJs)......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load annotation file (SJs) ......" << endl << endl; 	
	string annotation_file_path = optionInfo->annotation_file_path; // junction files
	ifstream annotation_ifs(annotation_file_path.c_str());
	Annotation_Info* annotationInfo = new Annotation_Info();
	if(annotation_provided_bool)
	{
		//annotation_ifs.open(annotation_file_path);
		annotationInfo->initiateAndReadAnnotationFile(indexInfo, annotation_ifs);
	}
	/////////////////////////////////////   finish loading annotation  /////////////////////////////////////
 	
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

	#ifdef CAL_TIME
	read_file_end = clock();
	double read_file_time = (double)(read_file_end - read_file_begin)/CLOCKS_PER_SEC;
	log_ofs << "read_file cpu time = " << read_file_time << endl;
	#endif
	//////////////////////////////////////////////////
	HeaderSection_Info* headerSectionInfo = new HeaderSection_Info(indexInfo);	
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
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
		// insert SJ into SJmap
		string tmpInputJuncFile = optionInfo->spliceJunctionAlignInferHash_file_path;
		//vector<string> inputJUNCfileVec;
		//inputJUNCfileVec.push_back(tmpInputJuncFile);
		alignInferJunctionHashInfo->insertJuncFromJuncFile(tmpInputJuncFile, indexInfo);
		cout << "start to output SJ map" << endl;
	}
	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       finish LOADing INDEX           ////////////////////////////////////////////////////////////////
	
	Stats_Info* statsInfo = new Stats_Info();
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
	string tmpAlignCompleteRead = outputDirStr + "/phase1_output/completePair/completePair.sam";
	ofstream tmpAlignCompleteRead_ofs(tmpAlignCompleteRead.c_str());
	string tmpAlignCompleteRead_alignInfo = outputDirStr + "/phase1_output/completePair/completePair.sam_alignInfo";
	ofstream tmpAlignCompleteRead_alignInfo_ofs(tmpAlignCompleteRead_alignInfo.c_str());

	#ifdef CHECK_MULTI
	string tmpAlignCompleteRead_multi = outputDirStr + "/phase1_output/completePair/completePair.sam.multi";
	ofstream tmpAlignCompleteRead_multi_ofs(tmpAlignCompleteRead_multi.c_str());
	string tmpAlignCompleteRead_unique = outputDirStr + "/phase1_output/completePair/completePair.sam.unique";
	ofstream tmpAlignCompleteRead_unique_ofs(tmpAlignCompleteRead_unique.c_str());
	#endif

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
	string tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile 
		= outputDirStr + "/phase1_output/repeat_region/bothEndsUnmapped_mappedToRepeatRegion.sam"; 
	ofstream tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs(
		tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());

	string mkdirOutputCommand_tmpAlignBothEndsUnmapped = mkdirOutputCommand_phase1 + "/bothEndsUnmapped";
	system(mkdirOutputCommand_tmpAlignBothEndsUnmapped.c_str());
	string tmpAlignBothEndsUnmapped = outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped.sam";
	ofstream tmpAlignBothEndsUnmapped_ofs(tmpAlignBothEndsUnmapped.c_str());
	string tmpAlignBothEndsUnmapped_lowScore = outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped_lowScore.sam";
	ofstream tmpAlignBothEndsUnmapped_lowScore_ofs(tmpAlignBothEndsUnmapped_lowScore.c_str());

	string mkdirOutputCommand_tmpAlignIncompletePair = mkdirOutputCommand_phase1 + "/incomplete";
	system(mkdirOutputCommand_tmpAlignIncompletePair.c_str());
	string tmpAlignIncompletePair = outputDirStr + "/phase1_output/incomplete/incomplete.alignInfo"; 
	ofstream tmpAlignIncompletePair_ofs(tmpAlignIncompletePair.c_str());
	string tmpAlignIncompletePair_SAM = outputDirStr + "/phase1_output/incomplete/incompletePair.sam"; 
	ofstream tmpAlignIncompletePair_SAM_ofs(tmpAlignIncompletePair_SAM.c_str());	

	string tmpIntermediateJunctionFile = outputDirStr + "/phase2_output/inter.junc";

	ifstream inputRead_ifs(InputReadFile.c_str());
	ifstream inputRead_PE_ifs(InputReadFile_PE.c_str());
	if(SE_or_PE_bool)
		inputRead_PE_ifs.close();

	vector< RepeatRegion_Info* > repeatRegionInfoVec;

	for(int tmp = 0; tmp < threads_num-1; tmp++)
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
	
	int tmpInputReadNumInBatchArray_phase1 = ReadNumInReadArray_Phase1;
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
			io_stage_phase1(inputRead_ifs,
				inputRead_PE_ifs,
				readArrayQueue,
				resultArrayQueue,
				endOfFile_bool,
				endOfProcessing_bool,
				tmpInputReadNumInBatchArray_phase1,
				tmpInputTimeWeight_phase1,
				tmpOutputTimeWeigth_phase1,

				log_ofs,
				readPreProcessInfo,
				tmpAlignCompleteRead_ofs,
				tmpAlignIncompletePair_ofs,
				tmpAlignOneEndUnmapped_ofs,
				tmpAlignBothEndsUnmapped_ofs,
				tmpAlignBothEndsUnmapped_lowScore_ofs,
				tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
				tmpAlignIncompletePair_SAM_ofs, input_log_ofs_phase1, output_log_ofs_phase1,
				fasta_or_fastq_bool,
				SE_or_PE_bool);
// #pragma omp section
// 			test(2);			
#pragma omp section
			process_stage_phase1(
				readArrayQueue,
				resultArrayQueue,
				endOfFile_bool,
				endOfProcessing_bool,	
				threads_num-1, 
				sa, 
				lcpCompress,
				childTab, 
				chrom,
				verifyChild, 
				indexInfo,
				preIndexMapLengthArray, 
				preIndexIntervalStartArray,
				preIndexIntervalEndArray,
				repeatRegionInfoVec,
				Do_cirRNA, 
				Do_extendHeadTail_phase1,
				annotation_provided_bool, 
				Do_annotation_only_bool, 
				annotationInfo,	
				outputDirectlyBool_Phase1Only,
				Do_Phase1_Only,	
				statsInfo, 
				fasta_or_fastq_bool, mapping_log_ofs_phase1,//, mapping_log_ofs_vec
				checkQualSeqForReadSegSeq,
				SE_or_PE_bool);
	  	}
	}

	//log_ofs << "perfectMatch_pair #: " << perfectMatch_pair << endl;
	repeatRegionFile_ofs << "Repeat Region Info: size = " << repeatRegionInfoVec.size() << endl;
	for(int tmpThread = 0; tmpThread < threads_num-1; tmpThread++)
	{
		repeatRegionInfoVec[tmpThread]->outputRepeatRegion(tmpThread+1, indexInfo, sa, 100, repeatRegionFile_ofs);
		//repeatRegionInfoVec[tmpThread]->outputRepeatRegion(esaInfo, tmpThread+1, indexInfo, //sa, 
		//	100, repeatRegionFile_ofs);
	}

	//settings_log_ofs << "readTotalNum: " << readTotalNum << endl;

	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	tmpAlignCompleteRead_ofs.close();
	tmpAlignOneEndUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs.close();
	tmpAlignIncompletePair_SAM_ofs.close();
	if(Do_Phase1_Only)
	{
		tmpAlignIncompletePair_ofs.close();
	}

	//fclose(fp_in);

	free(preIndexMapLengthArray); 
	free(preIndexIntervalStartArray); 
	free(preIndexIntervalEndArray);
	free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
	free(childTab);
	free(verifyChild);

	free(chrom);
	//esaInfo->freeMem();
	//delete esaInfo;
	
	#ifdef CAL_TIME
	overall_end = clock();
	#endif
	
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... 1st mapping process ends ......" << endl << endl; 

	log_ofs << endl << "**********************************" << endl << "**********************************";
	runtime_log_ofs << endl << "**********************************" << endl << "**********************************";

	#ifdef CAL_TIME
	//if(threads_num == 1)
	//{	
		double overall_time = (double)(overall_end - overall_begin)/CLOCKS_PER_SEC;
		double input_time = (double)input_cost/CLOCKS_PER_SEC;
		double align_time = (double)align_cost/CLOCKS_PER_SEC;
		double ouput_time = (double)output_cost/CLOCKS_PER_SEC;

		double getReadInfo_time = (double)getReadInfo_cost/CLOCKS_PER_SEC;
		double segMap_time = (double)segMap_cost/CLOCKS_PER_SEC;
		double getPath_time = (double)getPath_cost/CLOCKS_PER_SEC;
		double fixGap_time = (double)fixGap_cost/CLOCKS_PER_SEC;

		double fixGap_1_time = (double)fixGap_1_cost/CLOCKS_PER_SEC;
		double fixGap_2_time = (double)fixGap_2_cost/CLOCKS_PER_SEC;

		double getChrLocation_time = (double)getChrLocation_cost/CLOCKS_PER_SEC;

		double getPEalignInfo_time = (double)getPEalignInfo_cost/CLOCKS_PER_SEC;
		double selectBestAlign_time = (double)selectBestAlign_cost/CLOCKS_PER_SEC;
		double getSamFormat_time = (double)getSamFormat_cost/CLOCKS_PER_SEC;
		double freeMem_time = (double)freeMem_cost/CLOCKS_PER_SEC;

		double getInterval_time = (double)getInterval_cost/CLOCKS_PER_SEC;
		double getFirstInterval_time = (double)getFirstInterval_cost/CLOCKS_PER_SEC;
		double getLcp_1_time = (double)getLcp_cost_1/CLOCKS_PER_SEC;
		double getLcp_2_time = (double)getLcp_cost_2/CLOCKS_PER_SEC;
		double getLcp_3_time = (double)getLcp_cost_3/CLOCKS_PER_SEC;

		double searchPrefix_time = (double)searchPrefix_cost/CLOCKS_PER_SEC;
		double searchInSA_time = (double)searchInSA_cost/CLOCKS_PER_SEC;		
		double insert2SegInfo_time = (double)insert2SegInfo_cost/CLOCKS_PER_SEC;

		double readFromFile_time = (double)readFromFile_cost/CLOCKS_PER_SEC;
		double readPreProcess_time = (double)readPreProcess_cost/CLOCKS_PER_SEC;

		double score_string_time = (double)score_string_cost/CLOCKS_PER_SEC;
		double checkTwoStringMatch_time = (double)checkTwoStringMatch_cost/CLOCKS_PER_SEC;

		log_ofs << endl << "overall_time = " << overall_time << endl;
		log_ofs << endl << "input_time = " << input_time << endl;
		log_ofs << endl << "align_time = " << align_time << endl << endl;
		
		log_ofs << endl << "getReadInfo_time = " << getReadInfo_time << endl;
		log_ofs << endl << "segMap_time = " << segMap_time << endl;
		log_ofs << endl << "getPath_time = " << getPath_time << endl;
		log_ofs << endl << "fixGap_time = " << fixGap_time << endl;
		
		log_ofs << endl << "fixGap_1_time = " << fixGap_1_time << endl;
		log_ofs << endl << "fixGap_2_time = " << fixGap_2_time << endl;		

		log_ofs << endl << "getChrLocation_time = " << getChrLocation_time << endl;

		log_ofs << endl << "getPEalignInfo_time = " << getPEalignInfo_time << endl;
		log_ofs << endl << "selectBestAlign_time = " << selectBestAlign_time << endl;
		log_ofs << endl << "getSamFormat_time = " << getSamFormat_time << endl;
		log_ofs << endl << "freeMem_time = " << freeMem_time << endl;
		log_ofs << endl << endl << "ouput_time = " << ouput_time << endl;

		log_ofs << endl << "getInterval_time = " << getInterval_time << endl;
		log_ofs << endl << "getFirstInterval_time = " << getFirstInterval_time << endl;

		log_ofs << endl << "getLcp_time_1 = " << getLcp_1_time << endl;
		log_ofs << endl << "getLcp_time_2 = " << getLcp_2_time << endl;
		log_ofs << endl << "getLcp_time_3 = " << getLcp_3_time << endl;

		log_ofs << endl << "searchPrefix_time = " << searchPrefix_time << endl;
		log_ofs << endl << "searchInSA_time = " << searchInSA_time << endl;
		log_ofs << endl << "insert2SegInfo_time = " << insert2SegInfo_time << endl;

		log_ofs << endl << "readFromFile_time = " << readFromFile_time << endl;
		log_ofs << endl << "readPreProcess_time = " << readPreProcess_time << endl;		

		log_ofs << endl << "score_string_time = " << score_string_time << endl;
		log_ofs << endl << "checkTwoStringMatch_time = " << checkTwoStringMatch_time << endl;

		if(threads_num != 1)
			log_ofs << "!!! set thread=1 to get the time-cost information ... !!!" << endl;
	#endif

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

		int tmpInputReadNumInBatchArray_fixOneEndUnmapped = ReadNumInReadArray_FixOneEndUnmapped;
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
				process_stage_fixOneEndUnmapped(
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
					SE_or_PE_bool);
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
	//tmpAlignInfoForDebugFile_oneEndMapped_ofs.close();

	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl;  
	log_ofs << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... fixing oneEndUnmapped reads ends ......" << endl << endl; 
	
	log_ofs << endl << "**********************************************************************************" << endl << endl;
	runtime_log_ofs << endl << "**********************************************************************************" << endl << endl;


	///////////////////////////////////////  merging incomplete alignment files ///////////////////////////////////////////
	if(outputDirectlyBool_Phase1Only)
	{
		log_ofs << "start to merge incomplete alignment files ..." << endl;
		string cat_cmd_tmpAlignIncompletePair = "cat";
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			cat_cmd_tmpAlignIncompletePair = cat_cmd_tmpAlignIncompletePair + " " + tmpAlignIncompletePair + "." + int_to_str(tmp);
		}
		cat_cmd_tmpAlignIncompletePair = cat_cmd_tmpAlignIncompletePair + " " + tmpAlignIncompletePair;
		cat_cmd_tmpAlignIncompletePair = cat_cmd_tmpAlignIncompletePair + " > " + tmpAlignIncompletePair + ".all";
		system(cat_cmd_tmpAlignIncompletePair.c_str());
		log_ofs << "finish merging incomplete alignment files ..." << endl;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Sam 2 Junc   ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc starts ......" << endl << endl; 

	// new alignment2junc results:
	if(DoSam2JuncBool)
	{
		string juncfile = tmpIntermediateJunctionFile;
		log_ofs << "... initiate align2juncInfo ..." << endl;
		AlignmentToJunc_Info* align2juncInfo = new AlignmentToJunc_Info();
		int chromNum = indexInfo->returnChromNum();
		align2juncInfo->initiateAlignmentToJuncInfo(chromNum);
		log_ofs << "start to insert SJs into SJmap" << endl;
		vector<string> tmpAlignmentFileVec;
		tmpAlignmentFileVec.push_back(tmpAlignCompleteRead);
		tmpAlignmentFileVec.push_back(tmpAlignIncompletePair_SAM);
		tmpAlignmentFileVec.push_back(OutputSamFile_oneEndMapped);
		align2juncInfo->insertJuncFromAlignmentFileVec(tmpAlignmentFileVec, indexInfo);
		log_ofs << "finish inserting SJs into SJ map" << endl;
		align2juncInfo->outputSJmapVec(juncfile, indexInfo);
		delete align2juncInfo;
		int maxReadBaseNumInPathStructure = 30;
		alignInferJunctionHashInfo->insertJuncFromAlignmentFileVec(
			tmpAlignmentFileVec, indexInfo, maxReadBaseNumInPathStructure);
		string juncfile_alignInferHash = juncfile + ".alignInferHash";
		alignInferJunctionHashInfo->outputAlignInferInfoHashInfo(
			indexInfo, juncfile_alignInferHash);		
	}



	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... Sam 2 Junc ends ......" << endl << endl; 
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
	
	string OutputSamFile_fixHeadTail_complete_pair_alignInfo = outputDirStr + "/phase2_output/fixHeadTail_complete_pair.sam_alignInfo";
	ofstream OutputSamFile_fixHeadTail_complete_pair_alignInfo_ofs(OutputSamFile_fixHeadTail_complete_pair_alignInfo.c_str());

	string OutputSamFile_fixHeadTail_incomplete_pair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_ofs(OutputSamFile_fixHeadTail_incomplete_pair.c_str());
	
	string OutputSamFile_fixHeadTail_incomplete_pair_alignInfo = outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam_alignInfo";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_alignInfo_ofs(OutputSamFile_fixHeadTail_incomplete_pair_alignInfo.c_str());	

	string OutputSamFile_fixHeadTail_complete_unpair = outputDirStr + "/phase2_output/fixHeadTail_complete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_unpair_ofs(OutputSamFile_fixHeadTail_complete_unpair.c_str());

	string OutputSamFile_fixHeadTail_incomplete_unpair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_unpair_ofs(OutputSamFile_fixHeadTail_incomplete_unpair.c_str());	

	string OutputSamFile_fixHeadTail_pair_lowScore = outputDirStr + "/phase2_output/fixHeadTail_pair_lowScore.sam";
	ofstream OutputSamFile_fixHeadTail_pair_lowScore_ofs(OutputSamFile_fixHeadTail_pair_lowScore.c_str());	

	// start to read splice junction
	int junctionNum = 0;

	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	
	string InputSpliceJunction = tmpIntermediateJunctionFile;

	////////  take annotation as reference for remapping //////////////
	//InputSpliceJunction = "/data/homes/lxauky/GroundTruth2sam/resutls/sim1_test1/simulated_reads_test1_twoEnds.sam.noRandom.junc";
	settings_log_ofs << "InputSpliceJunction: " << InputSpliceJunction << endl; 
	log_ofs << "InputSpliceJunction: " << InputSpliceJunction << endl; 
	cout << "InputSpliceJunction: " << InputSpliceJunction << endl;
	///////////////////////////////////////////////////////////////////

	//if(annotation_provided_bool && Do_annotation_only_bool)
	//	InputSpliceJunction = annotation_file_path;	
	
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
			// loading SJs generated from aligned reads.
			FILE *fp_spliceJunction = fopen(InputSpliceJunction.c_str(), "r");
			fgets(entry, sizeof(entry), fp_spliceJunction);
			while(!feof(fp_spliceJunction))
			{
				fgets(entry, sizeof(entry), fp_spliceJunction);
				if(feof(fp_spliceJunction))
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
			fclose(fp_spliceJunction);
		}
		log_ofs << "after inserting SJs generated from alignments, junctionNum = " << junctionNum << endl;
		log_ofs << "start to insert SJs generated from annotation (if provided) " << endl;
		settings_log_ofs << "after inserting SJs generated from alignments, junctionNum = " << junctionNum << endl;
		settings_log_ofs << "start to insert SJs generated from annotation (if provided) " << endl;
		cout << "after inserting SJs generated from alignments, junctionNum = " << junctionNum << endl;
		cout << "start to insert SJs generated from annotation (if provided) " << endl;
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

		log_ofs << "After inserting SJs generated from alignments and annotation, junctionNum = " << junctionNum << endl;
		log_ofs << "finish building spliceJunction Hash" << endl;		
		log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		settings_log_ofs << "After inserting SJs generated from alignments and annotation, junctionNum = " << junctionNum << endl;
		settings_log_ofs << "finish building spliceJunction Hash" << endl;		
		settings_log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		cout << "After inserting SJs generated from alignments and annotation, junctionNum = " << junctionNum << endl;
		cout << "finish building spliceJunction Hash" << endl;		
		cout << "start doing remapping on unfixed head/tail alignments" << endl;
		string headTailSoftClippingFile = tmpAlignIncompletePair;// + ".all";

		//if(outputDirectlyBool_Phase1Only)
		//	headTailSoftClippingFile += ".all";
		AlignInfoInput_Array_Queue* fixHeadTailAlignInfoInputQueue = new AlignInfoInput_Array_Queue();
		Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue = new Result_FixHeadTail_Array_Queue();
		if(outputDirectlyBool_Phase1Only)
			headTailSoftClippingFile += ".all";
		ifstream inputUnfixedHeadTailRecord_ifs(headTailSoftClippingFile.c_str());
		//int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;
		bool endOfFile_bool = false;
		bool endOfProcessing_bool = false;

		int tmpInputReadNumInBatchArray_fixHeadTail = ReadNumInReadArray_FixHeadTail;
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
				process_stage_fixHeadTail(
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
					SE_or_PE_bool
					);
			}
		}
		fixHeadTailAlignInfoInputQueue->free();
		delete fixHeadTailAlignInfoInputQueue;
		fixHeadTailResultQueue->free();
		delete fixHeadTailResultQueue;
		inputUnfixedHeadTailRecord_ifs.close();
	}
	delete SJ;

	OutputSamFile_fixHeadTail_complete_pair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_ofs.close();
	OutputSamFile_fixHeadTail_complete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_pair_lowScore_ofs.close();
	
	nowtime = time(NULL);
	local = localtime(&nowtime);

	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	runtime_log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  

	statsInfo->getPhase1Stats();
	statsInfo->getFixUnpairedStats();
	statsInfo->getFixHeadTailStats();
	//statsInfo->outputAllStats(log_ofs, readTotalNum);
	statsInfo->outputAllStats(stats_ofs, Do_Phase1_Only);

	statsInfo->outputFinalStats(stats_ofs, Do_Phase1_Only);//, readTotalNum);	

	nowtime = time(NULL);
	local = localtime(&nowtime);

	cout << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  
	runtime_log_ofs << endl << "[" << asctime(local) << "... start to prepare for final output files ......" << endl << endl ;  

	string finalOutputSam = outputDirStr + "/output.sam";
	if(Do_Phase1_Only)
	{
		string cat_cmd = "cat "
			+ tmpHeadSectionInfo
			+ " " + tmpAlignCompleteRead  
			//+ " " + tmpAlignOneEndUnmapped
			+ " " + tmpAlignIncompletePair_SAM 
			+ " " + tmpAlignOneEndUnmapped
			+ " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile
			+ " " + tmpAlignBothEndsUnmapped_lowScore
			+ " " + tmpAlignBothEndsUnmapped
			+ " > " + finalOutputSam;
		system(cat_cmd.c_str()); 
	}
	else
	{
		string cat_cmd = "cat " 
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
		remove(tmpAlignCompleteRead_alignInfo.c_str());
		remove(OutputSamFile_oneEndMapped_alignInfo.c_str());
		remove(tmpAlignCompleteRead.c_str());
		remove(tmpHeadSectionInfo.c_str());

		remove(tmpAlignBothEndsUnmapped_lowScore.c_str());
		remove(OutputSamFile_fixHeadTail_pair_lowScore.c_str());
		remove(OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore.c_str());
	}

	annotation_ifs.close();
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